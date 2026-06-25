# Utilities for implementing a phase trombone with a custom SciBmad transport map.
#
# A trombone changes the betatron phase advance while preserving the local
# Twiss coordinates used to build the map. The phase inputs dnu_x and dnu_y are
# phase advances in radians, not tune units.

function _trombone_pz_coeff(x)
    try
        return x[6]
    catch
        return zero(x)
    end
end

function _trombone_beamline_index(element, ring)
    for (idx, candidate) in enumerate(ring.line)
        candidate === element && return idx
    end

    error("Could not find trombone element in ring.line.")
end

function trombone_map(v, q::Nothing, p)
    dnu_x = p[1]
    dnu_y = p[2]
    beta_x = p[3]
    alpha_x = p[4]
    beta_y = p[5]
    alpha_y = p[6]
    eta_x = p[7]
    etap_x = p[8]
    eta_y = p[9]
    etap_y = p[10]

    A11 = cos(dnu_x) + alpha_x * sin(dnu_x)
    A12 = beta_x * sin(dnu_x)
    A21 = -(1 + alpha_x^2) / beta_x * sin(dnu_x)
    A22 = cos(dnu_x) - alpha_x * sin(dnu_x)

    Dx = (1 - A11) * eta_x - A12 * etap_x
    Dpx = -A21 * eta_x + (1 - A22) * etap_x

    Tx = -A11 * Dpx + A21 * Dx
    Tpx = -A12 * Dpx + A22 * Dx

    x = A11 * v[1] + A12 * v[2] + Dx * v[6]
    px = A21 * v[1] + A22 * v[2] + Dpx * v[6]

    B11 = cos(dnu_y) + alpha_y * sin(dnu_y)
    B12 = beta_y * sin(dnu_y)
    B21 = -(1 + alpha_y^2) / beta_y * sin(dnu_y)
    B22 = cos(dnu_y) - alpha_y * sin(dnu_y)

    Dy = (1 - B11) * eta_y - B12 * etap_y
    Dpy = -B21 * eta_y + (1 - B22) * etap_y

    Ty = -B11 * Dpy + B21 * Dy
    Tpy = -B12 * Dpy + B22 * Dy

    y = B11 * v[3] + B12 * v[4] + Dy * v[6]
    py = B21 * v[3] + B22 * v[4] + Dpy * v[6]

    # Keep the map symplectic to first order by including the path-length terms.
    z = Tx * v[1] + Tpx * v[2] + Ty * v[3] + Tpy * v[4] + v[5]
    pz = v[6]

    return ((x, px, y, py, z, pz), q)
end

function trombone_params(element, ring, twiss_table, dnu_x, dnu_y)
    idx = _trombone_beamline_index(element, ring)

    return (
        dnu_x,
        dnu_y,
        twiss_table.beta_1[idx],
        twiss_table.alpha_1[idx],
        twiss_table.beta_2[idx],
        twiss_table.alpha_2[idx],
        _trombone_pz_coeff(twiss_table.orbit_x[idx]),
        _trombone_pz_coeff(twiss_table.orbit_px[idx]),
        _trombone_pz_coeff(twiss_table.orbit_y[idx]),
        _trombone_pz_coeff(twiss_table.orbit_py[idx]),
    )
end

function attach_trombone!(element, ring, twiss_table, dnu_x, dnu_y)
    element.transport_map = trombone_map
    element.transport_map_params = trombone_params(element, ring, twiss_table, dnu_x, dnu_y)
    return element
end
