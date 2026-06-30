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

function _trombone_twiss_row(element, ring, twiss_table)
    idx = _trombone_beamline_index(element, ring)

    if hasproperty(twiss_table, :beamline_index)
        row = findfirst(==(idx), twiss_table.beamline_index)
        row !== nothing && return row
    end

    return idx
end

function _trombone_apply_6x6(M, v)
    return (
        M[1, 1] * v[1] + M[1, 2] * v[2] + M[1, 3] * v[3] + M[1, 4] * v[4] + M[1, 5] * v[5] + M[1, 6] * v[6],
        M[2, 1] * v[1] + M[2, 2] * v[2] + M[2, 3] * v[3] + M[2, 4] * v[4] + M[2, 5] * v[5] + M[2, 6] * v[6],
        M[3, 1] * v[1] + M[3, 2] * v[2] + M[3, 3] * v[3] + M[3, 4] * v[4] + M[3, 5] * v[5] + M[3, 6] * v[6],
        M[4, 1] * v[1] + M[4, 2] * v[2] + M[4, 3] * v[3] + M[4, 4] * v[4] + M[4, 5] * v[5] + M[4, 6] * v[6],
        M[5, 1] * v[1] + M[5, 2] * v[2] + M[5, 3] * v[3] + M[5, 4] * v[4] + M[5, 5] * v[5] + M[5, 6] * v[6],
        M[6, 1] * v[1] + M[6, 2] * v[2] + M[6, 3] * v[3] + M[6, 4] * v[4] + M[6, 5] * v[5] + M[6, 6] * v[6],
    )
end

function _trombone_check_6x6_modes(H, B)
    length(H) >= 2 || error("De Moivre trombone needs at least two normal modes in H.")
    length(B) >= 2 || error("De Moivre trombone needs at least two normal modes in B.")
    size(H[1]) == (6, 6) || error("De Moivre trombone expects 6x6 H matrices.")
    size(B[1]) == (6, 6) || error("De Moivre trombone expects 6x6 B matrices.")
    return nothing
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
    idx = _trombone_twiss_row(element, ring, twiss_table)

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

function attach_trombone!(element, ring, dnu_x, dnu_y; twiss_kwargs...)
    tw = twiss(ring; twiss_kwargs...)
    return attach_trombone!(element, ring, tw.table, dnu_x, dnu_y)
end

function de_moivre_trombone_matrix(H, B, dnu_1, dnu_2; sin_sign=1)
    _trombone_check_6x6_modes(H, B)

    return one(H[1]) +
        (cos(dnu_1) - 1) * H[1] + sin_sign * sin(dnu_1) * B[1] +
        (cos(dnu_2) - 1) * H[2] + sin_sign * sin(dnu_2) * B[2]
end

function de_moivre_trombone_map(v, q, p)
    dnu_1 = p[1]
    dnu_2 = p[2]
    H = p[3]
    B = p[4]
    sin_sign = p[5]

    M = de_moivre_trombone_matrix(H, B, dnu_1, dnu_2; sin_sign=sin_sign)
    return (_trombone_apply_6x6(M, v), q)
end

function de_moivre_trombone_params(element, ring, twiss_table, dnu_1, dnu_2; sin_sign=1)
    row = _trombone_twiss_row(element, ring, twiss_table)
    H = Tuple(twiss_table.H[row])
    B = Tuple(twiss_table.B[row])
    _trombone_check_6x6_modes(H, B)
    return (dnu_1, dnu_2, H, B, sin_sign)
end

function attach_de_moivre_trombone!(element, ring, twiss_table, dnu_1, dnu_2; sin_sign=1)
    element.transport_map = de_moivre_trombone_map
    element.transport_map_params = de_moivre_trombone_params(
        element,
        ring,
        twiss_table,
        dnu_1,
        dnu_2;
        sin_sign=sin_sign,
    )
    return element
end

function attach_de_moivre_trombone!(element, ring, dnu_1, dnu_2; sin_sign=1, twiss_kwargs...)
    tw = twiss(ring; de_moivre=true, twiss_kwargs...)
    return attach_de_moivre_trombone!(element, ring, tw.table, dnu_1, dnu_2; sin_sign=sin_sign)
end

const attach_coupled_trombone! = attach_de_moivre_trombone!

function attach_de_moivre_trombones!(ring, trombones; sin_sign=1, twiss_kwargs...)
    tw = twiss(ring; de_moivre=true, twiss_kwargs...)

    for (element, dnu_1, dnu_2) in trombones
        attach_de_moivre_trombone!(element, ring, tw.table, dnu_1, dnu_2; sin_sign=sin_sign)
    end

    return ring
end

const attach_coupled_trombones! = attach_de_moivre_trombones!
