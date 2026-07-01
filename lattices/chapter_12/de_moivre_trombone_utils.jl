# Phase trombone utilities using SciBmad's De Moivre normal-form matrices.
#
# This version does not use alpha/beta Twiss parameters. It builds the
# trombone directly from H_i = A I_i A^-1 and B_i = A J_i A^-1.

function _de_moivre_trombone_beamline_index(element, ring)
    for (idx, candidate) in enumerate(ring.line)
        candidate === element && return idx
    end

    error("Could not find trombone element in ring.line.")
end

function _de_moivre_trombone_twiss_row(element, ring, twiss_table)
    idx = _de_moivre_trombone_beamline_index(element, ring)

    if hasproperty(twiss_table, :beamline_index)
        row = findfirst(==(idx), twiss_table.beamline_index)
        row !== nothing && return row
    end

    return idx
end

function _de_moivre_trombone_check_modes(H, B)
    length(H) >= 2 || error("De Moivre trombone needs at least two normal modes in H.")
    length(B) >= 2 || error("De Moivre trombone needs at least two normal modes in B.")
    size(H[1]) == (6, 6) || error("De Moivre trombone expects 6x6 H matrices.")
    size(B[1]) == (6, 6) || error("De Moivre trombone expects 6x6 B matrices.")
    return nothing
end

function de_moivre_trombone_matrix(H, B, dphi_1, dphi_2; sin_sign=1)
    _de_moivre_trombone_check_modes(H, B)

    return one(H[1]) +
        (cos(dphi_1) - 1) * H[1] + sin_sign * sin(dphi_1) * B[1] +
        (cos(dphi_2) - 1) * H[2] + sin_sign * sin(dphi_2) * B[2]
end

function _de_moivre_trombone_apply_6x6(M, v)
    return (
        M[1, 1] * v[1] + M[1, 2] * v[2] + M[1, 3] * v[3] + M[1, 4] * v[4] + M[1, 5] * v[5] + M[1, 6] * v[6],
        M[2, 1] * v[1] + M[2, 2] * v[2] + M[2, 3] * v[3] + M[2, 4] * v[4] + M[2, 5] * v[5] + M[2, 6] * v[6],
        M[3, 1] * v[1] + M[3, 2] * v[2] + M[3, 3] * v[3] + M[3, 4] * v[4] + M[3, 5] * v[5] + M[3, 6] * v[6],
        M[4, 1] * v[1] + M[4, 2] * v[2] + M[4, 3] * v[3] + M[4, 4] * v[4] + M[4, 5] * v[5] + M[4, 6] * v[6],
        M[5, 1] * v[1] + M[5, 2] * v[2] + M[5, 3] * v[3] + M[5, 4] * v[4] + M[5, 5] * v[5] + M[5, 6] * v[6],
        M[6, 1] * v[1] + M[6, 2] * v[2] + M[6, 3] * v[3] + M[6, 4] * v[4] + M[6, 5] * v[5] + M[6, 6] * v[6],
    )
end

function de_moivre_trombone_map(v, q, params)
    dphi_1 = params[1]
    dphi_2 = params[2]
    H = params[3]
    B = params[4]
    sin_sign = params[5]

    M = de_moivre_trombone_matrix(H, B, dphi_1, dphi_2; sin_sign=sin_sign)
    return (_de_moivre_trombone_apply_6x6(M, v), q)
end

function de_moivre_trombone_params(element, ring, twiss_table, dphi_1, dphi_2; sin_sign=1)
    row = _de_moivre_trombone_twiss_row(element, ring, twiss_table)
    H = Tuple(twiss_table.H[row])
    B = Tuple(twiss_table.B[row])
    _de_moivre_trombone_check_modes(H, B)
    return (dphi_1, dphi_2, H, B, sin_sign)
end

function attach_de_moivre_trombone!(element, ring, twiss_table, dphi_1, dphi_2; sin_sign=1)
    element.transport_map = de_moivre_trombone_map
    element.transport_map_params = de_moivre_trombone_params(
        element,
        ring,
        twiss_table,
        dphi_1,
        dphi_2;
        sin_sign=sin_sign,
    )
    return element
end

function attach_de_moivre_trombone!(element, ring, dphi_1, dphi_2; sin_sign=1, twiss_kwargs...)
    tw = twiss(ring; de_moivre=true, twiss_kwargs...)
    return attach_de_moivre_trombone!(element, ring, tw.table, dphi_1, dphi_2; sin_sign=sin_sign)
end

function attach_de_moivre_trombones!(ring, trombones; sin_sign=1, twiss_kwargs...)
    tw = twiss(ring; de_moivre=true, twiss_kwargs...)

    for (element, dphi_1, dphi_2) in trombones
        attach_de_moivre_trombone!(element, ring, tw.table, dphi_1, dphi_2; sin_sign=sin_sign)
    end

    return ring
end

const attach_coupled_trombone! = attach_de_moivre_trombone!
const attach_coupled_trombones! = attach_de_moivre_trombones!
