using LinearAlgebra: inv

# Phase trombone utilities using SciBmad's normalizing map `a`.
#
# This version assumes twiss(...; normalizing_map=true) stores a normalizing map
# at each element. The trombone matrix is built as A * R(dphi_1, dphi_2) * A^-1.

function load_ring(path; make_trombones=false, kwargs...)
    ring = include(path)
    if make_trombones
        make_trombones!(ring; kwargs...)
    end
    return ring
end

function Trombone(; dnux=0.0, dnuy=0.0, kwargs...)
    return LineElement(;
        kind="Trombone",
        L=0.0,
        transport_map_params=(dnux, dnuy),
        kwargs...,
    )
end

function _normalizing_map_property(x, name, default=nothing)
    hasproperty(x, name) || return default
    try
        return getproperty(x, name)
    catch
        return default
    end
end

function _normalizing_map_trombone_beamline_index(element, ring)
    for (idx, candidate) in enumerate(ring.line)
        candidate === element && return idx
    end

    element_index = _normalizing_map_property(element, :beamline_index)
    if element_index isa Integer && firstindex(ring.line) <= element_index <= lastindex(ring.line)
        return element_index
    end

    element_name = _normalizing_map_property(element, :name)
    if element_name !== nothing
        row = findfirst(
            candidate -> _normalizing_map_property(candidate, :name) == element_name,
            ring.line,
        )
        row !== nothing && return row
    end

    error("Could not find trombone element in ring.line.")
end

function _normalizing_map_trombone_ring_element(element, ring)
    return ring.line[_normalizing_map_trombone_beamline_index(element, ring)]
end

function _normalizing_map_trombone_twiss_row(element, ring, twiss_table)
    idx = _normalizing_map_trombone_beamline_index(element, ring)

    if hasproperty(twiss_table, :beamline_index)
        row = findfirst(==(idx), twiss_table.beamline_index)
        row !== nothing && return row
    end

    return idx
end

function _normalizing_map_twiss(ring; twiss_kwargs...)
    try
        return twiss(ring; normalizing_map=true, twiss_kwargs...)
    catch err
        if err isa MethodError
            error(
                "twiss(...; normalizing_map=true) is not supported by this SciBmad version. " *
                "Update SciBmad to a version with the normalizing_map keyword, or pass a " *
                "twiss_table that already contains an `a` column."
            )
        end
        rethrow()
    end
end

function _normalizing_map_from_table(twiss_table, row)
    if hasproperty(twiss_table, :a)
        return twiss_table.a[row]
    end

    error(
        "Twiss table has no `a` column. Compute twiss with normalizing_map=true."
    )
end

function _normalizing_map_linear_matrix(a)
    if a isa AbstractMatrix
        return Matrix(a)
    end

    try
        return Matrix(jacobian(a, SciBmad.NNF.VARS_CPARAM))
    catch
        try
            return Matrix(jacobian(a))
        catch err
            error("Could not convert the normalizing map `a` into a 6x6 matrix: $err")
        end
    end
end

function _normalizing_map_check_6x6(A)
    size(A) == (6, 6) || error("Normalizing-map trombone expects a 6x6 A matrix.")
    return nothing
end

function normalizing_map_trombone_rotation(dphi_1, dphi_2, sample; sin_sign=1)
    c1 = cos(dphi_1)
    s1 = sin_sign * sin(dphi_1)
    c2 = cos(dphi_2)
    s2 = sin_sign * sin(dphi_2)

    one_entry = one(sample) * one(c1)
    zero_entry = zero(one_entry)
    R = fill(zero_entry, 6, 6)

    for i in 1:6
        R[i, i] = one_entry
    end

    R[1, 1] = c1
    R[1, 2] = s1
    R[2, 1] = -s1
    R[2, 2] = c1

    R[3, 3] = c2
    R[3, 4] = s2
    R[4, 3] = -s2
    R[4, 4] = c2

    return R
end

function normalizing_map_trombone_matrix(A, dphi_1, dphi_2; sin_sign=1)
    _normalizing_map_check_6x6(A)
    R = normalizing_map_trombone_rotation(dphi_1, dphi_2, A[1, 1]; sin_sign=sin_sign)
    return A * R * inv(A)
end

function normalizing_map_trombone_matrix(A, Ainv, dphi_1, dphi_2; sin_sign=1)
    _normalizing_map_check_6x6(A)
    _normalizing_map_check_6x6(Ainv)
    R = normalizing_map_trombone_rotation(dphi_1, dphi_2, A[1, 1]; sin_sign=sin_sign)
    return A * R * Ainv
end

function _normalizing_map_trombone_apply_6x6(M, v)
    return (
        M[1, 1] * v[1] + M[1, 2] * v[2] + M[1, 3] * v[3] + M[1, 4] * v[4] + M[1, 5] * v[5] + M[1, 6] * v[6],
        M[2, 1] * v[1] + M[2, 2] * v[2] + M[2, 3] * v[3] + M[2, 4] * v[4] + M[2, 5] * v[5] + M[2, 6] * v[6],
        M[3, 1] * v[1] + M[3, 2] * v[2] + M[3, 3] * v[3] + M[3, 4] * v[4] + M[3, 5] * v[5] + M[3, 6] * v[6],
        M[4, 1] * v[1] + M[4, 2] * v[2] + M[4, 3] * v[3] + M[4, 4] * v[4] + M[4, 5] * v[5] + M[4, 6] * v[6],
        M[5, 1] * v[1] + M[5, 2] * v[2] + M[5, 3] * v[3] + M[5, 4] * v[4] + M[5, 5] * v[5] + M[5, 6] * v[6],
        M[6, 1] * v[1] + M[6, 2] * v[2] + M[6, 3] * v[3] + M[6, 4] * v[4] + M[6, 5] * v[5] + M[6, 6] * v[6],
    )
end

function normalizing_map_trombone_map(A; sin_sign=1)
    _normalizing_map_check_6x6(A)
    Ainv = inv(A)

    return function (v, q, params)
        dnux = params[1]
        dnuy = params[2]
        M = normalizing_map_trombone_matrix(A, Ainv, dnux, dnuy; sin_sign=sin_sign)
        return (_normalizing_map_trombone_apply_6x6(M, v), q)
    end
end

function normalizing_map_trombone_A(element, ring, twiss_table)
    row = _normalizing_map_trombone_twiss_row(element, ring, twiss_table)
    A = _normalizing_map_linear_matrix(_normalizing_map_from_table(twiss_table, row))
    _normalizing_map_check_6x6(A)
    return A
end

function normalizing_map_trombone_params(dnux, dnuy)
    return (dnux, dnuy)
end

function normalizing_map_trombone_params(element)
    params = _normalizing_map_property(element, :transport_map_params)
    params === nothing && error("Trombone element has no transport_map_params. Expected (dnux, dnuy).")
    length(params) == 2 || error("Trombone element expects exactly two transport_map_params: (dnux, dnuy).")
    return params
end

function mark_trombone!(element, dnux, dnuy)
    element.kind = "Trombone"
    element.transport_map_params = normalizing_map_trombone_params(dnux, dnuy)
    return element
end

function mark_trombone!(element, ring, dnux, dnuy)
    mark_trombone!(element, dnux, dnuy)

    idx = _normalizing_map_trombone_beamline_index(element, ring)
    ring_element = ring.line[idx]
    mark_trombone!(ring_element, dnux, dnuy)
    return ring_element
end

function normalizing_map_trombone_elements(ring)
    return filter(element -> _normalizing_map_property(element, :kind) == "Trombone", ring.line)
end

function _set_normalizing_map_trombone!(element, A, dnux, dnuy; sin_sign=1)
    element.transport_map = normalizing_map_trombone_map(A; sin_sign=sin_sign)
    element.transport_map_params = normalizing_map_trombone_params(dnux, dnuy)
    return element
end

function make_trombone!(element, ring, twiss_table, dnux, dnuy; sin_sign=1)
    ring_element = _normalizing_map_trombone_ring_element(element, ring)
    A = normalizing_map_trombone_A(ring_element, ring, twiss_table)
    _set_normalizing_map_trombone!(ring_element, A, dnux, dnuy; sin_sign=sin_sign)

    if ring_element !== element
        _set_normalizing_map_trombone!(element, A, dnux, dnuy; sin_sign=sin_sign)
    end

    return ring_element
end

function make_trombone!(element, ring, twiss_table; sin_sign=1)
    dnux, dnuy = normalizing_map_trombone_params(element)
    return make_trombone!(element, ring, twiss_table, dnux, dnuy; sin_sign=sin_sign)
end

function make_trombone!(element, ring, dnux, dnuy; sin_sign=1, twiss_kwargs...)
    tw = _normalizing_map_twiss(ring; twiss_kwargs...)
    return make_trombone!(element, ring, tw.table, dnux, dnuy; sin_sign=sin_sign)
end

function make_trombone!(element, ring; sin_sign=1, twiss_kwargs...)
    tw = _normalizing_map_twiss(ring; twiss_kwargs...)
    return make_trombone!(element, ring, tw.table; sin_sign=sin_sign)
end

function _make_trombone_from_spec!(ring, twiss_table, element; sin_sign=1)
    return make_trombone!(element, ring, twiss_table; sin_sign=sin_sign)
end

function _make_trombone_from_spec!(ring, twiss_table, spec::Tuple; sin_sign=1)
    length(spec) == 3 || error("Trombone spec tuple must be (element, dnux, dnuy).")
    element, dnux, dnuy = spec
    return make_trombone!(element, ring, twiss_table, dnux, dnuy; sin_sign=sin_sign)
end

function make_trombones!(ring, trombones; sin_sign=1, twiss_kwargs...)
    tw = _normalizing_map_twiss(ring; twiss_kwargs...)

    for trombone in trombones
        _make_trombone_from_spec!(ring, tw.table, trombone; sin_sign=sin_sign)
    end

    return ring
end

function make_trombones!(ring; sin_sign=1, twiss_kwargs...)
    return make_trombones!(ring, normalizing_map_trombone_elements(ring); sin_sign=sin_sign, twiss_kwargs...)
end
