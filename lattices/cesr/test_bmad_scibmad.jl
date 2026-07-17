#!/usr/bin/env julia

"""
Element-by-element CESR comparison between a Tao/Bmad reference export and
SciBmad.  Bmad matrices are evaluated on the Bmad reference trajectory.  Each
SciBmad element is independently linearized about the same entrance point with
GTPSA, so a discrepancy is assigned to the element that creates it instead of
being hidden inside a one-turn product.

Run from the repository root with:

    julia --project=. lattices/cesr/test_bmad_scibmad.jl

Useful options:

    --parse-only          validate the Bmad archive without loading SciBmad
    --max-elements=N      compare only the first N tracking elements
    --strict              fail if the numerical tolerances are exceeded
    --reference=PATH      use another .tar.gz archive or raw fallback text file
    --csv=PATH            choose the per-element CSV output path
"""

using LinearAlgebra
using Printf
using SciBmad
using GTPSA
import Beamlines

const HERE = @__DIR__
const DEFAULT_REFERENCE = joinpath(HERE, "bmad_reference_output.tar.gz")
const DEFAULT_CSV = joinpath(HERE, "bmad_scibmad_comparison.csv")
const NUMBER_RE = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"

include(joinpath(HERE, "cesr.jl"))

struct BmadElementMap
    index::Int
    name::String
    key::String
    s_start::Float64
    s_end::Float64
    R::Matrix{Float64}
    vec0::Vector{Float64}
    orbit_out::Vector{Float64}
    symplectic_error::Float64
end

parse_number(s) = parse(Float64, replace(String(s), 'D' => 'E', 'd' => 'e'))
numbers(s) = parse_number.(getproperty.(collect(eachmatch(NUMBER_RE, s)), :match))

function read_reference(path::AbstractString)
    isfile(path) || error("Bmad reference file does not exist: $path")
    if endswith(lowercase(path), ".tar.gz") || endswith(lowercase(path), ".tgz")
        member = "bmad_reference_output/tao_fallback.txt"
        try
            return read(pipeline(`tar -xOzf $path $member`), String)
        catch err
            error("Could not read $member from $path with tar: $(sprint(showerror, err))")
        end
    end
    return read(path, String)
end

function parse_element_block(block::AbstractString)
    idx_match = match(r"(?m)^\s*Element #\s*(\d+)\s*$", block)
    idx_match === nothing && return nothing
    index = parse(Int, idx_match.captures[1])

    name_match = match(r"(?m)^\s*Element Name:\s*(.*?)\s*$", block)
    key_match = match(r"(?m)^\s*Key:\s*(.*?)\s*$", block)
    s_match = match(r"(?m)^\s*S_start, S:\s*([^,]+),\s*([^\r\n]+)", block)
    name_match === nothing && error("Element $index has no name")
    key_match === nothing && error("Element $index has no key")
    s_match === nothing && error("Element $index has no longitudinal positions")

    lines = split(block, '\n')
    matrix_header = findfirst(line -> occursin("Transfer Matrix", line), lines)
    matrix_header === nothing && error("Element $index has no transfer matrix")
    rows = Vector{Vector{Float64}}()
    for line in @view lines[(matrix_header + 1):end]
        isempty(strip(line)) && continue
        vals = numbers(line)
        if occursin(':', line) && length(vals) >= 7
            push!(rows, vals[1:7])
            length(rows) == 6 && break
        end
    end
    length(rows) == 6 || error("Element $index has only $(length(rows)) matrix rows")
    R = reduce(vcat, permutedims.(getindex.(rows, Ref(1:6))))
    vec0 = [row[7] for row in rows]

    orbit = Float64[]
    for label in ("X", "Y", "Z")
        m = match(Regex("(?m)^\\s*" * label * ":\\s*(" * NUMBER_RE.pattern * ")\\s+(" * NUMBER_RE.pattern * ")"), block)
        m === nothing && error("Element $index has no $label tracking coordinates")
        # Tao labels these columns Position[mm] and Momentum[1E-3].
        push!(orbit, 1e-3 * parse_number(m.captures[1]))
        push!(orbit, 1e-3 * parse_number(m.captures[2]))
    end

    symp_match = match(r"Mat symplectic error:\s*([^\]]+)", block)
    symp = symp_match === nothing ? NaN : parse_number(symp_match.captures[1])

    return BmadElementMap(
        index,
        strip(name_match.captures[1]),
        strip(key_match.captures[1]),
        parse_number(s_match.captures[1]),
        parse_number(s_match.captures[2]),
        R,
        vec0,
        orbit,
        symp,
    )
end

function parse_bmad_reference(text::AbstractString)
    result = Dict{Int,BmadElementMap}()
    command_re = r"(?m)^[ \t]*Tao>\s*show ele\s+[^\r\n]+\r?\n"
    matches = collect(eachmatch(command_re, text))
    for (j, command) in enumerate(matches)
        first_byte = command.offset + ncodeunits(command.match)
        last_byte = j == length(matches) ? ncodeunits(text) : matches[j + 1].offset - 1
        first_byte > last_byte && continue
        element = parse_element_block(SubString(text, first_byte, last_byte))
        element === nothing || (result[element.index] = element)
    end
    return result
end

function validate_reference(reference::Dict{Int,BmadElementMap})
    expected = Set(0:869)
    actual = Set(keys(reference))
    missing = setdiff(expected, actual)
    isempty(missing) || error("Expected Bmad indices 0:869; missing=$(sort!(collect(missing)))")
    beginning = reference[0]
    maximum(abs, beginning.R - I(6)) <= 1e-14 || error("Bmad BEGINNING map is not identity")
    maximum(abs, beginning.vec0) <= 1e-14 || error("Bmad BEGINNING kick is not zero")
    return nothing
end

csv_quote(x) = '"' * replace(string(x), '"' => "\"\"") * '"'

function element_kind(element)
    hasproperty(element, :kind) && return string(getproperty(element, :kind))
    return string(nameof(typeof(element)))
end

function element_length(element)
    hasproperty(element, :L) || return 0.0
    return Float64(getproperty(element, :L))
end

function matrix_worst(A::AbstractMatrix)
    location = argmax(abs.(A))
    return maximum(abs, A), location[1], location[2]
end

function write_csv(path, rows)
    columns = (
        :index, :name, :bmad_key, :scibmad_kind, :bmad_length_m,
        :scibmad_length_m, :length_error_m, :local_R_max, :local_R_fro,
        :local_R_row, :local_R_col, :vec0_max, :orbit_out_max,
        :bmad_affine_residual, :cumulative_R_max,
    )
    open(path, "w") do io
        println(io, join(string.(columns), ','))
        for row in rows
            values = map(column -> getproperty(row, column), columns)
            println(io, join(csv_quote.(values), ','))
        end
    end
end

function compare(reference; max_elements::Union{Nothing,Int}=nothing, csv_path=DEFAULT_CSV)
    ring = load_cesr()
    lattice_report = (
        kicker_elements=count(element -> element.kind == "Kicker", ring.line),
        nonzero_kicker_elements=count(
            element -> element.kind == "Kicker" &&
                       (!iszero(element.Kn0L) || !iszero(element.Ks0L)),
            ring.line,
        ),
        wiggler_elements=count(element -> element.kind == "Wiggler", ring.line),
        solenoid_quadrupole_overlaps=count(
            element -> element.kind == "Solenoid" &&
                       !iszero(element.Ksol) && !iszero(element.Kn1),
            ring.line,
        ),
    )
    n = isnothing(max_elements) ? length(ring.line) : min(max_elements, length(ring.line))
    length(ring.line) == 869 || error("Expected 869 SciBmad elements, got $(length(ring.line))")

    cumulative_bmad = Matrix{Float64}(I, 6, 6)
    cumulative_scibmad = Matrix{Float64}(I, 6, 6)
    rows = NamedTuple[]
    sizehint!(rows, n)

    for i in 1:n
        bmad = reference[i]
        bmad_in = reference[i - 1].orbit_out
        element = Beamlines.deepcopy_no_beamline(ring.line[i])
        one_element_line = Beamlines.Beamline(
            [element];
            p_over_q_ref=ring.p_over_q_ref,
            species_ref=ring.species_ref,
        )
        map = WigglerModels.gtpsa_transport_map(
            one_element_line;
            v0=bmad_in,
            order=1,
        )
        R_scibmad = Matrix(GTPSA.jacobian(map.v))
        out_scibmad = GTPSA.scalar.(map.v)
        vec0_scibmad = out_scibmad - R_scibmad * bmad_in

        delta_R = R_scibmad - bmad.R
        local_R_max, worst_row, worst_col = matrix_worst(delta_R)
        cumulative_bmad = bmad.R * cumulative_bmad
        cumulative_scibmad = R_scibmad * cumulative_scibmad
        bmad_affine_out = bmad.R * bmad_in + bmad.vec0

        push!(rows, (
            index=i,
            name=bmad.name,
            bmad_key=bmad.key,
            scibmad_kind=element_kind(element),
            bmad_length_m=bmad.s_end - bmad.s_start,
            scibmad_length_m=element_length(element),
            length_error_m=abs(element_length(element) - (bmad.s_end - bmad.s_start)),
            local_R_max=local_R_max,
            local_R_fro=norm(delta_R),
            local_R_row=worst_row,
            local_R_col=worst_col,
            vec0_max=maximum(abs, vec0_scibmad - bmad.vec0),
            orbit_out_max=maximum(abs, out_scibmad - bmad.orbit_out),
            bmad_affine_residual=maximum(abs, bmad_affine_out - bmad.orbit_out),
            cumulative_R_max=maximum(abs, cumulative_scibmad - cumulative_bmad),
        ))
        (i == 1 || i % 100 == 0 || i == n) && @printf("Compared %d/%d elements\n", i, n)
    end

    write_csv(csv_path, rows)
    return rows, lattice_report
end

function print_summary(rows, reference, csv_path)
    isempty(rows) && return
    worst_local = sort(rows; by=row -> row.local_R_max, rev=true)[1:min(15, length(rows))]
    worst_orbit = argmax(getproperty.(rows, :orbit_out_max))
    worst_cumulative = argmax(getproperty.(rows, :cumulative_R_max))
    max_affine = maximum(getproperty.(rows, :bmad_affine_residual))
    max_length = maximum(getproperty.(rows, :length_error_m))

    println("\nBmad reference: 869 tracking elements plus BEGINNING")
    println("Compared:       $(length(rows)) SciBmad elements")
    println("CSV report:     $csv_path")
    @printf("Bmad print/affine consistency: %.3e\n", max_affine)
    @printf("Maximum length mismatch:       %.3e m\n", max_length)
    @printf("Worst exit-orbit mismatch:     %.3e at #%d %s\n",
        rows[worst_orbit].orbit_out_max, rows[worst_orbit].index, rows[worst_orbit].name)
    @printf("Worst cumulative R mismatch:   %.3e at #%d %s\n",
        rows[worst_cumulative].cumulative_R_max,
        rows[worst_cumulative].index,
        rows[worst_cumulative].name)
    println("\nLargest local matrix discrepancies:")
    println(" index  element                    Bmad key          SciBmad kind       max|dR|      entry")
    for row in worst_local
        @printf(" %5d  %-25s  %-16s  %-16s  %10.3e  R%d%d\n",
            row.index, first(row.name, min(length(row.name), 25)), row.bmad_key,
            row.scibmad_kind, row.local_R_max, row.local_R_row, row.local_R_col)
    end
end

function option_value(args, prefix, default)
    hit = findfirst(arg -> startswith(arg, prefix), args)
    return hit === nothing ? default : split(args[hit], '='; limit=2)[2]
end

function main(args=ARGS)
    reference_path = abspath(option_value(args, "--reference=", DEFAULT_REFERENCE))
    csv_path = abspath(option_value(args, "--csv=", DEFAULT_CSV))
    parse_only = "--parse-only" in args
    strict = "--strict" in args
    max_text = option_value(args, "--max-elements=", "")
    max_elements = isempty(max_text) ? nothing : parse(Int, max_text)

    reference = parse_bmad_reference(read_reference(reference_path))
    validate_reference(reference)
    println("Validated Bmad reference tracking indices 0:869 from $reference_path")
    extra = sort!(collect(setdiff(Set(keys(reference)), Set(0:869))))
    isempty(extra) || println("Reference also contains non-tracking lord indices: $extra")
    parse_only && return 0

    rows, lattice_report = compare(reference; max_elements, csv_path)
    print_summary(rows, reference, csv_path)
    println("Consolidated CESR lattice: $lattice_report")

    if strict
        # The matrix tolerance is limited by Tao's seven-decimal matrix printout.
        maximum(getproperty.(rows, :length_error_m)) <= 1.1e-6 || error("Element lengths differ")
        maximum(getproperty.(rows, :local_R_max)) <= 1e-6 || error("Local R tolerance exceeded")
        maximum(getproperty.(rows, :orbit_out_max)) <= 1e-7 || error("Orbit tolerance exceeded")
    end
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
