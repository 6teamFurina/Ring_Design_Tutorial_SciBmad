using LinearAlgebra
using Statistics

function chapter14_find_file(parts...)
    candidates = [
        joinpath(pwd(), parts...),
        joinpath(pwd(), "Ring_Design_Tutorial_SciBmad", parts...),
        joinpath(dirname(pwd()), "Ring_Design_Tutorial_SciBmad", parts...),
    ]

    for candidate in candidates
        isfile(candidate) && return candidate
    end

    error("Could not find file: " * joinpath(parts...))
end

function chapter14_remove_bmad_comment(line)
    return strip(split(line, "!"; limit=2)[1])
end

function chapter14_join_bmad_continuations(raw_lines)
    joined = String[]
    buffer = ""

    for raw in raw_lines
        line = chapter14_remove_bmad_comment(raw)
        isempty(line) && continue

        buffer = isempty(buffer) ? line : buffer * " " * line

        # This covers the multi-line element definitions used in the tutorial files.
        if !endswith(buffer, ",")
            push!(joined, buffer)
            buffer = ""
        end
    end

    !isempty(buffer) && push!(joined, buffer)
    return joined
end

function chapter14_parse_bmad_number(expr)
    # Tutorial lengths are simple numeric expressions, for example:
    # ((1.241+5.855+0.609)-2*3.41)/3
    try
        return Float64(Base.invokelatest(eval, Meta.parse(expr)))
    catch
        return 0.0
    end
end

function chapter14_parse_bmad_layout(path; root_line="RING", bpm_after=("B", "BH", "DB"))
    element_specs = Dict{String, NamedTuple}()
    line_defs = Dict{String, Vector{String}}()

    for line in chapter14_join_bmad_continuations(readlines(path))
        line_match = match(r"^(\w+):\s*line\s*=\s*\((.*)\)\s*$"i, line)
        if line_match !== nothing
            name = line_match.captures[1]
            tokens = strip.(split(line_match.captures[2], ","))
            line_defs[name] = filter(!isempty, tokens)
            continue
        end

        occursin(r"overlay"i, line) && continue
        element_match = match(r"^(\w+):\s*([A-Za-z]+)", line)
        element_match === nothing && continue

        name = element_match.captures[1]
        kind = element_match.captures[2]
        length_match = match(r"\bL\s*=\s*([^,]+)"i, line)
        L = length_match === nothing ? 0.0 : chapter14_parse_bmad_number(length_match.captures[1])
        element_specs[name] = (; kind, L)
    end

    return chapter14_expand_bmad_layout(element_specs, line_defs, root_line; bpm_after)
end

function chapter14_expand_bmad_token(token, line_defs)
    token = strip(token)
    mult = match(r"^(\d+)\s*\*\s*(\w+)$", token)

    if mult !== nothing
        n = parse(Int, mult.captures[1])
        name = mult.captures[2]
        expanded = String[]
        for _ in 1:n
            append!(expanded, chapter14_expand_bmad_token(name, line_defs))
        end
        return expanded
    end

    if haskey(line_defs, token)
        expanded = String[]
        for child in line_defs[token]
            append!(expanded, chapter14_expand_bmad_token(child, line_defs))
        end
        return expanded
    end

    return [token]
end

function chapter14_expand_bmad_layout(element_specs, line_defs, root_line; bpm_after=("B", "BH", "DB"))
    tokens = chapter14_expand_bmad_token(root_line, line_defs)

    rows = NamedTuple[]
    ch_names = String[]
    ch_s = Float64[]
    cv_names = String[]
    cv_s = Float64[]
    bpm_names = String[]
    bpm_s = Float64[]

    s = 0.0
    ch_count = 0
    cv_count = 0
    bpm_count = 0

    for token in tokens
        spec = get(element_specs, token, (; kind="Unknown", L=0.0))
        s += spec.L

        if token == "CH"
            ch_count += 1
            push!(ch_names, "CH_$(ch_count)")
            push!(ch_s, s)
        elseif token == "CV"
            cv_count += 1
            push!(cv_names, "CV_$(cv_count)")
            push!(cv_s, s)
        end

        if token in bpm_after
            bpm_count += 1
            push!(bpm_names, "BPM_$(bpm_count)")
            push!(bpm_s, s)
        end

        push!(rows, (; name=token, kind=spec.kind, L=spec.L, s=s))
    end

    return (; rows, tokens, ch_names, ch_s, cv_names, cv_s, bpm_names, bpm_s, circumference=s)
end

function chapter14_sawtooth_orbit(s; circumference, amplitude=4.0e-4, n_rf=6, ripple=5.0e-5)
    phase = mod(n_rf * s / circumference, 1.0)
    ramp = amplitude * (2phase - 1)
    betatron_ripple = ripple * sin(2pi * 7.3 * s / circumference)
    return ramp + betatron_ripple
end

function chapter14_downstream_phase_advance(s_bpm, s_corrector; circumference, tune)
    ds = s_bpm >= s_corrector ? s_bpm - s_corrector : s_bpm - s_corrector + circumference
    return 2pi * tune * ds / circumference
end

function chapter14_beta_model(s; circumference, base=25.0, modulation=0.35, harmonic=12.0)
    return base * (1 + modulation * sin(2pi * harmonic * s / circumference))
end

function chapter14_horizontal_response_matrix(bpm_s, ch_s; circumference, tune=54.23)
    response = zeros(length(bpm_s), length(ch_s))
    denom = 2sin(pi * tune)

    for i in eachindex(bpm_s)
        beta_i = chapter14_beta_model(bpm_s[i]; circumference)
        for j in eachindex(ch_s)
            beta_j = chapter14_beta_model(ch_s[j]; circumference, harmonic=12.0, modulation=0.25)
            dphi = chapter14_downstream_phase_advance(bpm_s[i], ch_s[j]; circumference, tune)
            response[i, j] = sqrt(beta_i * beta_j) / denom * cos(dphi - pi * tune)
        end
    end

    return response
end

function chapter14_optimize_horizontal_correctors(bpm_s, ch_s, x_bpm; circumference, tune=54.23, regularization=1e-4)
    R = chapter14_horizontal_response_matrix(bpm_s, ch_s; circumference, tune)

    # Regularization keeps the fit from using unnecessarily large corrector kicks.
    A = [R; sqrt(regularization) * I(size(R, 2))]
    b = [-x_bpm; zeros(size(R, 2))]

    kicks = A \ b
    corrected_x = x_bpm + R * kicks

    return (; kicks, corrected_x, response=R)
end

function chapter14_format_float_vector(values; per_line=5)
    rows = String[]
    for start in 1:per_line:length(values)
        stop = min(start + per_line - 1, length(values))
        push!(rows, "    " * join((repr(Float64(v)) for v in values[start:stop]), ", "))
    end
    return "[\n" * join(rows, ",\n") * "\n]"
end

function chapter14_write_optimized_ring(path, optimized_ring)
    content = """
    # Generated from the Chapter 14.2 horizontal corrector optimization.
    # Include chapter14_sawtooth_ring_before.jl before this file.

    const CH14_SAWTOOTH_RING_AFTER = (
        layout = CH14_RING0_LAYOUT,
        ch_kicks = $(chapter14_format_float_vector(optimized_ring.ch_kicks)),
        x_bpm = $(chapter14_format_float_vector(optimized_ring.x_bpm)),
    )
    """

    write(path, content)
    return path
end
