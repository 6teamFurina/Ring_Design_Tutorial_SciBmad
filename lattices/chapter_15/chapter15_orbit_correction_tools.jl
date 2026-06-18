using LinearAlgebra
using Random
using Statistics

function chapter15_find_file(parts...)
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

function chapter15_remove_bmad_comment(line)
    return strip(split(line, "!"; limit=2)[1])
end

function chapter15_join_bmad_continuations(raw_lines)
    joined = String[]
    buffer = ""

    for raw in raw_lines
        line = chapter15_remove_bmad_comment(raw)
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

function chapter15_parse_bmad_number(expr)
    # Tutorial lengths are simple numeric expressions, for example:
    # ((1.241+5.855+0.609)-2*3.41)/3
    try
        return Float64(Base.invokelatest(eval, Meta.parse(expr)))
    catch
        return 0.0
    end
end

function chapter15_parse_bmad_layout(path; root_line="RING", bpm_after=("B", "BH", "DB"))
    element_specs = Dict{String, NamedTuple}()
    line_defs = Dict{String, Vector{String}}()

    for line in chapter15_join_bmad_continuations(readlines(path))
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
        L = length_match === nothing ? 0.0 : chapter15_parse_bmad_number(length_match.captures[1])
        element_specs[name] = (; kind, L)
    end

    return chapter15_expand_bmad_layout(element_specs, line_defs, root_line; bpm_after)
end

function chapter15_expand_bmad_token(token, line_defs)
    token = strip(token)
    mult = match(r"^(\d+)\s*\*\s*(\w+)$", token)

    if mult !== nothing
        n = parse(Int, mult.captures[1])
        name = mult.captures[2]
        expanded = String[]
        for _ in 1:n
            append!(expanded, chapter15_expand_bmad_token(name, line_defs))
        end
        return expanded
    end

    if haskey(line_defs, token)
        expanded = String[]
        for child in line_defs[token]
            append!(expanded, chapter15_expand_bmad_token(child, line_defs))
        end
        return expanded
    end

    return [token]
end

function chapter15_expand_bmad_layout(element_specs, line_defs, root_line; bpm_after=("B", "BH", "DB"))
    tokens = chapter15_expand_bmad_token(root_line, line_defs)

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

function chapter15_sawtooth_orbit(s; circumference, amplitude=4.0e-4, n_rf=6, ripple=5.0e-5)
    phase = mod(n_rf * s / circumference, 1.0)
    ramp = amplitude * (2phase - 1)
    betatron_ripple = ripple * sin(2pi * 7.3 * s / circumference)
    return ramp + betatron_ripple
end

function chapter15_downstream_phase_advance(s_bpm, s_corrector; circumference, tune)
    ds = s_bpm >= s_corrector ? s_bpm - s_corrector : s_bpm - s_corrector + circumference
    return 2pi * tune * ds / circumference
end

function chapter15_beta_model(s; circumference, base=25.0, modulation=0.35, harmonic=12.0)
    return base * (1 + modulation * sin(2pi * harmonic * s / circumference))
end

function chapter15_alpha_model(s; circumference, base=25.0, modulation=0.35, harmonic=12.0)
    dbeta_ds = base * modulation * (2pi * harmonic / circumference) * cos(2pi * harmonic * s / circumference)
    return -0.5 * dbeta_ds
end

function chapter15_horizontal_response_matrix(bpm_s, ch_s; circumference, tune=54.23)
    response = zeros(length(bpm_s), length(ch_s))
    denom = 2sin(pi * tune)

    for i in eachindex(bpm_s)
        beta_i = chapter15_beta_model(bpm_s[i]; circumference)
        for j in eachindex(ch_s)
            beta_j = chapter15_beta_model(ch_s[j]; circumference, harmonic=12.0, modulation=0.25)
            dphi = chapter15_downstream_phase_advance(bpm_s[i], ch_s[j]; circumference, tune)
            response[i, j] = sqrt(beta_i * beta_j) / denom * cos(dphi - pi * tune)
        end
    end

    return response
end

function chapter15_phase_space_response(obs_s, kicker_s; circumference, tune=54.23, plane=:x)
    response = zeros(2 * length(obs_s), length(kicker_s))
    base = plane == :y ? 19.0 : 25.0
    modulation = plane == :y ? 0.28 : 0.35
    harmonic = plane == :y ? 10.0 : 12.0

    for i in eachindex(obs_s)
        beta_i = chapter15_beta_model(obs_s[i]; circumference, base, modulation, harmonic)
        alpha_i = chapter15_alpha_model(obs_s[i]; circumference, base, modulation, harmonic)
        for j in eachindex(kicker_s)
            kicker_s[j] <= obs_s[i] || continue

            beta_j = chapter15_beta_model(kicker_s[j]; circumference, base, modulation=0.25, harmonic)
            phase = 2pi * tune * (obs_s[i] - kicker_s[j]) / circumference
            response[2i - 1, j] = sqrt(beta_i * beta_j) * sin(phase)
            response[2i, j] = sqrt(beta_j / beta_i) * (cos(phase) - alpha_i * sin(phase))
        end
    end

    return response
end

function chapter15_local_bump_solution(kicker_s, ip_s, close_s, initial_ip_phase_space; circumference, tune=54.23, plane=:x)
    obs_s = [ip_s, close_s]
    response = chapter15_phase_space_response(obs_s, kicker_s; circumference, tune, plane)
    target_change = [-initial_ip_phase_space[1], -initial_ip_phase_space[2], 0.0, 0.0]
    kicks = response \ target_change
    achieved_change = response * kicks
    residual = achieved_change - target_change
    final_ip_phase_space = initial_ip_phase_space + achieved_change[1:2]
    closure_change = achieved_change[3:4]

    return (; kicks, response, achieved_change, residual, final_ip_phase_space, closure_change)
end

function chapter15_linear_interpolate(s_data, y_data, s0)
    idx = searchsortedlast(s_data, s0)
    idx = clamp(idx, 1, length(s_data) - 1)
    s1 = s_data[idx]
    s2 = s_data[idx + 1]
    y1 = y_data[idx]
    y2 = y_data[idx + 1]
    weight = (s0 - s1) / (s2 - s1)
    return (1 - weight) * y1 + weight * y2
end

function chapter15_local_slope(s_data, y_data, s0)
    idx = searchsortedlast(s_data, s0)
    idx = clamp(idx, 1, length(s_data) - 1)
    return (y_data[idx + 1] - y_data[idx]) / (s_data[idx + 1] - s_data[idx])
end

function chapter15_optimize_horizontal_correctors(bpm_s, ch_s, x_bpm; circumference, tune=54.23, regularization=1e-4)
    R = chapter15_horizontal_response_matrix(bpm_s, ch_s; circumference, tune)

    # Regularization keeps the fit from using unnecessarily large corrector kicks.
    A = [R; sqrt(regularization) * I(size(R, 2))]
    b = [-x_bpm; zeros(size(R, 2))]

    kicks = A \ b
    corrected_x = x_bpm + R * kicks

    return (; kicks, corrected_x, response=R)
end

function chapter15_orbit_stability_check(x_bpm; max_allowed=5.0e-3)
    max_abs = maximum(abs, x_bpm)
    rms = sqrt(mean(abs2, x_bpm))
    penalty = max(0.0, max_abs / max_allowed - 1.0)^2
    return (; stable=max_abs <= max_allowed, max_abs, rms, penalty)
end

function chapter15_single_pass_response_matrix(bpm_s, kicker_s; circumference, tune=54.23)
    response = zeros(length(bpm_s), length(kicker_s))

    for i in eachindex(bpm_s)
        beta_i = chapter15_beta_model(bpm_s[i]; circumference)
        for j in eachindex(kicker_s)
            if bpm_s[i] >= kicker_s[j]
                beta_j = chapter15_beta_model(kicker_s[j]; circumference, harmonic=12.0, modulation=0.25)
                phase = 2pi * tune * (bpm_s[i] - kicker_s[j]) / circumference
                response[i, j] = sqrt(beta_i * beta_j) * sin(phase)
            end
        end
    end

    return response
end

function chapter15_optimize_single_pass_horizontal_correctors(bpm_s, ch_s, x_bpm; circumference, tune=54.23, regularization=1e-3)
    R = chapter15_single_pass_response_matrix(bpm_s, ch_s; circumference, tune)
    A = [R; sqrt(regularization) * I(size(R, 2))]
    b = [-x_bpm; zeros(size(R, 2))]

    kicks = A \ b
    corrected_x = x_bpm + R * kicks

    return (; kicks, corrected_x, response=R)
end

function chapter15_format_float_vector(values; per_line=5)
    rows = String[]
    for start in 1:per_line:length(values)
        stop = min(start + per_line - 1, length(values))
        push!(rows, "    " * join((repr(Float64(v)) for v in values[start:stop]), ", "))
    end
    return "[\n" * join(rows, ",\n") * "\n]"
end

function chapter15_write_optimized_ring(path, optimized_ring)
    content = """
    # Generated from the Chapter 15.2 horizontal corrector optimization.
    # Include chapter15_sawtooth_ring_before.jl before this file.

    const CH14_SAWTOOTH_RING_AFTER = (
        layout = CH14_RING0_LAYOUT,
        ch_kicks = $(chapter15_format_float_vector(optimized_ring.ch_kicks)),
        x_bpm = $(chapter15_format_float_vector(optimized_ring.x_bpm)),
    )
    """

    write(path, content)
    return path
end

function chapter15_format_string_vector(values; per_line=4)
    rows = String[]
    for start in 1:per_line:length(values)
        stop = min(start + per_line - 1, length(values))
        push!(rows, "    " * join((repr(String(v)) for v in values[start:stop]), ", "))
    end
    return "[\n" * join(rows, ",\n") * "\n]"
end

function chapter15_write_local_bump_solution(path, local_bump)
    content = """
    # Generated from the Chapter 15.4 local IP bump exercise.
    # Include chapter15_sawtooth_ring_before.jl before using this with the ring layout.

    const CH14_LOCAL_BUMP_SOLUTION = (
        ip_s = $(repr(Float64(local_bump.ip_s))),
        close_s = $(repr(Float64(local_bump.close_s))),
        horizontal_kicker_names = $(chapter15_format_string_vector(local_bump.horizontal_kicker_names)),
        vertical_kicker_names = $(chapter15_format_string_vector(local_bump.vertical_kicker_names)),
        horizontal_kicker_s = $(chapter15_format_float_vector(local_bump.horizontal_kicker_s)),
        vertical_kicker_s = $(chapter15_format_float_vector(local_bump.vertical_kicker_s)),
        horizontal_kicks = $(chapter15_format_float_vector(local_bump.horizontal_kicks)),
        vertical_kicks = $(chapter15_format_float_vector(local_bump.vertical_kicks)),
        horizontal_initial_ip = $(chapter15_format_float_vector(local_bump.horizontal_initial_ip; per_line=2)),
        vertical_initial_ip = $(chapter15_format_float_vector(local_bump.vertical_initial_ip; per_line=2)),
        horizontal_final_ip = $(chapter15_format_float_vector(local_bump.horizontal_final_ip; per_line=2)),
        vertical_final_ip = $(chapter15_format_float_vector(local_bump.vertical_final_ip; per_line=2)),
        horizontal_closure_change = $(chapter15_format_float_vector(local_bump.horizontal_closure_change; per_line=2)),
        vertical_closure_change = $(chapter15_format_float_vector(local_bump.vertical_closure_change; per_line=2)),
        horizontal_residual = $(chapter15_format_float_vector(local_bump.horizontal_residual)),
        vertical_residual = $(chapter15_format_float_vector(local_bump.vertical_residual)),
    )
    """

    write(path, content)
    return path
end

function chapter15_truncated_gaussian(n; sigma=1.0e-5, cutoff=3.0, seed=14)
    values = zeros(n)
    rng = Random.MersenneTwister(seed)

    for i in 1:n
        while true
            candidate = sigma * randn(rng)
            if abs(candidate) <= cutoff * sigma
                values[i] = candidate
                break
            end
        end
    end

    return values
end

function chapter15_quadrupole_misalignment_orbit(bpm_s, quad_s, quad_k1l, quad_x_offsets; circumference, tune=54.23)
    response = chapter15_horizontal_response_matrix(bpm_s, quad_s; circumference, tune)
    equivalent_kicks = quad_k1l .* quad_x_offsets
    orbit = response * equivalent_kicks
    return (; orbit, equivalent_kicks, response)
end

function chapter15_quadrupole_single_pass_orbit(bpm_s, quad_s, quad_k1l, quad_x_offsets; circumference, tune=54.23)
    response = chapter15_single_pass_response_matrix(bpm_s, quad_s; circumference, tune)
    equivalent_kicks = quad_k1l .* quad_x_offsets
    orbit = response * equivalent_kicks
    return (; orbit, equivalent_kicks, response)
end
