using Printf

mutable struct TutorialBmadElement
    name::String
    kind::String
    attrs::Dict{String, Float64}
end

struct TutorialBmadOverlayTerm
    element::String
    attribute::String
    expression::String
end

struct TutorialBmadLattice
    elements::Dict{String, TutorialBmadElement}
    lines::Dict{String, Vector{String}}
    overlays::Dict{String, Vector{TutorialBmadOverlayTerm}}
    root_line::String
    p0c::Union{Nothing, Float64}
    particle::String
    radiation_damping_on::Bool
end

function tutorial_bmad_remove_comment(line)
    return strip(split(line, "!"; limit=2)[1])
end

function tutorial_bmad_join_continuations(raw_lines)
    joined = String[]
    buffer = ""

    for raw in raw_lines
        line = tutorial_bmad_remove_comment(raw)
        isempty(line) && continue

        buffer = isempty(buffer) ? line : buffer * " " * line
        if !endswith(buffer, ",")
            push!(joined, buffer)
            buffer = ""
        end
    end

    !isempty(buffer) && push!(joined, buffer)
    return joined
end

function tutorial_bmad_number(expr; variables=Dict{String, Float64}())
    clean = strip(expr)
    for (name, value) in variables
        clean = replace(clean, Regex("\\b$(name)\\b", "i") => "($(repr(value)))")
    end

    return Float64(Base.invokelatest(eval, Meta.parse(clean)))
end

function tutorial_bmad_attr_dict(text)
    attrs = Dict{String, Float64}()

    for part in split(text, ",")
        occursin("=", part) || continue
        key, value = split(part, "="; limit=2)
        key = lowercase(strip(key))

        try
            attrs[key] = tutorial_bmad_number(value)
        catch
            # Non-numeric attributes such as particle names are handled separately.
        end
    end

    return attrs
end

function tutorial_bmad_parse_overlay_terms(text)
    terms = TutorialBmadOverlayTerm[]
    body_match = match(r"\{([^}]*)\}", text)
    body_match === nothing && return terms

    for term in split(body_match.captures[1], ",")
        m = match(r"^\s*(\w+)\[(\w+)\]\s*:\s*(.+?)\s*$"i, term)
        m === nothing && continue
        push!(terms, TutorialBmadOverlayTerm(m.captures[1], lowercase(m.captures[2]), strip(m.captures[3])))
    end

    return terms
end

function tutorial_bmad_parse(path; root_line="RING")
    elements = Dict{String, TutorialBmadElement}()
    lines = Dict{String, Vector{String}}()
    overlays = Dict{String, Vector{TutorialBmadOverlayTerm}}()
    p0c = nothing
    particle = "electron"
    radiation_damping_on = false
    selected_line = root_line

    for line in tutorial_bmad_join_continuations(readlines(path))
        parameter_match = match(r"^parameter\[(\w+)\]\s*=\s*(.+)$"i, line)
        if parameter_match !== nothing
            key = lowercase(parameter_match.captures[1])
            value = strip(parameter_match.captures[2])
            if key == "p0c"
                p0c = tutorial_bmad_number(value)
            elseif key == "particle"
                particle = lowercase(value)
            end
            continue
        end

        radiation_match = match(r"^bmad_com\[radiation_damping_on\]\s*=\s*([TF])$"i, line)
        if radiation_match !== nothing
            radiation_damping_on = uppercase(radiation_match.captures[1]) == "T"
            continue
        end

        use_match = match(r"^use\s*,\s*(\w+)$"i, line)
        if use_match !== nothing
            selected_line = use_match.captures[1]
            continue
        end

        assignment_match = match(r"^(\w+)\[(\w+)\]\s*=\s*(.+)$"i, line)
        if assignment_match !== nothing
            name = assignment_match.captures[1]
            attr = lowercase(assignment_match.captures[2])
            value = tutorial_bmad_number(assignment_match.captures[3])

            if haskey(elements, name)
                elements[name].attrs[attr] = value
            elseif haskey(overlays, name)
                for term in overlays[name]
                    haskey(elements, term.element) || continue
                    elements[term.element].attrs[term.attribute] =
                        tutorial_bmad_number(term.expression; variables=Dict(attr => value))
                end
            end
            continue
        end

        line_match = match(r"^(\w+):\s*line\s*=\s*\((.*)\)\s*$"i, line)
        if line_match !== nothing
            name = line_match.captures[1]
            tokens = strip.(split(line_match.captures[2], ","))
            lines[name] = filter(!isempty, tokens)
            continue
        end

        overlay_match = match(r"^(\w+):\s*overlay\s*=\s*(.+)$"i, line)
        if overlay_match !== nothing
            overlays[overlay_match.captures[1]] = tutorial_bmad_parse_overlay_terms(overlay_match.captures[2])
            continue
        end

        element_match = match(r"^(\w+):\s*([A-Za-z]+)\b\s*,?\s*(.*)$"i, line)
        if element_match !== nothing
            name = element_match.captures[1]
            kind = lowercase(element_match.captures[2])
            attrs = tutorial_bmad_attr_dict(element_match.captures[3])
            elements[name] = TutorialBmadElement(name, kind, attrs)
        end
    end

    return TutorialBmadLattice(elements, lines, overlays, selected_line, p0c, particle, radiation_damping_on)
end

function tutorial_bmad_expand_token(token, lines)
    token = strip(token)
    mult = match(r"^(\d+)\s*\*\s*(\w+)$", token)

    if mult !== nothing
        n = parse(Int, mult.captures[1])
        name = mult.captures[2]
        expanded = String[]
        for _ in 1:n
            append!(expanded, tutorial_bmad_expand_token(name, lines))
        end
        return expanded
    end

    if haskey(lines, token)
        expanded = String[]
        for child in lines[token]
            append!(expanded, tutorial_bmad_expand_token(child, lines))
        end
        return expanded
    end

    return [token]
end

function tutorial_bmad_expanded_names(lat::TutorialBmadLattice)
    return tutorial_bmad_expand_token(lat.root_line, lat.lines)
end

function tutorial_bmad_float(value)
    return @sprintf("%.17e", Float64(value))
end

function tutorial_bmad_get(attrs, key, default=0.0)
    return get(attrs, lowercase(key), default)
end

function tutorial_bmad_instance_name(base_name, kind, counts)
    idx = get!(counts, base_name, 0) + 1
    counts[base_name] = idx

    if lowercase(kind) == "hkicker" && base_name == "CH"
        return "CH_$(idx)"
    elseif lowercase(kind) == "vkicker" && base_name == "CV"
        return "CV_$(idx)"
    elseif lowercase(kind) == "marker" && base_name == "BPM"
        return "BPM_$(idx)"
    elseif idx == 1
        return base_name
    end

    return "$(base_name)__$(idx)"
end

function tutorial_bmad_constructor(spec::TutorialBmadElement, instance_name; circumference=nothing)
    attrs = spec.attrs
    name_arg = "name=$(repr(instance_name))"
    L = tutorial_bmad_get(attrs, "l", 0.0)
    L_arg = "L=$(tutorial_bmad_float(L))"

    if spec.kind == "drift"
        return "Drift($name_arg, $L_arg)"
    elseif spec.kind == "quadrupole"
        return "Quadrupole($name_arg, $L_arg, Kn1=$(tutorial_bmad_float(tutorial_bmad_get(attrs, "k1"))))"
    elseif spec.kind == "sextupole"
        return "Sextupole($name_arg, $L_arg, Kn2=$(tutorial_bmad_float(tutorial_bmad_get(attrs, "k2"))))"
    elseif spec.kind == "sbend"
        g = tutorial_bmad_get(attrs, "g")
        e1 = tutorial_bmad_get(attrs, "e1")
        e2 = tutorial_bmad_get(attrs, "e2")
        return "SBend($name_arg, $L_arg, g_ref=$(tutorial_bmad_float(g)), Kn0=$(tutorial_bmad_float(g)), e1=$(tutorial_bmad_float(e1)), e2=$(tutorial_bmad_float(e2)))"
    elseif spec.kind == "rfcavity"
        parts = ["RFCavity($name_arg", L_arg]
        haskey(attrs, "voltage") && push!(parts, "voltage=$(tutorial_bmad_float(attrs["voltage"]))")
        if haskey(attrs, "harmon") && circumference !== nothing && circumference > 0
            push!(parts, "rf_frequency=$(tutorial_bmad_float(attrs["harmon"] * 299792458.0 / circumference))")
        end
        return join(parts, ", ") * ")"
    elseif spec.kind == "hkicker"
        return "HKicker($name_arg, $L_arg)"
    elseif spec.kind == "vkicker"
        return "VKicker($name_arg, $L_arg)"
    elseif spec.kind == "kicker"
        return "Kicker($name_arg, $L_arg)"
    elseif spec.kind == "marker"
        return "Marker($name_arg)"
    end

    return "Marker($name_arg)"
end

function tutorial_bmad_write_scibmad(path, lat::TutorialBmadLattice; const_name="RING", source="")
    expanded = tutorial_bmad_expanded_names(lat)
    counts = Dict{String, Int}()
    constructors = String[]
    unknown = String[]
    circumference = sum(
        haskey(lat.elements, token) ? tutorial_bmad_get(lat.elements[token].attrs, "l", 0.0) : 0.0
        for token in expanded
    )

    for token in expanded
        if !haskey(lat.elements, token)
            push!(unknown, token)
            continue
        end

        spec = lat.elements[token]
        instance_name = tutorial_bmad_instance_name(token, spec.kind, counts)
        push!(constructors, "    " * tutorial_bmad_constructor(spec, instance_name; circumference))
    end

    species = lowercase(lat.particle) == "electron" ? "electron" : lat.particle
    beamline_kwargs = ["species_ref=Species($(repr(species)))"]
    lat.p0c !== nothing && push!(beamline_kwargs, "pc_ref=$(tutorial_bmad_float(lat.p0c))")

    lines = String[]
    push!(lines, "# Generated by lattices/common/bmad_to_scibmad.jl.")
    !isempty(source) && push!(lines, "# Source: $(source)")
    push!(lines, "using SciBmad")
    push!(lines, "")
    push!(lines, "const $(const_name)_METADATA = (")
    push!(lines, "    root_line = $(repr(lat.root_line)),")
    push!(lines, "    n_expanded_elements = $(length(expanded)),")
    push!(lines, "    n_written_elements = $(length(constructors)),")
    push!(lines, "    circumference = $(tutorial_bmad_float(circumference)),")
    push!(lines, "    unknown_tokens = $(repr(unknown)),")
    push!(lines, "    radiation_damping_on_in_bmad = $(lat.radiation_damping_on),")
    push!(lines, ")")
    push!(lines, "")
    push!(lines, "const $(const_name) = Beamline([")
    push!(lines, join(constructors, ",\n"))
    push!(lines, "]; $(join(beamline_kwargs, ", ")))")
    push!(lines, "")

    write(path, join(lines, "\n"))
    return (; path, expanded, unknown)
end

function tutorial_bmad_to_scibmad(input_path, output_path; const_name="RING", source="")
    lat = tutorial_bmad_parse(input_path)
    return tutorial_bmad_write_scibmad(output_path, lat; const_name, source)
end
