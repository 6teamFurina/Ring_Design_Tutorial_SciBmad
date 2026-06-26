using SciBmad

function chapter15_thin_kick_map(v, q::Nothing, p)
    hkick, vkick = p
    return ((v[1], v[2] + hkick, v[3], v[4] + vkick, v[5], v[6]), q)
end

function chapter15_layout_event_rows(layout; ch_kicks=nothing, cv_kicks=nothing, quad_k1l=nothing, quad_scale=1.0)
    events = NamedTuple[]

    for (i, s) in pairs(layout.bpm_s)
        push!(events, (; s, priority=30, name=layout.bpm_names[i], kind=:bpm, strength=0.0))
    end

    for (i, s) in pairs(layout.ch_s)
        kick = isnothing(ch_kicks) ? 0.0 : ch_kicks[i]
        push!(events, (; s, priority=20, name=layout.ch_names[i], kind=:ch, strength=kick))
    end

    for (i, s) in pairs(layout.cv_s)
        kick = isnothing(cv_kicks) ? 0.0 : cv_kicks[i]
        push!(events, (; s, priority=21, name=layout.cv_names[i], kind=:cv, strength=kick))
    end

    strengths = quad_scale .* (isnothing(quad_k1l) ? layout.quad_k1l : quad_k1l)
    for (i, s) in pairs(layout.quad_s)
        push!(events, (; s, priority=10, name=layout.quad_names[i], kind=:quad, strength=strengths[i]))
    end

    return sort(events, by = row -> (row.s, row.priority, row.name))
end

function chapter15_event_element(event)
    if event.kind == :bpm
        return Marker(name=event.name)
    elseif event.kind == :ch
        ele = Kicker(name=event.name)
        ele.transport_map = chapter15_thin_kick_map
        ele.transport_map_params = (event.strength, 0.0)
        return ele
    elseif event.kind == :cv
        ele = Kicker(name=event.name)
        ele.transport_map = chapter15_thin_kick_map
        ele.transport_map_params = (0.0, event.strength)
        return ele
    elseif event.kind == :quad
        return Multipole(name=event.name, Kn1L=event.strength)
    end

    error("Unsupported Chapter 15 layout event kind: $(event.kind)")
end

function chapter15_layout_to_beamline(
    layout;
    ch_kicks=nothing,
    cv_kicks=nothing,
    quad_k1l=nothing,
    quad_scale=1.0,
    species_ref=Species("electron"),
    E_ref=18e9,
    drift_prefix="CH15_D",
    min_drift=1e-12,
)
    events = chapter15_layout_event_rows(layout; ch_kicks, cv_kicks, quad_k1l, quad_scale)
    elements = LineElement[]
    s_now = 0.0
    drift_count = 0

    for event in events
        ds = event.s - s_now
        if ds > min_drift
            drift_count += 1
            push!(elements, Drift(name="$(drift_prefix)_$(drift_count)", L=ds))
            s_now = event.s
        elseif ds < -min_drift
            error("Layout events are not sorted in increasing s.")
        end

        push!(elements, chapter15_event_element(event))
    end

    tail = layout.circumference - s_now
    if tail > min_drift
        drift_count += 1
        push!(elements, Drift(name="$(drift_prefix)_$(drift_count)", L=tail))
    end

    return Beamline(elements; species_ref, E_ref)
end

function chapter15_beamline_twiss_table(ring; kwargs...)
    return twiss(ring; kwargs...).table
end

function chapter15_find_stable_layout_twiss(
    layout;
    candidate_quad_scales=vcat(10.0 .^ range(-4, 0, length=17), -10.0 .^ range(-4, 0, length=17)),
    kwargs...
)
    last_error = nothing

    for quad_scale in candidate_quad_scales
        ring = chapter15_layout_to_beamline(layout; quad_scale, kwargs...)
        try
            table = chapter15_beamline_twiss_table(ring)
            return (; quad_scale, ring, table)
        catch err
            last_error = err
        end
    end

    error("No stable layout Twiss found for the supplied candidate scales. Last error: $(last_error)")
end
