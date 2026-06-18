using SciBmad
using CairoMakie
using Printf

const tutorial_root = normpath(joinpath(@__DIR__, "..", ".."))

const CH13_SPECIES = Species("electron")
const CH13_E_REF = 10e6
const CH13_N_DRIFT = 10
const CH13_N_BEND = 40
const CH13_N_QUAD = 40

function thin_kick_map(v, q::Nothing, p)
    hkick, vkick = p
    return ((v[1], v[2] + hkick, v[3], v[4] + vkick, v[5], v[6]), q)
end

function thin_kick(name; hkick=0.0, vkick=0.0)
    marker = Marker(name=name)
    marker.transport_map = thin_kick_map
    marker.transport_map_params = (hkick, vkick)
    return marker
end

function chapter14_drift_slices()
    return [Drift(name=@sprintf("D_%02d", i), L=0.5 / CH13_N_DRIFT) for i in 1:CH13_N_DRIFT]
end

function chapter14_bend_slices()
    return [
        SBend(
            name=@sprintf("B_%02d", i),
            L=0.5 / CH13_N_BEND,
            e1=(i == 1 ? 0.1 : 0.0),
            g_ref=1.0,
            Kn0=1.001,
        )
        for i in 1:CH13_N_BEND
    ]
end

function chapter14_quad_slices()
    return [Quadrupole(name=@sprintf("Q_%02d", i), L=0.6 / CH13_N_QUAD, Kn1=0.23) for i in 1:CH13_N_QUAD]
end

function build_chapter14_lattice(; b_vkick=0.0, q_hkick=0.0, q_vkick=0.0)
    b_kick = thin_kick("B_KICK"; vkick=b_vkick)
    q_kick = thin_kick("Q_KICK"; hkick=q_hkick, vkick=q_vkick)

    return Beamline(
        vcat(chapter14_drift_slices(), [b_kick], chapter14_bend_slices(), [q_kick], chapter14_quad_slices());
        species_ref=CH13_SPECIES,
        E_ref=CH13_E_REF,
    )
end

mutable struct LatticeState
    design
    model
    base
end

function chapter14_state()
    design = build_chapter14_lattice()
    return LatticeState(design, deepcopy(design), deepcopy(design))
end

function chapter14_figure27_state()
    state = chapter14_state()
    change_kick!(state.model, "B_KICK"; dvkick=-0.0005)
    set_kick!(state.model, "Q_KICK"; hkick=0.001)
    state.base = deepcopy(state.model)
    set_kick!(state.model, "Q_KICK"; vkick=5e-4)
    return state
end

function element_by_name(beamline, name)
    for element in beamline.line
        string(element.name) == name && return element
    end
    error("Element named $(name) was not found.")
end

function set_kick!(beamline, name; hkick=nothing, vkick=nothing)
    element = element_by_name(beamline, name)
    old_hkick, old_vkick = element.transport_map_params
    new_hkick = isnothing(hkick) ? old_hkick : hkick
    new_vkick = isnothing(vkick) ? old_vkick : vkick
    element.transport_map_params = (new_hkick, new_vkick)
    return beamline
end

function change_kick!(beamline, name; dhkick=0.0, dvkick=0.0)
    element = element_by_name(beamline, name)
    old_hkick, old_vkick = element.transport_map_params
    element.transport_map_params = (old_hkick + dhkick, old_vkick + dvkick)
    return beamline
end

function reference_bunch(beamline, coordinates)
    coordinates = reshape(coordinates, 1, 6)

    if hasproperty(beamline, :p_over_q_ref)
        return Bunch(coordinates; species=beamline.species_ref, p_over_q_ref=beamline.p_over_q_ref)
    elseif hasproperty(beamline, :R_ref)
        return Bunch(coordinates; species=beamline.species_ref, R_ref=beamline.R_ref)
    else
        return Bunch(coordinates; species=beamline.species_ref)
    end
end

function track_orbit(beamline; v0=zeros(6))
    bunch = reference_bunch(beamline, copy(v0))
    history = zeros(length(beamline.line) + 1, 6)
    history[1, :] .= bunch.coords.v[1, :]

    s = zeros(length(beamline.line) + 1)
    for (i, element) in enumerate(beamline.line)
        track!(bunch, element)
        history[i + 1, :] .= bunch.coords.v[1, :]
        s[i + 1] = s[i] + element.L
    end

    return (; s, x=history[:, 1], px=history[:, 2], y=history[:, 3], py=history[:, 4], history)
end

function orbit_difference(a, b)
    return (; s=a.s, x=a.x .- b.x, px=a.px .- b.px, y=a.y .- b.y, py=a.py .- b.py)
end

function readable_ylims!(ax, orbit; minimum_span_mm=0.08)
    values = 1e3 .* vcat(orbit.x, orbit.y)
    lo, hi = extrema(values)
    center = 0.5 * (lo + hi)
    span = max(hi - lo, minimum_span_mm)
    ylims!(ax, center - 0.55 * span, center + 0.55 * span)
end

function add_orbit_panel!(fig, row, title, orbit; color_x=:royalblue3, color_y=:darkorange2)
    ax = Axis(fig[row, 1], xlabel="s [m]", ylabel="orbit [mm]", title=title)
    lines!(ax, orbit.s, 1e3 .* orbit.x; color=color_x, linewidth=2, label="x")
    scatter!(ax, orbit.s, 1e3 .* orbit.x; color=color_x, markersize=8)
    lines!(ax, orbit.s, 1e3 .* orbit.y; color=color_y, linewidth=2, label="y")
    scatter!(ax, orbit.s, 1e3 .* orbit.y; color=color_y, markersize=8)
    xlims!(ax, -0.02, orbit.s[end] + 0.02)
    readable_ylims!(ax, orbit)
    axislegend(ax; position=:lt)
    return ax
end
