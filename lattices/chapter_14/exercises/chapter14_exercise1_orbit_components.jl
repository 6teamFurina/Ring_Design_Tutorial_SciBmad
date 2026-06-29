include(joinpath(@__DIR__, "..", "chapter14_common.jl"))

function add_mixed_middle_panel!(fig, row, title, model_orbit, design_orbit)
    ax = Axis(fig[row, 1], xlabel="s [m]", ylabel="orbit [mm]", title=title)

    lines!(ax, model_orbit.s, 1e3 .* model_orbit.x; color=:royalblue3, linewidth=2, label="model x")
    scatter!(ax, model_orbit.s, 1e3 .* model_orbit.x; color=:royalblue3, markersize=8)

    lines!(ax, design_orbit.s, 1e3 .* design_orbit.x; color=:darkorange2, linewidth=2, label="design x")
    scatter!(ax, design_orbit.s, 1e3 .* design_orbit.x; color=:darkorange2, markersize=8)

    combined = (; s=model_orbit.s, x=model_orbit.x, y=design_orbit.x)
    readable_ylims!(ax, combined)
    xlims!(ax, -0.02, model_orbit.s[end] + 0.02)
    axislegend(ax; position=:lt)
    return ax
end

function plot_exercise1_orbit_components(state)
    design_orbit = track_orbit(state.design)
    model_orbit = track_orbit(state.model)
    base_orbit = track_orbit(state.base)

    fig = Figure(size=(760, 780))
    ax1 = add_orbit_panel!(fig, 1, "Orbit [model - design]", orbit_difference(model_orbit, design_orbit))
    ax2 = add_mixed_middle_panel!(fig, 2, "Orbit [model x and design x]", model_orbit, design_orbit)
    ax3 = add_orbit_panel!(fig, 3, "Orbit [model - base]", orbit_difference(model_orbit, base_orbit))
    linkxaxes!(ax1, ax2, ax3)
    return fig
end

state = chapter14_figure27_state()
exercise1_figure = plot_exercise1_orbit_components(state)
display(exercise1_figure)
