include(joinpath(@__DIR__, "..", "chapter14_common.jl"))

using GTPSA
using LinearAlgebra

const TUNE_DESCRIPTOR = Descriptor(6, 2, 2, 1)
const DK_TUNE = params(TUNE_DESCRIPTOR)
const ZERO6 = [0, 0, 0, 0, 0, 0]

function tps_const(x)
    try
        return x[ZERO6]
    catch
        return x
    end
end

optics_table(tw) = hasproperty(tw, :table) ? tw.table : tw
parameter_gradient(x) = GTPSA.gradient(x, include_params=true)[7:end]

function sliced_drift(name, L, n)
    return [Drift(name=@sprintf("%s_%02d", name, i), L=L / n) for i in 1:n]
end

function sliced_quad(name, L, Kn1, n)
    return [Quadrupole(name=@sprintf("%s_%02d", name, i), L=L / n, Kn1=Kn1) for i in 1:n]
end

function build_tune_ring(k; knobs=nothing, n_quad=12, n_drift=24)
    strength(i) = isnothing(knobs) ? k[i] : k[i] + knobs[i]

    elements = vcat(
        sliced_quad("QF", 0.25, strength(1), n_quad),
        sliced_drift("D1", 1.00, n_drift),
        sliced_quad("QD", 0.25, strength(2), n_quad),
        sliced_drift("D2", 1.00, n_drift),
    )

    return Beamline(elements; species_ref=CH13_SPECIES, E_ref=CH13_E_REF)
end

function ring_tunes(k; knobs=nothing, descriptor=nothing, constants=true)
    ring = build_tune_ring(k; knobs=knobs)
    tw = isnothing(descriptor) ? twiss(ring) : twiss(ring, GTPSA_descriptor=descriptor)
    table = optics_table(tw)
    tunes = [table.phi_1[end], table.phi_2[end]]
    return constants ? tps_const.(tunes) : tunes
end

function tune_residual(k, target)
    return ring_tunes(k) - target
end

function tune_residual_with_knobs(k, target)
    return ring_tunes(k; knobs=DK_TUNE, descriptor=TUNE_DESCRIPTOR, constants=false) .- target
end

function tune_jacobian(k, target)
    residual = tune_residual_with_knobs(k, target)
    return vcat((parameter_gradient(r)' for r in residual)...)
end

function match_tunes(k0, target; max_iter=10, tolerance=1e-9)
    k = copy(k0)

    for iteration in 1:max_iter
        r = tune_residual(k, target)
        tunes = r + target
        @printf(
            "iteration %2d: Qx = %.9f, Qy = %.9f, residual norm = %.3e\n",
            iteration,
            tunes[1],
            tunes[2],
            norm(r),
        )

        norm(r) < tolerance && break

        J = tune_jacobian(k, target)
        step = -(J \ r)
        step_norm = norm(step)
        step_norm > 0.25 && (step .*= 0.25 / step_norm)
        k .+= step
    end

    return k
end

function plot_beta_beating(k_design, k_model)
    design_table = optics_table(twiss(build_tune_ring(k_design)))
    model_table = optics_table(twiss(build_tune_ring(k_model)))

    s = tps_const.(model_table.s)
    beta_x_design = tps_const.(design_table.beta_1)
    beta_y_design = tps_const.(design_table.beta_2)
    beta_x_model = tps_const.(model_table.beta_1)
    beta_y_model = tps_const.(model_table.beta_2)

    beating_x = 100 .* (beta_x_model .- beta_x_design) ./ beta_x_design
    beating_y = 100 .* (beta_y_model .- beta_y_design) ./ beta_y_design
    tune_difference = ring_tunes(k_model) - ring_tunes(k_design)

    fig = Figure(size=(800, 620))
    ax = Axis(
        fig[1, 1],
        xlabel="s [m]",
        ylabel="beta beating [%]",
        title="Exercise 3: beta(model) - beta(design)",
    )
    lines!(ax, s, beating_x; color=:royalblue3, linewidth=2, label="beta_1")
    lines!(ax, s, beating_y; color=:darkorange2, linewidth=2, label="beta_2")
    axislegend(ax; position=:rt)

    ax2 = Axis(
        fig[2, 1],
        ylabel="tune difference",
        title="Tune(model) - tune(design)",
        xticks=([1, 2], ["Qx", "Qy"]),
    )
    barplot!(ax2, [1, 2], tune_difference; color=[:royalblue3, :darkorange2])
    hlines!(ax2, [0.0]; color=:gray45, linestyle=:dash)

    return fig
end

K_DESIGN = [0.30, -0.30]
TARGET_TUNES = [0.08, 0.14]

println("Design tunes: ", ring_tunes(K_DESIGN))
K_MODEL = match_tunes(K_DESIGN, TARGET_TUNES)
println("\nMatched strengths: ", K_MODEL)
println("Matched tunes:     ", ring_tunes(K_MODEL))

exercise3_figure = plot_beta_beating(K_DESIGN, K_MODEL)
display(exercise3_figure)
