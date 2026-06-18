include(joinpath(@__DIR__, "chapter14_common.jl"))

using LinearAlgebra

function tps_const(x)
    try
        return x[[0, 0, 0, 0, 0, 0]]
    catch
        return x
    end
end

optics_table(tw) = hasproperty(tw, :table) ? tw.table : tw

function sliced_drift(name, L, n)
    return [Drift(name=@sprintf("%s_%02d", name, i), L=L / n) for i in 1:n]
end

function sliced_quad(name, L, Kn1, n)
    return [Quadrupole(name=@sprintf("%s_%02d", name, i), L=L / n, Kn1=Kn1) for i in 1:n]
end

function build_tune_ring(kqf, kqd; n_quad=12, n_drift=24)
    elements = vcat(
        sliced_quad("QF", 0.25, kqf, n_quad),
        sliced_drift("D1", 1.00, n_drift),
        sliced_quad("QD", 0.25, kqd, n_quad),
        sliced_drift("D2", 1.00, n_drift),
    )

    return Beamline(elements; species_ref=CH13_SPECIES, E_ref=CH13_E_REF)
end

function ring_tunes(k)
    table = optics_table(twiss(build_tune_ring(k[1], k[2])))
    return tps_const.([table.phi_1[end], table.phi_2[end]])
end

function finite_difference_jacobian(f, x; step=1e-5)
    r0 = f(x)
    J = zeros(length(r0), length(x))

    for j in eachindex(x)
        xp = copy(x)
        xm = copy(x)
        xp[j] += step
        xm[j] -= step
        J[:, j] = (f(xp) - f(xm)) / (2step)
    end

    return J
end

function match_tunes(k0, target; max_iter=10, tolerance=1e-9)
    k = copy(k0)
    residual(k) = ring_tunes(k) - target

    for iteration in 1:max_iter
        r = residual(k)
        tunes = r + target
        @printf(
            "iteration %2d: Qx = %.9f, Qy = %.9f, residual norm = %.3e\n",
            iteration,
            tunes[1],
            tunes[2],
            norm(r),
        )

        norm(r) < tolerance && break

        J = finite_difference_jacobian(residual, k)
        step = -(J \ r)
        step_norm = norm(step)
        step_norm > 0.25 && (step .*= 0.25 / step_norm)
        k .+= step
    end

    return k
end

function plot_beta_beating(k_design, k_model)
    design_table = optics_table(twiss(build_tune_ring(k_design[1], k_design[2])))
    model_table = optics_table(twiss(build_tune_ring(k_model[1], k_model[2])))

    s = tps_const.(model_table.s)
    beta_x_design = tps_const.(design_table.beta_1)
    beta_y_design = tps_const.(design_table.beta_2)
    beta_x_model = tps_const.(model_table.beta_1)
    beta_y_model = tps_const.(model_table.beta_2)

    beating_x = 100 .* (beta_x_model .- beta_x_design) ./ beta_x_design
    beating_y = 100 .* (beta_y_model .- beta_y_design) ./ beta_y_design

    fig = Figure(size=(780, 420))
    ax = Axis(
        fig[1, 1],
        xlabel="s [m]",
        ylabel="beta beating [%]",
        title="Exercise 3: beta(model) - beta(design)",
    )
    lines!(ax, s, beating_x; color=:royalblue3, linewidth=2, label="beta_1")
    lines!(ax, s, beating_y; color=:darkorange2, linewidth=2, label="beta_2")
    axislegend(ax; position=:rt)
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
