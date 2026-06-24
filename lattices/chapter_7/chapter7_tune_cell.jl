# Chapter 7: Tune Cell
#
# Insert a special tune cell into the 2 o'clock straight section and vary six
# quadrupole-family strengths to reach Qx = 54.08 and Qy = 54.14. Four families
# match the modified straight section to the neighboring arc, while two
# families set the optics and phase advance of all FODO cells in the straight.

using SciBmad
using GTPSA
using LinearAlgebra
using Printf

const tutorial_root = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(tutorial_root, "lattices", "chapter_5", "chapter5_ring_definition.jl"))
include(joinpath(tutorial_root, "lattices", "chapter_6", "chapter6_IR_solution.jl"))

const C5 = Chapter5Ring
const TARGET_TUNES = [54.08, 54.14]
const TUNE_DESCRIPTOR = Descriptor(6, 2, 6, 1)
const DK_TUNE = params(TUNE_DESCRIPTOR)

# The six independent knobs are ordered as:
# QFF2_2, QDF2_2, QFF3_2, QDF3_2, QFSS_2, QDSS_2.
const K_START = [
    C5.quad_strength[:QFF2],
    C5.quad_strength[:QDF2],
    C5.quad_strength[:QFF3],
    C5.quad_strength[:QDF3],
    C5.quad_strength[:QFSS],
    C5.quad_strength[:QDSS],
]

function special_quad_strength(k; knobs=nothing)
    strength(i) = isnothing(knobs) ? k[i] : k[i] + knobs[i]
    return Dict(
        :QFF2_2 => strength(1),
        :QDF2_2 => strength(2),
        :QFF3_2 => strength(3),
        :QDF3_2 => strength(4),
        :QFSS_2 => strength(5),
        :QDSS_2 => strength(6),
    )
end

function make_tune_element(name::Symbol, k; knobs=nothing)
    strengths = special_quad_strength(k; knobs=knobs)
    haskey(strengths, name) &&
        return Quadrupole(name=String(name), L=C5.L_quad, Kn1=strengths[name])
    return C5.make_element(name)
end

make_tune_elements(line, k; knobs=nothing) = [make_tune_element(name, k; knobs=knobs) for name in line]

# All physical QFSS_2 magnets share k[5], and all QDSS_2 magnets share k[6].
const FODOSSF_2 = [:QFSS_2, :D1, :DB, :D2, :QDSS_2, :D1, :DB, :D2]
const FODOSSR_2 = [:QFSS_2, :D2, :DB, :D1, :QDSS_2, :D2, :DB, :D1]

# The left match contains four independent matching families. The right match
# reuses the same strengths in mirror-reversed order, so it adds no new knobs.
const ARC_TO_SSF_2 = [
    :QF, :D1, :BH, :D2, :QD, :D1, :BH, :D2,
    :QFF1, :D1, :BH, :D2, :QDF1, :D1, :BH, :D2,
    :QFF2_2, :D1, :DB, :D2, :QDF2_2, :D1, :DB, :D2,
    :QFF3_2, :D1, :DB, :D2, :QDF3_2, :D1, :DB, :D2,
]

const SS_TO_ARCR_2 = [
    :QFSS_2, :D2, :DB, :D1, :QDF3_2, :D2, :DB, :D1,
    :QFF3_2, :D2, :DB, :D1, :QDF2_2, :D2, :DB, :D1,
    :QFF2_2, :D2, :BH, :D1, :QDF1, :D2, :BH, :D1,
    :QFF1, :D2, :BH, :D1, :QD, :D2, :BH, :D1,
]

function build_IPF_elements()
    return [
        Quadrupole(name="QEF1", L=0.5, Kn1=K_QEF1),
        C5.make_element(:D1), C5.make_element(:DB), C5.make_element(:D2),
        Quadrupole(name="QEF2", L=0.5, Kn1=K_QEF2),
        C5.make_element(:D1), C5.make_element(:DB), C5.make_element(:D2),
        Drift(name="DEF1", L=20.46),
        Quadrupole(name="QEF3", L=1.6, Kn1=K_QEF3),
        Drift(name="DEF2", L=3.76),
        Quadrupole(name="QEF4", L=1.2, Kn1=K_QEF4),
        Drift(name="DEF3", L=5.8),
        Marker(name="IP6"),
    ]
end

function build_IPR_elements()
    return [
        Drift(name="DER3", L=5.3),
        Quadrupole(name="QER4", L=1.8, Kn1=K_QER4),
        Drift(name="DER2", L=0.5),
        Quadrupole(name="QER3", L=1.4, Kn1=K_QER3),
        Drift(name="DER1", L=23.82),
        Quadrupole(name="QER2", L=0.5, Kn1=K_QER2),
        C5.make_element(:D2), C5.make_element(:DB), C5.make_element(:D1),
        Quadrupole(name="QER1", L=0.5, Kn1=K_QER1),
        C5.make_element(:D2), C5.make_element(:DB), C5.make_element(:D1),
    ]
end

function build_ring_with_tune_cell(k; knobs=nothing)
    # The tune cell occupies the straight between the end of sextant 1 and the
    # beginning of sextant 3. The low-beta IR from Chapter 6 remains at IP6.
    sextant1_tune = vcat(
        C5.make_elements(C5.repeat_line(C5.FODOSSF, 4)),
        C5.make_elements(C5.SS_TO_ARCF),
        C5.make_elements(C5.repeat_line(C5.FODOAF, 20)),
        make_tune_elements(ARC_TO_SSF_2, k; knobs=knobs),
        make_tune_elements(C5.repeat_line(FODOSSF_2, 4), k; knobs=knobs),
    )

    sextant3_tune = vcat(
        make_tune_elements(C5.repeat_line(FODOSSR_2, 4), k; knobs=knobs),
        make_tune_elements(SS_TO_ARCR_2, k; knobs=knobs),
        C5.make_elements(C5.repeat_line(C5.FODOAR, 20)),
        C5.make_elements(C5.ARC_TO_SSR),
        C5.make_elements(C5.repeat_line(C5.FODOSSR, 4)),
    )

    sextant5_ir = vcat(
        C5.make_elements(C5.repeat_line(C5.FODOSSF, 4)),
        C5.make_elements(C5.SS_TO_ARCF),
        C5.make_elements(C5.repeat_line(C5.FODOAF, 20)),
        C5.make_elements(C5.ARC_TO_SSF),
        C5.make_elements(C5.FODOSSF),
        build_IPF_elements(),
    )

    sextant7_ir = vcat(
        build_IPR_elements(),
        C5.make_elements(C5.FODOSSR),
        C5.make_elements(C5.SS_TO_ARCR),
        C5.make_elements(C5.repeat_line(C5.FODOAR, 20)),
        C5.make_elements(C5.ARC_TO_SSR),
        C5.make_elements(C5.repeat_line(C5.FODOSSR, 4)),
    )

    elements = vcat(
        sextant1_tune,
        sextant3_tune,
        sextant5_ir,
        sextant7_ir,
        C5.make_elements(C5.SEXTANT9),
        C5.make_elements(C5.SEXTANT11),
    )

    ring = Beamline(elements; species_ref=C5.species_ref, E_ref=C5.E_ref)
    return ring, elements
end

optics_table(tw) = hasproperty(tw, :table) ? tw.table : tw

const zero6 = [0, 0, 0, 0, 0, 0]

function tps_const(x)
    try
        return x[zero6]
    catch
        return x
    end
end

parameter_gradient(x) = GTPSA.gradient(x, include_params=true)[7:end]

function straight_cell_phase_advances(k)
    # Read the phase advance across one center FODO period from the complete
    # stable ring. Calling twiss on an isolated no-bend cell can leave the
    # longitudinal normal form degenerate even when its transverse motion is
    # stable.
    ring, elements = build_ring_with_tune_cell(k)
    table = optics_table(twiss(ring))
    qf_indices = findall(ele -> ele.name == "QFSS_2", elements)
    row4 = qf_indices[4] + 1
    row5 = qf_indices[5] + 1
    return 360 .* tps_const.([
        table.phi_1[row5] - table.phi_1[row4],
        table.phi_2[row5] - table.phi_2[row4],
    ])
end

function tune_cell_metrics(k; knobs=nothing, descriptor=nothing, constants=true)
    ring, elements = build_ring_with_tune_cell(k; knobs=knobs)
    tw = isnothing(descriptor) ? twiss(ring) : twiss(ring, GTPSA_descriptor=descriptor)
    table = optics_table(tw)

    # The fourth and fifth QFSS_2 magnets are one FODO period apart at the
    # center. Equal Twiss values at these locations enforce local periodicity.
    qf_indices = findall(ele -> ele.name == "QFSS_2", elements)
    length(qf_indices) == 9 || error("Expected nine QFSS_2 elements.")
    row4 = qf_indices[4] + 1
    row5 = qf_indices[5] + 1

    periodicity = [
        (table.beta_1[row4] - table.beta_1[row5]) / table.beta_1[row5],
        table.alpha_1[row4] - table.alpha_1[row5],
        (table.beta_2[row4] - table.beta_2[row5]) / table.beta_2[row5],
        table.alpha_2[row4] - table.alpha_2[row5],
    ]
    tunes = [table.phi_1[end], table.phi_2[end]]

    if constants
        periodicity = tps_const.(periodicity)
        tunes = tps_const.(tunes)
    end

    return (periodicity=periodicity, tunes=tunes, table=table, elements=elements)
end

function tune_cell_residual(k)
    metrics = tune_cell_metrics(k)
    return vcat(metrics.periodicity, metrics.tunes - TARGET_TUNES)
end

function tune_cell_residual_with_knobs(k)
    metrics = tune_cell_metrics(
        k;
        knobs=DK_TUNE,
        descriptor=TUNE_DESCRIPTOR,
        constants=false,
    )
    return vcat(metrics.periodicity, metrics.tunes .- TARGET_TUNES)
end

function tune_cell_residual_jacobian(k)
    residual = tune_cell_residual_with_knobs(k)
    return vcat((parameter_gradient(r)' for r in residual)...)
end

function damped_least_squares(
    f,
    x0;
    jacobian=tune_cell_residual_jacobian,
    maxiter=60,
    tol=1e-11,
    lambda0=1e-3,
)
    x = copy(x0)
    lambda = lambda0

    for iter in 1:maxiter
        r = f(x)
        merit_now = sum(abs2, r)
        J = jacobian(x)
        step = -(J' * J + lambda * I) \ (J' * r)
        trial = x + step
        merit_trial = sum(abs2, f(trial))

        @printf(
            "iter %2d  merit = %.6e  step = %.3e  lambda = %.3e\n",
            iter, merit_now, norm(step), lambda,
        )

        if merit_trial < merit_now
            x = trial
            lambda /= 10
        else
            lambda *= 10
        end

        norm(r) < tol && break
        norm(step) < tol && break
    end
    return x
end

initial = tune_cell_metrics(K_START)
println("Initial tune-cell optics:")
@printf("  Qx = %.12f, Qy = %.12f\n", initial.tunes...)
@printf("  straight-cell phase advances = %.9f deg, %.9f deg\n", straight_cell_phase_advances(K_START)...)
println("  periodicity residual = ", initial.periodicity)

K_TUNE = damped_least_squares(tune_cell_residual, K_START)
final = tune_cell_metrics(K_TUNE)

println("\nOptimized tune-cell strengths:")
for (name, value) in zip(
    ["K_QFF2_2", "K_QDF2_2", "K_QFF3_2", "K_QDF3_2", "K_QFSS_2", "K_QDSS_2"],
    K_TUNE,
)
    @printf("  %-10s = %+.15f\n", name, value)
end

println("\nFinal checks:")
@printf("  Qx = %.12f  target = %.12f\n", final.tunes[1], TARGET_TUNES[1])
@printf("  Qy = %.12f  target = %.12f\n", final.tunes[2], TARGET_TUNES[2])
@printf("  straight-cell phase advances = %.9f deg, %.9f deg\n", straight_cell_phase_advances(K_TUNE)...)
println("  periodicity residual = ", final.periodicity)
println("  complete residual    = ", tune_cell_residual(K_TUNE))
