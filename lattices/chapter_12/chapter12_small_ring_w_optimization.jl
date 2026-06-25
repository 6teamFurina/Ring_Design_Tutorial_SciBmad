using SciBmad
using LinearAlgebra
using Printf

include(joinpath(@__DIR__, "trombone_utils.jl"))

if !isdefined(Main, :SmallRingControls)
    @kwdef mutable struct SmallRingControls
        sf5::Float64 = 0.0
        sd5::Float64 = 0.0
        sf7::Float64 = 0.0
        sd7::Float64 = 0.0

        dnux_mlrf_6::Float64 = 0.0
        dnux_mlrr_6::Float64 = 0.0
        dnuy_mlrf_6::Float64 = 0.0
        dnuy_mlrr_6::Float64 = 0.0
    end
end

if !isdefined(Main, :SMALL_RING_CONTROLS)
    const SMALL_RING_CONTROLS = SmallRingControls()
end

small_dnu_ip4_x() = -(SMALL_RING_CONTROLS.dnux_mlrf_6 + SMALL_RING_CONTROLS.dnux_mlrr_6)
small_dnu_ip4_y() = -(SMALL_RING_CONTROLS.dnuy_mlrf_6 + SMALL_RING_CONTROLS.dnuy_mlrr_6)

function small_ring_cell(label, bend_angle; sf=0.0, sd=0.0)
    return SciBmad.LineElement[
        Quadrupole(name="QF_$label", L=0.25, Kn1=0.8),
        Drift(name="D1_$label", L=0.20),
        Sextupole(name="SF_$label", L=0.10, Kn2=sf),
        Drift(name="D2_$label", L=0.20),
        SBend(name="B_$label", L=1.0, angle=bend_angle),
        Drift(name="D3_$label", L=0.20),
        Sextupole(name="SD_$label", L=0.10, Kn2=sd),
        Drift(name="D4_$label", L=0.20),
        Quadrupole(name="QD_$label", L=0.25, Kn1=-0.8),
        Drift(name="D5_$label", L=0.40),
    ]
end

function build_small_ring_for_w()
    ip6 = Marker(name="ip6")
    mlrr_6 = Marker(name="mlrr_6")
    marc_end = Marker(name="marc_end")
    ip8 = Marker(name="ip8")
    ip10 = Marker(name="ip10")
    ip12 = Marker(name="ip12")
    ip2 = Marker(name="ip2")
    ip4 = Marker(name="ip4")
    mlrf_6 = Marker(name="mlrf_6")

    bend_angle = 2pi / 12
    elements = SciBmad.LineElement[
        ip6,
        Drift(name="IR6_R", L=0.35),
        mlrr_6,
    ]

    # Arc 7 is immediately after IP6 in the ESR ring and contains one optimized
    # sextupole pair in this compact teaching model.
    append!(elements, small_ring_cell("7A", bend_angle; sf=DefExpr(() -> 60.0 + SMALL_RING_CONTROLS.sf7), sd=DefExpr(() -> -60.0 + SMALL_RING_CONTROLS.sd7)))
    append!(elements, small_ring_cell("7B", bend_angle; sf=DefExpr(() -> 60.0 - SMALL_RING_CONTROLS.sf7), sd=DefExpr(() -> -60.0 - SMALL_RING_CONTROLS.sd7)))
    push!(elements, marc_end)

    # The remaining arcs keep the ring closed and provide a similar clock-order
    # layout to the full ESR model.
    push!(elements, Drift(name="IR8_L", L=0.25), ip8)
    append!(elements, small_ring_cell("9A", bend_angle; sf=50.0, sd=-50.0))
    append!(elements, small_ring_cell("9B", bend_angle; sf=50.0, sd=-50.0))
    push!(elements, ip10)
    append!(elements, small_ring_cell("11A", bend_angle; sf=50.0, sd=-50.0))
    append!(elements, small_ring_cell("11B", bend_angle; sf=50.0, sd=-50.0))
    push!(elements, ip12)
    append!(elements, small_ring_cell("1A", bend_angle; sf=50.0, sd=-50.0))
    append!(elements, small_ring_cell("1B", bend_angle; sf=50.0, sd=-50.0))
    push!(elements, ip2)
    append!(elements, small_ring_cell("3A", bend_angle; sf=50.0, sd=-50.0))
    append!(elements, small_ring_cell("3B", bend_angle; sf=50.0, sd=-50.0))
    push!(elements, ip4)

    # Arc 5 is immediately before IP6 in the ESR ring.  The two cells use
    # opposite changes, mimicking the paired overlay style of the Bmad example.
    append!(elements, small_ring_cell("5A", bend_angle; sf=DefExpr(() -> 60.0 + SMALL_RING_CONTROLS.sf5), sd=DefExpr(() -> -60.0 + SMALL_RING_CONTROLS.sd5)))
    append!(elements, small_ring_cell("5B", bend_angle; sf=DefExpr(() -> 60.0 - SMALL_RING_CONTROLS.sf5), sd=DefExpr(() -> -60.0 - SMALL_RING_CONTROLS.sd5)))
    push!(elements, mlrf_6, Drift(name="IR6_L", L=0.35), ip6)

    ring = Beamline(elements; species_ref=Species("electron"), E_ref=3e9)
    tw = twiss(ring)

    attach_trombone!(
        mlrr_6,
        ring,
        tw.table,
        DefExpr(() -> SMALL_RING_CONTROLS.dnux_mlrr_6),
        DefExpr(() -> SMALL_RING_CONTROLS.dnuy_mlrr_6),
    )
    attach_trombone!(
        mlrf_6,
        ring,
        tw.table,
        DefExpr(() -> SMALL_RING_CONTROLS.dnux_mlrf_6),
        DefExpr(() -> SMALL_RING_CONTROLS.dnuy_mlrf_6),
    )
    attach_trombone!(
        ip4,
        ring,
        tw.table,
        DefExpr(small_dnu_ip4_x),
        DefExpr(small_dnu_ip4_y),
    )

    return ring, Dict(
        :ip6 => ip6,
        :mlrr_6 => mlrr_6,
        :marc_end => marc_end,
        :ip8 => ip8,
        :ip4 => ip4,
        :mlrf_6 => mlrf_6,
    )
end

if !isdefined(Main, :small_ring_12_4)
    small_ring_12_4, small_ring_markers = build_small_ring_for_w()
end

const SMALL_TPS_ZERO = [0, 0, 0, 0, 0, 0]
const SMALL_TPS_PZ = [0, 0, 0, 0, 0, 1]

small_tps_constant(x) = x[SMALL_TPS_ZERO]
small_tps_pz_coefficient(x) = x[SMALL_TPS_PZ]

function small_w_ab(beta, alpha)
    beta0 = small_tps_constant(beta)
    alpha0 = small_tps_constant(alpha)
    dbeta_dpz = small_tps_pz_coefficient(beta)
    dalpha_dpz = small_tps_pz_coefficient(alpha)

    A = dalpha_dpz - alpha0 / beta0 * dbeta_dpz
    B = dbeta_dpz / beta0
    return [A, B]
end

function small_marker_index(marker, ring)
    for (idx, candidate) in enumerate(ring.line)
        candidate === marker && return idx
    end
    error("Marker is not in this ring.")
end

const SMALL_W_TARGET_MARKERS = [:ip6, :marc_end]
const SMALL_W_TARGET_INDICES = [
    small_marker_index(small_ring_markers[name], small_ring_12_4)
    for name in SMALL_W_TARGET_MARKERS
]

function apply_small_w_knobs!(x)
    SMALL_RING_CONTROLS.sf5 = x[1]
    SMALL_RING_CONTROLS.sd5 = x[2]
    SMALL_RING_CONTROLS.sf7 = x[3]
    SMALL_RING_CONTROLS.sd7 = x[4]

    SMALL_RING_CONTROLS.dnux_mlrf_6 = x[5]
    SMALL_RING_CONTROLS.dnux_mlrr_6 = x[6]
    SMALL_RING_CONTROLS.dnuy_mlrf_6 = x[7]
    SMALL_RING_CONTROLS.dnuy_mlrr_6 = x[8]
end

function small_w_residual(x; ring=small_ring_12_4)
    apply_small_w_knobs!(x)
    tw = twiss(ring; GTPSA_descriptor=Descriptor(6, 2))
    table = tw.table
    residual = Float64[]

    for idx in SMALL_W_TARGET_INDICES
        append!(residual, small_w_ab(table.beta_1[idx], table.alpha_1[idx]))
        append!(residual, small_w_ab(table.beta_2[idx], table.alpha_2[idx]))
    end

    return residual
end

function small_w_report(x)
    residual = small_w_residual(x)
    @printf("total residual norm = %.6g\n", norm(residual))

    for k in eachindex(SMALL_W_TARGET_MARKERS)
        offset = 4 * (k - 1) + 1
        wa = hypot(residual[offset], residual[offset + 1])
        wb = hypot(residual[offset + 2], residual[offset + 3])
        @printf("  %-8s  Wa = %10.6g   Wb = %10.6g\n", String(SMALL_W_TARGET_MARKERS[k]), wa, wb)
    end

    return residual
end

function small_finite_difference_jacobian(x, residual0, steps)
    jacobian = zeros(length(residual0), length(x))

    for j in eachindex(x)
        x_step = copy(x)
        x_step[j] += steps[j]
        jacobian[:, j] = (small_w_residual(x_step) .- residual0) ./ steps[j]
    end

    return jacobian
end

function optimize_small_w_function(; n_iterations=6)
    x = zeros(8)
    finite_difference_steps = [0.005, 0.005, 0.005, 0.005, 0.0005, 0.0005, 0.0005, 0.0005]

    for iter in 0:n_iterations
        residual = small_w_residual(x)
        @printf("iter %d: residual norm = %.6g\n", iter, norm(residual))
        @printf("  knobs = %s\n", repr(round.(x; digits=6)))

        iter == n_iterations && break

        jacobian = small_finite_difference_jacobian(x, residual, finite_difference_steps)
        gradient = jacobian' * residual
        direction = -gradient / norm(gradient)

        # Sextupole strengths are much larger-scale variables than trombone
        # phase shifts, so use different teaching step sizes. This also keeps
        # the arc-5/arc-7 sextupoles visibly involved in the optimization.
        direction[1:4] .*= 100.0
        direction[5:8] .*= 0.02

        best_x = copy(x)
        best_norm = norm(residual)

        for scale in (1.0, 0.7, 0.5, 0.25, 0.1)
            trial_x = x .+ scale .* direction
            trial_norm = norm(small_w_residual(trial_x))
            @printf("  trial scale %.2f -> %.6g\n", scale, trial_norm)

            if trial_norm < best_norm
                best_x = trial_x
                break
            end
        end

        x = best_x
    end

    println("\nFinal small-ring W report:")
    small_w_report(x)

    return x
end

if abspath(PROGRAM_FILE) == @__FILE__
    @printf("small ring elements: %d\n", length(small_ring_12_4.line))
    println("Initial small-ring W report:")
    small_w_report(zeros(8))
    optimize_small_w_function()
end
