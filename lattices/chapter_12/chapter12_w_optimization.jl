using SciBmad
using LinearAlgebra
using Printf

# Compatibility helpers for the ESR lattice file generated from the Bmad model.
# Newer SciBmad versions use PhaseReference, while this generated file uses PhaseRef.
if !isdefined(Main, :PhaseRef)
    const PhaseRef = PhaseReference
end

if !isdefined(Main, :BeamlineChildRef)
    struct BeamlineChildRef
        beamline_index::Int
    end
end

# Older examples used findchildren(element, ring).  For this tutorial we only
# need the beamline index of a marker already present in ring.line.
findchildren(element, ring) = [
    BeamlineChildRef(i) for (i, candidate) in enumerate(ring.line) if candidate === element
]

const WORKSPACE_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const ESR_OPT_FILE = joinpath(WORKSPACE_ROOT, "esr-da-opt.jl")

# esr-da-opt.jl loads the full ESR ring, defines the CONTROLS object, applies
# chromaticity-compensating sextupole expressions, and installs trombone maps.
if !isdefined(Main, :ring_12_4)
    ring_12_4 = include(ESR_OPT_FILE)
end

const TPS_ZERO = [0, 0, 0, 0, 0, 0]
const TPS_PZ = [0, 0, 0, 0, 0, 1]

tps_constant(x) = x[TPS_ZERO]
tps_pz_coefficient(x) = x[TPS_PZ]

function w_ab(beta, alpha)
    beta0 = tps_constant(beta)
    alpha0 = tps_constant(alpha)
    dbeta_dpz = tps_pz_coefficient(beta)
    dalpha_dpz = tps_pz_coefficient(alpha)

    # These are the two components whose quadrature is the W function.
    A = dalpha_dpz - alpha0 / beta0 * dbeta_dpz
    B = dbeta_dpz / beta0

    return [A, B]
end

function marker_index(marker_name, ring)
    marker = getfield(Main, marker_name)
    return findchildren(marker, ring)[1].beamline_index
end

# The original Bmad 12.4 example targets IP6 and END_7.  The ESR Julia lattice
# does not contain a marker literally named END_7, so marc_end is used as the
# end-of-arc-7 marker between IP6 and IP8.
const W_TARGET_MARKERS = [:ip6, :marc_end]
const W_TARGET_INDICES = [marker_index(name, ring_12_4) for name in W_TARGET_MARKERS]
const CONTROL_NAMES = fieldnames(typeof(CONTROLS))

function zero_controls!()
    for name in CONTROL_NAMES
        setfield!(CONTROLS, name, 0.0)
    end
end

function apply_w_knobs!(x)
    zero_controls!()

    # Sextupole knobs in arcs 5 and 7.  These are the SciBmad/ESR analog of the
    # OF_5, OD_5, OF_7, and OD_7 variables in the Bmad tutorial.
    CONTROLS.dksf1_5 = x[1]
    CONTROLS.dksf2_5 = x[2]
    CONTROLS.dksd1_5 = x[3]
    CONTROLS.dksd2_5 = x[4]
    CONTROLS.dksf1_7 = x[5]
    CONTROLS.dksf2_7 = x[6]
    CONTROLS.dksd1_7 = x[7]
    CONTROLS.dksd2_7 = x[8]

    # Phase-trombone knobs around the 6 o'clock region.  The compensating IP4
    # trombone defined in esr-da-opt.jl keeps the total tune shift balanced.
    CONTROLS.dnux_mlrf_6 = x[9]
    CONTROLS.dnux_mlrr_6 = x[10]
    CONTROLS.dnuy_mlrf_6 = x[11]
    CONTROLS.dnuy_mlrr_6 = x[12]
end

function w_residual(x; ring=ring_12_4)
    apply_w_knobs!(x)
    tw = twiss(ring; GTPSA_descriptor=Descriptor(6, 2))
    table = tw.table
    residual = Float64[]

    for idx in W_TARGET_INDICES
        append!(residual, w_ab(table.beta_1[idx], table.alpha_1[idx]))
        append!(residual, w_ab(table.beta_2[idx], table.alpha_2[idx]))
    end

    return residual
end

function w_report(x; ring=ring_12_4)
    residual = w_residual(x; ring)
    @printf("total residual norm = %.6g\n", norm(residual))

    for k in eachindex(W_TARGET_MARKERS)
        offset = 4 * (k - 1) + 1
        wa = hypot(residual[offset], residual[offset + 1])
        wb = hypot(residual[offset + 2], residual[offset + 3])
        @printf("  %-8s  Wa = %10.6g   Wb = %10.6g\n", String(W_TARGET_MARKERS[k]), wa, wb)
    end

    return residual
end

function finite_difference_jacobian(x, residual0, steps)
    jacobian = zeros(length(residual0), length(x))

    for j in eachindex(x)
        x_step = copy(x)
        x_step[j] += steps[j]
        jacobian[:, j] = (w_residual(x_step) .- residual0) ./ steps[j]
    end

    return jacobian
end

function optimize_w_function(; n_iterations=4)
    x = zeros(12)

    # Sextupole knobs are in Kn2 units; trombone knobs are radians.
    finite_difference_steps = vcat(fill(0.005, 8), fill(0.0005, 4))

    for iter in 0:n_iterations
        residual = w_residual(x)
        @printf("iter %d: residual norm = %.6g\n", iter, norm(residual))
        @printf("  knobs = %s\n", repr(round.(x; digits=6)))

        iter == n_iterations && break

        jacobian = finite_difference_jacobian(x, residual, finite_difference_steps)
        gradient = jacobian' * residual

        # A conservative steepest-descent step is robust for this teaching
        # example. A production run would normally use a proper optimizer.
        direction = -gradient / norm(gradient)
        direction[1:8] .*= 0.12
        direction[9:12] .*= 0.02

        best_x = copy(x)
        best_norm = norm(residual)

        for scale in (1.0, 0.7, 0.5, 0.25, 0.1)
            trial_x = x .+ scale .* direction
            trial_norm = norm(w_residual(trial_x))
            @printf("  trial scale %.2f -> %.6g\n", scale, trial_norm)

            if trial_norm < best_norm
                best_x = trial_x
                break
            end
        end

        x = best_x
    end

    println("\nFinal W report:")
    w_report(x)

    return x
end

if abspath(PROGRAM_FILE) == @__FILE__
    optimize_w_function()
end
