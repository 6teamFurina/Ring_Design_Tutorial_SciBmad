# Chapter 9, Exercise 1: Set the synchrotron tune by varying RF voltage.
#
# The synchrotron tune Qz is the longitudinal analogue of the transverse
# betatron tunes Qx and Qy. It is obtained from the eigenvalues of the
# one-turn linear map in the longitudinal phase space (z, pz).
#
# This script:
#   1. builds the Chapter 5 ring with four FODORF cells at 10 o'clock;
#   2. calculates Qz from the GTPSA one-turn tracking Jacobian;
#   3. uses bisection to find the common RF0 voltage that gives Qz = 0.05.

using SciBmad
using GTPSA
using DifferentiationInterface
import DifferentiationInterface as DI
using LinearAlgebra
using Printf

const tutorial_root = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(tutorial_root, "lattices", "chapter_5", "chapter5_ring_definition.jl"))

const C5 = Chapter5Ring
const TARGET_QZ = 0.05

# Chapter 9 RF-cell dimensions. The 0.3 m spacing replaces the negative
# drift produced by the inconsistent expression in the original PDF.
const RF_L = 4.017
const RF_HARMON = 7560
const DRF_L = 0.3
const RF_VOLTAGE_START = 68.0 / 18.0 * 1e6

function make_FODORF_elements(voltage)
    rf() = RFCavity(
        name="RF0",
        L=RF_L,
        harmon=RF_HARMON,
        voltage=voltage,
    )
    drf() = Drift(name="DRF", L=DRF_L)

    return [
        C5.make_element(:QFSS),
        drf(), rf(), drf(), rf(), drf(),
        C5.make_element(:QDSS),
        drf(), rf(), drf(), rf(), drf(),
    ]
end

repeat_FODORF(n, voltage) =
    reduce(vcat, (make_FODORF_elements(voltage) for _ in 1:n))

function build_ring_with_rf(voltage)
    # The four RF cells straddle the SEXTANT9/SEXTANT11 boundary, matching
    # the placement specified in the original Chapter 9 example.
    sextant9_rf = vcat(
        C5.make_elements(C5.repeat_line(C5.FODOSSF, 4)),
        C5.make_elements(C5.SS_TO_ARCF),
        C5.make_elements(C5.repeat_line(C5.FODOAF, 20)),
        C5.make_elements(C5.ARC_TO_SSF),
        C5.make_elements(C5.repeat_line(C5.FODOSSF, 2)),
        repeat_FODORF(2, voltage),
    )

    sextant11_rf = vcat(
        repeat_FODORF(2, voltage),
        C5.make_elements(C5.repeat_line(C5.FODOSSR, 2)),
        C5.make_elements(C5.SS_TO_ARCR),
        C5.make_elements(C5.repeat_line(C5.FODOAR, 20)),
        C5.make_elements(C5.ARC_TO_SSR),
        C5.make_elements(C5.repeat_line(C5.FODOSSR, 4)),
    )

    elements = vcat(
        C5.make_elements(C5.SEXTANT1),
        C5.make_elements(C5.SEXTANT3),
        C5.make_elements(C5.SEXTANT5),
        C5.make_elements(C5.SEXTANT7),
        sextant9_rf,
        sextant11_rf,
    )
    return Beamline(elements; species_ref=C5.species_ref, E_ref=C5.E_ref)
end

function track_one_turn(v0, ring)
    v = similar(v0)
    v .= v0
    bunch = Bunch(v; species=ring.species_ref, p_over_q_ref=ring.p_over_q_ref)
    track!(bunch, ring)
    return copy(bunch.coords.v)
end

function transfer_matrix_gtpsa(ring; x0=zeros(6))
    prep = DI.prepare_jacobian(
        track_one_turn,
        AutoGTPSA(),
        x0,
        DI.Constant(ring),
    )
    return DI.jacobian(
        track_one_turn,
        prep,
        AutoGTPSA(),
        x0,
        DI.Constant(ring),
    )
end

function longitudinal_one_turn_matrix(ring)
    # The RF phase used here makes the phase-space origin the synchronous
    # closed orbit. GTPSA gives the local 2x2 longitudinal map directly:
    #
    #       [ z_final  ]       [ z_initial  ]
    #       [ pz_final ] = Mz * [ pz_initial ].
    M = transfer_matrix_gtpsa(ring)
    return M[[5, 6], [5, 6]]
end

function synchrotron_tune(voltage)
    ring = build_ring_with_rf(voltage)
    Mz = longitudinal_one_turn_matrix(ring)

    # For a stable oscillation, the eigenvalues are exp(+/- i*2*pi*Qz).
    # Using acos(trace(Mz)/2) returns the tune in the interval [0, 0.5].
    cos_phase = clamp(tr(Mz) / 2, -1.0, 1.0)
    Qz = acos(cos_phase) / (2pi)
    stable = abs(tr(Mz) / 2) < 1
    return Qz, Mz, stable
end

function bracket_target(target; voltage_start=RF_VOLTAGE_START, growth=2.0, max_steps=30)
    # Start from a small positive voltage and expand the upper bound until
    # the calculated tune is above the target.
    low = max(voltage_start / 100, 1.0)
    high = voltage_start
    q_low = first(synchrotron_tune(low))
    q_high = first(synchrotron_tune(high))

    for _ in 1:max_steps
        q_low <= target <= q_high && return low, high
        low, q_low = high, q_high
        high *= growth
        q_high = first(synchrotron_tune(high))
    end
    error("Could not bracket Qz = $target. Last result: Qz = $q_high at V = $high V.")
end

function optimize_rf_voltage(target=TARGET_QZ; tune_tol=1e-10, max_iter=80)
    low, high = bracket_target(target)

    for iteration in 1:max_iter
        voltage = (low + high) / 2
        qz, _, stable = synchrotron_tune(voltage)

        @printf("iter %2d: voltage = %.9e V, Qz = %.12f\n", iteration, voltage, qz)
        stable || error("Longitudinal map became unstable at V = $voltage V.")

        abs(qz - target) < tune_tol && return voltage
        qz < target ? (low = voltage) : (high = voltage)
    end
    return (low + high) / 2
end

println("Initial RF setting:")
qz_start, Mz_start, stable_start = synchrotron_tune(RF_VOLTAGE_START)
ring_start = build_ring_with_rf(RF_VOLTAGE_START)
num_rf_cavities = count(element -> element.name == "RF0", ring_start.line)
@printf("  RF0.voltage = %.9e V per cavity\n", RF_VOLTAGE_START)
@printf("  number of RF0 cavities = %d\n", num_rf_cavities)
@printf("  total installed RF voltage = %.9e V\n", num_rf_cavities * RF_VOLTAGE_START)
@printf("  Qz = %.12f, stable = %s\n", qz_start, stable_start)
println("  longitudinal one-turn matrix =")
display(Mz_start)

println("\nOptimizing common RF0 voltage for Qz = $(TARGET_QZ):")
RF_VOLTAGE_OPT = optimize_rf_voltage()

qz_final, Mz_final, stable_final = synchrotron_tune(RF_VOLTAGE_OPT)
println("\nFinal result:")
@printf("  RF0.voltage = %.12f V per cavity\n", RF_VOLTAGE_OPT)
@printf("  total installed RF voltage = %.12f V\n", num_rf_cavities * RF_VOLTAGE_OPT)
@printf("  Qz = %.12f, target = %.12f\n", qz_final, TARGET_QZ)
@printf("  stable = %s\n", stable_final)
println("  longitudinal one-turn matrix =")
display(Mz_final)
