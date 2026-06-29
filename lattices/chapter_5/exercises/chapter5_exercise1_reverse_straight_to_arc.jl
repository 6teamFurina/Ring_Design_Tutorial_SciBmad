# Chapter 5, Exercise 1:
# Construct the reverse straight-section-to-arc connection SS_TO_ARCR.
#
# SS_TO_ARCR is the mirror image of the previously optimized forward
# arc-to-straight connection ARC_TO_SSF. Therefore, no new optimization is
# needed: copy the forward Chapter 2 and Chapter 3 quadrupole strengths and
# reverse the drift ordering from D1-...-D2 to D2-...-D1.

using SciBmad
using Printf

const tutorial_root = normpath(joinpath(@__DIR__, "..", "..", ".."))
const chapter1_solution = joinpath(tutorial_root, "lattices", "chapter_1", "chapter1_fodoF_solution.jl")
const chapter2_solution = joinpath(tutorial_root, "lattices", "chapter_2", "chapter2_dispsupF_solution.jl")
const chapter3_solution = joinpath(tutorial_root, "lattices", "chapter_3", "chapter3_mSSF_solution.jl")

isfile(chapter1_solution) ||
    error("Cannot find chapter1_fodoF_solution.jl. Run Chapter 1 first.")
isfile(chapter2_solution) ||
    error("Cannot find chapter2_dispsupF_solution.jl. Run Chapter 2 first.")
isfile(chapter3_solution) ||
    error("Cannot find chapter3_mSSF_solution.jl. Run Chapter 3 first.")

include(chapter1_solution)
include(chapter2_solution)
include(chapter3_solution)

const L_quad = 0.5
const D1_len = 0.609
const D2_len = 1.241
const DB_len = 5.855
const B_len = 6.86
const BH_angle = (pi / 132) / 2
const K_ss = 0.351957452649287

const species_ref = Species("electron")
const E_ref = 18e9

D1() = Drift(name="D1", L=D1_len)
D2() = Drift(name="D2", L=D2_len)
DB() = Drift(name="DB", L=DB_len)
BH() = SBend(name="BH", L=B_len, angle=BH_angle)

QFSS() = Quadrupole(name="QFSS", L=L_quad, Kn1=K_ss)
QD() = Quadrupole(name="QD", L=L_quad, Kn1=kQD_arc)
QFF1() = Quadrupole(name="QFF1", L=L_quad, Kn1=kQFF1)
QDF1() = Quadrupole(name="QDF1", L=L_quad, Kn1=kQDF1)
QFF2() = Quadrupole(name="QFF2", L=L_quad, Kn1=K_QFF2)
QDF2() = Quadrupole(name="QDF2", L=L_quad, Kn1=K_QDF2)
QFF3() = Quadrupole(name="QFF3", L=L_quad, Kn1=K_QFF3)
QDF3() = Quadrupole(name="QDF3", L=L_quad, Kn1=K_QDF3)

# Match the reverse straight section to the dispersion creator.
#
# This is the mirror of the forward matching section MSSF. Starting at the
# straight-section QF, the matching quadrupoles appear in reverse order.
MDCR = [
    QFSS(), D2(), DB(), D1(),
    QDF3(), D2(), DB(), D1(),
    QFF3(), D2(), DB(), D1(),
    QDF2(), D2(), DB(), D1(),
]

# Reverse dispersion creator.
#
# This is the mirror of the forward dispersion suppressor DISPSUPF. It starts
# with QFF2, which is shared with the end of MDCR, and gradually creates the
# periodic reverse-arc dispersion using four half-angle bends.
DISPCRER = [
    QFF2(), D2(), BH(), D1(),
    QDF1(), D2(), BH(), D1(),
    QFF1(), D2(), BH(), D1(),
    QD(), D2(), BH(), D1(),
]

# Reverse straight section to reverse arc.
SS_TO_ARCR_elements = vcat(MDCR, DISPCRER)
SS_TO_ARCR = Beamline(
    SS_TO_ARCR_elements;
    species_ref=species_ref,
    E_ref=E_ref,
)

println("Constructed SS_TO_ARCR:")
@printf("  MDCR elements:      %d\n", length(MDCR))
@printf("  DISPCRER elements:  %d\n", length(DISPCRER))
@printf("  Total elements:     %d\n", length(SS_TO_ARCR_elements))
@printf("  Half bends:         %d\n", 4)
println()
println("SS_TO_ARCR is ready to be placed before a reverse arc FODO line.")
