# Chapter 6, Exercise 5:
# Replace the 6 o'clock straight section of the Chapter 5 ring with the optimized interaction region.

using SciBmad
using Printf

const tutorial_root = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(tutorial_root, "lattices", "chapter_5", "chapter5_ring_definition.jl"))
include(joinpath(tutorial_root, "lattices", "chapter_6", "chapter6_IR_solution.jl"))

const C5 = Chapter5Ring

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
    # IP6 is already the final element of IPF, so IPR begins with DER3 here.
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

# Sextants 1, 3, 9, and 11 remain unchanged. Sextant 5 ends at IP6 through IPF, and Sextant 7 starts immediately downstream of IP6 through IPR.
SEXTANT5_IR = vcat(
    C5.make_elements(C5.repeat_line(C5.FODOSSF, 4)),
    C5.make_elements(C5.SS_TO_ARCF),
    C5.make_elements(C5.repeat_line(C5.FODOAF, 20)),
    C5.make_elements(C5.ARC_TO_SSF),
    C5.make_elements(C5.FODOSSF),
    build_IPF_elements(),
)

SEXTANT7_IR = vcat(
    build_IPR_elements(),
    C5.make_elements(C5.FODOSSR),
    C5.make_elements(C5.SS_TO_ARCR),
    C5.make_elements(C5.repeat_line(C5.FODOAR, 20)),
    C5.make_elements(C5.ARC_TO_SSR),
    C5.make_elements(C5.repeat_line(C5.FODOSSR, 4)),
)

ring_with_ir_elements = vcat(
    C5.make_elements(C5.SEXTANT1),
    C5.make_elements(C5.SEXTANT3),
    SEXTANT5_IR,
    SEXTANT7_IR,
    C5.make_elements(C5.SEXTANT9),
    C5.make_elements(C5.SEXTANT11),
)

ring_with_ir = Beamline(
    ring_with_ir_elements;
    species_ref=C5.species_ref,
    E_ref=C5.E_ref,
)

# The IR replaces six ordinary straight FODO cells with equal total length, so the ring circumference and total bend remain unchanged.
standard_fodo_length = 2C5.L_quad + 2(C5.D1_len + C5.DB_len + C5.D2_len)
ipf_length = 49.23
ipr_length = 49.23

function chapter5_element_length(name::Symbol)
    haskey(C5.quad_strength, name) && return C5.L_quad
    name == :D1 && return C5.D1_len
    name == :D2 && return C5.D2_len
    name == :DB && return C5.DB_len
    return C5.B_len
end

original_circumference = sum(chapter5_element_length(name) for name in C5.RING_NAMES)
new_circumference = original_circumference - 6standard_fodo_length + ipf_length + ipr_length

println("Constructed the Chapter 6 ring with the low-beta IR.")
@printf("  Original ring elements: %d\n", length(C5.RING_NAMES))
@printf("  Ring-with-IR elements:  %d\n", length(ring_with_ir_elements))
@printf("  Original circumference: %.9f m\n", original_circumference)
@printf("  New circumference:      %.9f m\n", new_circumference)
@printf("  Circumference change:   %+.3e m\n", new_circumference - original_circumference)
println("  IP6 marker count:       ", count(ele -> ele.name == "IP6", ring_with_ir_elements))

println("\nCalculating the periodic optics of the ring with the IR...")
tw_ring_with_ir = twiss(ring_with_ir)
println("Periodic full-ring Twiss calculation completed.")

# The Twiss table contains the entrance point followed by one row after each element. Therefore, the row immediately after the IP6 marker is index + 1.
ip_element_index = findfirst(ele -> ele.name == "IP6", ring_with_ir_elements)
ip_table_index = ip_element_index + 1

println("\nPeriodic full-ring Twiss at IP6:")
@printf(
    "  beta_1 = %.12f, alpha_1 = %.3e\n",
    tw_ring_with_ir.beta_1[ip_table_index],
    tw_ring_with_ir.alpha_1[ip_table_index],
)
@printf(
    "  beta_2 = %.12f, alpha_2 = %.3e\n",
    tw_ring_with_ir.beta_2[ip_table_index],
    tw_ring_with_ir.alpha_2[ip_table_index],
)
