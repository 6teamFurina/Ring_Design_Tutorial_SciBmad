# Chapter 5, Exercise 2:
# Gather the optimized sections and construct the complete six-sextant ring.

using SciBmad
using Printf

const tutorial_root = normpath(joinpath(@__DIR__, "..", "..", ".."))

# Gather the optimized forward and reverse quadrupole strengths.
solution_files = [
    joinpath("lattices", "chapter_1", "chapter1_fodoF_solution.jl"),
    joinpath("lattices", "chapter_1", "chapter1_fodoR_solution.jl"),
    joinpath("lattices", "chapter_2", "chapter2_dispsupF_solution.jl"),
    joinpath("lattices", "chapter_2", "chapter2_dispsupR_solution.jl"),
    joinpath("lattices", "chapter_3", "chapter3_mSSF_solution.jl"),
    joinpath("lattices", "chapter_3", "chapter3_mSSR_solution.jl"),
]

for file in solution_files
    solution_path = joinpath(tutorial_root, file)
    isfile(solution_path) || error("Cannot find $(file). Run the preceding chapters first.")
    include(solution_path)
end

const L_quad = 0.5
const D1_len = 0.609
const D2_len = 1.241
const DB_len = 5.855
const B_len = 6.86
const B_angle = pi / 132
const BH_angle = B_angle / 2
const K_ss = 0.351957452649287

const species_ref = Species("electron")
const E_ref = 18e9

# Store all quadrupole strengths under the names used in the original tutorial.
quad_strength = Dict(
    :QF => kQF_arc,
    :QD => kQD_arc,
    :QFR => kQF_arc_R,
    :QDR => kQD_arc_R,
    :QFSS => K_ss,
    :QDSS => -K_ss,
    :QFF1 => kQFF1,
    :QDF1 => kQDF1,
    :QFR1 => kQFR1,
    :QDR1 => kQDR1,
    :QFF2 => K_QFF2,
    :QDF2 => K_QDF2,
    :QFF3 => K_QFF3,
    :QDF3 => K_QDF3,
    :QFR2 => K_QFR2,
    :QDR2 => K_QDR2,
    :QFR3 => K_QFR3,
    :QDR3 => K_QDR3,
)

function make_element(name::Symbol)
    haskey(quad_strength, name) &&
        return Quadrupole(name=String(name), L=L_quad, Kn1=quad_strength[name])
    name == :D1 && return Drift(name="D1", L=D1_len)
    name == :D2 && return Drift(name="D2", L=D2_len)
    name == :DB && return Drift(name="DB", L=DB_len)
    name == :B && return SBend(name="B", L=B_len, angle=B_angle)
    name == :BH && return SBend(name="BH", L=B_len, angle=BH_angle)
    error("Unknown element symbol: $(name)")
end

make_elements(line) = [make_element(name) for name in line]
repeat_line(line, n) = reduce(vcat, (copy(line) for _ in 1:n))

# (2): Gather arc and straight-section FODO cells.

# Straight-section forward FODO.
FODOSSF = [:QFSS, :D1, :DB, :D2, :QDSS, :D1, :DB, :D2]

# Arc forward FODO.
FODOAF = [:QF, :D1, :B, :D2, :QD, :D1, :B, :D2]

# Arc reverse FODO.
FODOAR = [:QFR, :D2, :B, :D1, :QDR, :D2, :B, :D1]

# Straight-section reverse FODO.
FODOSSR = [:QFSS, :D2, :DB, :D1, :QDSS, :D2, :DB, :D1]

# (3): Gather dispersion suppression and matching-to-straight lines.

# Forward dispersion suppressor and forward matching section.
DISPSUPF = [
    :QF, :D1, :BH, :D2, :QD, :D1, :BH, :D2,
    :QFF1, :D1, :BH, :D2, :QDF1, :D1, :BH, :D2,
]
MSSF = [
    :QFF2, :D1, :DB, :D2, :QDF2, :D1, :DB, :D2,
    :QFF3, :D1, :DB, :D2, :QDF3, :D1, :DB, :D2,
]
ARC_TO_SSF = vcat(DISPSUPF, MSSF)

# Reverse dispersion suppressor and reverse matching section.
DISPSUPR = [
    :QFR, :D2, :BH, :D1, :QDR, :D2, :BH, :D1,
    :QFR1, :D2, :BH, :D1, :QDR1, :D2, :BH, :D1,
]
MSSR = [
    :QFR2, :D2, :DB, :D1, :QDR2, :D2, :DB, :D1,
    :QFR3, :D2, :DB, :D1, :QDR3, :D2, :DB, :D1,
]
ARC_TO_SSR = vcat(DISPSUPR, MSSR)

# (4): Gather matching-to-dispersion-creation and dispersion-creation lines.
#
# By mirror symmetry:
#   ARC_TO_SSF <-> SS_TO_ARCR
#   ARC_TO_SSR <-> SS_TO_ARCF

# Match the forward straight section to the forward dispersion creator.
# These quadrupole strengths are copied from the reverse matching section.
MDCF = [
    :QFSS, :D1, :DB, :D2, :QDR3, :D1, :DB, :D2,
    :QFR3, :D1, :DB, :D2, :QDR2, :D1, :DB, :D2,
]
DISPCREF = [
    :QFR2, :D1, :BH, :D2, :QDR1, :D1, :BH, :D2,
    :QFR1, :D1, :BH, :D2, :QDR, :D1, :BH, :D2,
]
SS_TO_ARCF = vcat(MDCF, DISPCREF)

# Match the reverse straight section to the reverse dispersion creator.
# These quadrupole strengths are copied from the forward matching section.
MDCR = [
    :QFSS, :D2, :DB, :D1, :QDF3, :D2, :DB, :D1,
    :QFF3, :D2, :DB, :D1, :QDF2, :D2, :DB, :D1,
]
DISPCRER = [
    :QFF2, :D2, :BH, :D1, :QDF1, :D2, :BH, :D1,
    :QFF1, :D2, :BH, :D1, :QD, :D2, :BH, :D1,
]
SS_TO_ARCR = vcat(MDCR, DISPCRER)

# (5): Build the ring.
SEXTANT1 = vcat(
    repeat_line(FODOSSF, 4), SS_TO_ARCF, repeat_line(FODOAF, 20),
    ARC_TO_SSF, repeat_line(FODOSSF, 4),
)
SEXTANT3 = vcat(
    repeat_line(FODOSSR, 4), SS_TO_ARCR, repeat_line(FODOAR, 20),
    ARC_TO_SSR, repeat_line(FODOSSR, 4),
)
SEXTANT5 = copy(SEXTANT1)
SEXTANT7 = copy(SEXTANT3)
SEXTANT9 = copy(SEXTANT1)
SEXTANT11 = copy(SEXTANT3)

RING_NAMES = vcat(SEXTANT1, SEXTANT3, SEXTANT5, SEXTANT7, SEXTANT9, SEXTANT11)
ring = Beamline(
    make_elements(RING_NAMES);
    species_ref=species_ref,
    E_ref=E_ref,
)

println("Constructed the complete ring:")
@printf("  Elements per forward sextant: %d\n", length(SEXTANT1))
@printf("  Elements per reverse sextant: %d\n", length(SEXTANT3))
@printf("  Elements in full ring:        %d\n", length(RING_NAMES))

# Check that the ring closes geometrically.
element_length(name::Symbol) =
    haskey(quad_strength, name) ? L_quad :
    name == :D1 ? D1_len :
    name == :D2 ? D2_len :
    name == :DB ? DB_len : B_len

element_angle(name::Symbol) = name == :B ? B_angle : name == :BH ? BH_angle : 0.0

struct FloorFrame
    x::Float64
    z::Float64
    theta::Float64
end

function advance_floor(frame::FloorFrame, L, angle)
    if abs(angle) < 1e-14
        return FloorFrame(
            frame.x + L * sin(frame.theta),
            frame.z + L * cos(frame.theta),
            frame.theta,
        )
    end

    radius = L / angle
    theta1 = frame.theta + angle
    return FloorFrame(
        frame.x + radius * (cos(frame.theta) - cos(theta1)),
        frame.z + radius * (sin(theta1) - sin(frame.theta)),
        theta1,
    )
end

function track_floor(line)
    frame = FloorFrame(0.0, 0.0, 0.0)
    for name in line
        frame = advance_floor(frame, element_length(name), element_angle(name))
    end
    return frame
end

frame = track_floor(RING_NAMES)
circumference = sum(element_length(name) for name in RING_NAMES)
total_bend = sum(element_angle(name) for name in RING_NAMES)

println("\nGeometrical closure:")
@printf("  Circumference:       %.9f m\n", circumference)
@printf("  Total bend:          %.15f rad\n", total_bend)
@printf("  Total bend - 2pi:    %+.3e rad\n", total_bend - 2pi)
@printf("  Final x:             %+.3e m\n", frame.x)
@printf("  Final z:             %+.3e m\n", frame.z)
@printf("  Position error:      %.3e m\n", hypot(frame.x, frame.z))

# Calculate the periodic optics of the complete ring.
println("\nCalculating full-ring periodic optics...")
tw_ring = twiss(ring)
println("Full-ring twiss calculation completed.")

twiss_table = hasproperty(tw_ring, :table) ? tw_ring.table : tw_ring
@printf("  Mode-1 tune: %.12f\n", twiss_table.phi_1[end])
@printf("  Mode-2 tune: %.12f\n", twiss_table.phi_2[end])
