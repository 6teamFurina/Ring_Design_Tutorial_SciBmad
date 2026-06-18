# Reusable Chapter 5 ring definition.
#
# Include this file from later chapters to rebuild the unmodified Chapter 5
# ring without running the Chapter 5 exercise checks or full-ring Twiss solve.

module Chapter5Ring

using SciBmad

tutorial_root = normpath(joinpath(@__DIR__, "..", ".."))

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

L_quad = 0.5
D1_len = 0.609
D2_len = 1.241
DB_len = 5.855
B_len = 6.86
B_angle = pi / 132
BH_angle = B_angle / 2
K_ss = 0.351957452649287

species_ref = Species("electron")
E_ref = 18e9

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

FODOAF = [:QF, :D1, :B, :D2, :QD, :D1, :B, :D2]
FODOAR = [:QFR, :D2, :B, :D1, :QDR, :D2, :B, :D1]
FODOSSF = [:QFSS, :D1, :DB, :D2, :QDSS, :D1, :DB, :D2]
FODOSSR = [:QFSS, :D2, :DB, :D1, :QDSS, :D2, :DB, :D1]

ARC_TO_SSF = [
    :QF, :D1, :BH, :D2, :QD, :D1, :BH, :D2,
    :QFF1, :D1, :BH, :D2, :QDF1, :D1, :BH, :D2,
    :QFF2, :D1, :DB, :D2, :QDF2, :D1, :DB, :D2,
    :QFF3, :D1, :DB, :D2, :QDF3, :D1, :DB, :D2,
]

SS_TO_ARCF = [
    :QFSS, :D1, :DB, :D2, :QDR3, :D1, :DB, :D2,
    :QFR3, :D1, :DB, :D2, :QDR2, :D1, :DB, :D2,
    :QFR2, :D1, :BH, :D2, :QDR1, :D1, :BH, :D2,
    :QFR1, :D1, :BH, :D2, :QDR, :D1, :BH, :D2,
]

ARC_TO_SSR = [
    :QFR, :D2, :BH, :D1, :QDR, :D2, :BH, :D1,
    :QFR1, :D2, :BH, :D1, :QDR1, :D2, :BH, :D1,
    :QFR2, :D2, :DB, :D1, :QDR2, :D2, :DB, :D1,
    :QFR3, :D2, :DB, :D1, :QDR3, :D2, :DB, :D1,
]

SS_TO_ARCR = [
    :QFSS, :D2, :DB, :D1, :QDF3, :D2, :DB, :D1,
    :QFF3, :D2, :DB, :D1, :QDF2, :D2, :DB, :D1,
    :QFF2, :D2, :BH, :D1, :QDF1, :D2, :BH, :D1,
    :QFF1, :D2, :BH, :D1, :QD, :D2, :BH, :D1,
]

SEXTANT1 = vcat(
    repeat_line(FODOSSF, 4),
    SS_TO_ARCF,
    repeat_line(FODOAF, 20),
    ARC_TO_SSF,
    repeat_line(FODOSSF, 4),
)

SEXTANT3 = vcat(
    repeat_line(FODOSSR, 4),
    SS_TO_ARCR,
    repeat_line(FODOAR, 20),
    ARC_TO_SSR,
    repeat_line(FODOSSR, 4),
)

SEXTANT5 = copy(SEXTANT1)
SEXTANT7 = copy(SEXTANT3)
SEXTANT9 = copy(SEXTANT1)
SEXTANT11 = copy(SEXTANT3)

RING_NAMES = vcat(SEXTANT1, SEXTANT3, SEXTANT5, SEXTANT7, SEXTANT9, SEXTANT11)
ring = Beamline(make_elements(RING_NAMES); species_ref=species_ref, E_ref=E_ref)

end
