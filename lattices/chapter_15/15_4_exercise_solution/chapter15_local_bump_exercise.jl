using CairoMakie
using Printf
using Statistics

include(joinpath(@__DIR__, "..", "chapter15_orbit_correction_tools.jl"))
include(joinpath(@__DIR__, "..", "15_1_before_optimization", "chapter15_sawtooth_ring_before.jl"))
include(joinpath(@__DIR__, "..", "15_2_after_optimization", "chapter15_sawtooth_ring_after.jl"))

ring0_layout = CH14_RING0_LAYOUT
bpm_s = ring0_layout.bpm_s
ch_s = ring0_layout.ch_s
cv_s = ring0_layout.cv_s
ch_names_14 = ring0_layout.ch_names
cv_names_14 = ring0_layout.cv_names
x_corrected = CH14_SAWTOOTH_RING_AFTER.x_bpm

ip_s = ring0_layout.circumference / 2
bump_close_s = min(ch_s[108], cv_s[108]) - 0.5

horizontal_bump_names = ["CH_104", "CH_105", "CH_106", "CH_107"]
vertical_bump_names = ["CV_104", "CV_105", "CV_106", "CV_107"]

horizontal_bump_indices = [findfirst(==(name), ch_names_14) for name in horizontal_bump_names]
vertical_bump_indices = [findfirst(==(name), cv_names_14) for name in vertical_bump_names]

horizontal_bump_s = ch_s[horizontal_bump_indices]
vertical_bump_s = cv_s[vertical_bump_indices]

horizontal_initial_ip = [
    chapter15_linear_interpolate(bpm_s, x_corrected, ip_s),
    chapter15_local_slope(bpm_s, x_corrected, ip_s),
]

# Representative measured vertical IP error for this exercise.
vertical_initial_ip = [-0.42e-3, 5.5e-6]

horizontal_local_bump = chapter15_local_bump_solution(
    horizontal_bump_s,
    ip_s,
    bump_close_s,
    horizontal_initial_ip;
    circumference=ring0_layout.circumference,
    tune=54.23,
    plane=:x,
)

vertical_local_bump = chapter15_local_bump_solution(
    vertical_bump_s,
    ip_s,
    bump_close_s,
    vertical_initial_ip;
    circumference=ring0_layout.circumference,
    tune=54.31,
    plane=:y,
)

println("Chapter 15.4 local IP bump")
@printf("IP s position: %.3f m\n", ip_s)
@printf("Bump closure point: %.3f m\n", bump_close_s)
println("Horizontal bump kickers: ", join(horizontal_bump_names, ", "))
println("Vertical bump kickers:   ", join(vertical_bump_names, ", "))
@printf("Initial IP x:  %.4f mm,  px: %.4f urad\n", 1e3 * horizontal_initial_ip[1], 1e6 * horizontal_initial_ip[2])
@printf("Final IP x:    %.4e mm, px: %.4e urad\n", 1e3 * horizontal_local_bump.final_ip_phase_space[1], 1e6 * horizontal_local_bump.final_ip_phase_space[2])
@printf("Initial IP y:  %.4f mm,  py: %.4f urad\n", 1e3 * vertical_initial_ip[1], 1e6 * vertical_initial_ip[2])
@printf("Final IP y:    %.4e mm, py: %.4e urad\n", 1e3 * vertical_local_bump.final_ip_phase_space[1], 1e6 * vertical_local_bump.final_ip_phase_space[2])
@printf("Horizontal closure residual: %.3e m, %.3e rad\n", horizontal_local_bump.closure_change[1], horizontal_local_bump.closure_change[2])
@printf("Vertical closure residual:   %.3e m, %.3e rad\n", vertical_local_bump.closure_change[1], vertical_local_bump.closure_change[2])

println("Horizontal local bump kicks:")
for (name, kick) in zip(horizontal_bump_names, horizontal_local_bump.kicks)
    @printf("  %-6s %+10.4f urad\n", name, 1e6 * kick)
end

println("Vertical local bump kicks:")
for (name, kick) in zip(vertical_bump_names, vertical_local_bump.kicks)
    @printf("  %-6s %+10.4f urad\n", name, 1e6 * kick)
end

local_bump_solution = (
    ip_s = ip_s,
    close_s = bump_close_s,
    horizontal_kicker_names = horizontal_bump_names,
    vertical_kicker_names = vertical_bump_names,
    horizontal_kicker_s = horizontal_bump_s,
    vertical_kicker_s = vertical_bump_s,
    horizontal_kicks = horizontal_local_bump.kicks,
    vertical_kicks = vertical_local_bump.kicks,
    horizontal_initial_ip = horizontal_initial_ip,
    vertical_initial_ip = vertical_initial_ip,
    horizontal_final_ip = horizontal_local_bump.final_ip_phase_space,
    vertical_final_ip = vertical_local_bump.final_ip_phase_space,
    horizontal_closure_change = horizontal_local_bump.closure_change,
    vertical_closure_change = vertical_local_bump.closure_change,
    horizontal_residual = horizontal_local_bump.residual,
    vertical_residual = vertical_local_bump.residual,
)

local_bump_solution_path = joinpath(@__DIR__, "chapter15_local_bump_solution.jl")
chapter15_write_local_bump_solution(local_bump_solution_path, local_bump_solution)

@printf("Saved local bump solution: %s\n", local_bump_solution_path)
