# Chapter 6, Exercise 6:
# Compare IPF matching with finite-difference quadrupole-strength step sizes
# of 1e-6 and 1e-4.

include(joinpath(@__DIR__, "..", "chapter6_matching_common.jl"))

function build_IPF_for_step_test(k)
    return Beamline(
        [
            Quadrupole(name="QEF1", L=0.5, Kn1=k[1]), D1(), DB(), D2(),
            Quadrupole(name="QEF2", L=0.5, Kn1=k[2]), D1(), DB(), D2(),
            Drift(name="DEF1", L=20.46),
            Quadrupole(name="QEF3", L=1.6, Kn1=k[3]),
            Drift(name="DEF2", L=3.76),
            Quadrupole(name="QEF4", L=1.2, Kn1=k[4]),
            Drift(name="DEF3", L=5.8),
            IP6(),
        ];
        species_ref=species_ref,
        E_ref=E_ref,
    )
end

input_a_step_test, input_b_step_test = periodic_twiss(build_forward_straight_fodo())

function ipf_step_test_residual(k)
    Mx, My = transverse_blocks(transfer_matrix_gtpsa(build_IPF_for_step_test(k)))
    output_a = propagate_twiss(input_a_step_test, Mx)
    output_b = propagate_twiss(input_b_step_test, My)
    return normalized_ip_residual(output_a, output_b)
end

function run_step_test(fd_step)
    K_start = [K_ss, -K_ss, K_ss, -K_ss]
    k, history, converged = damped_least_squares_with_history(
        ipf_step_test_residual,
        K_start;
        fd_step=fd_step,
        maxiter=100,
        tol=1e-12,
    )

    final_residual = ipf_step_test_residual(k)
    final_merit = sum(abs2, final_residual)
    accepted_steps = count(row -> row.accepted, history)

    return (
        fd_step=fd_step,
        strengths=k,
        residual=final_residual,
        merit=final_merit,
        iterations=length(history),
        accepted_steps=accepted_steps,
        converged=converged,
        history=history,
    )
end

result_1e6 = run_step_test(1e-6)
result_1e4 = run_step_test(1e-4)

# Compare the actual Jacobians used by the optimizer. Centered differences have
# second-order truncation error, so 1e-4 can still be sufficiently accurate for
# this SciBmad implementation.
K_start_comparison = [K_ss, -K_ss, K_ss, -K_ss]
J_start_1e6 = residual_jacobian(ipf_step_test_residual, K_start_comparison; fd_step=1e-6)
J_start_1e4 = residual_jacobian(ipf_step_test_residual, K_start_comparison; fd_step=1e-4)
relative_start_jacobian_difference = norm(J_start_1e4 - J_start_1e6) / norm(J_start_1e6)

J_final_1e6 = residual_jacobian(ipf_step_test_residual, result_1e6.strengths; fd_step=1e-6)
J_final_1e4 = residual_jacobian(ipf_step_test_residual, result_1e6.strengths; fd_step=1e-4)
relative_final_jacobian_difference = norm(J_final_1e4 - J_final_1e6) / norm(J_final_1e6)

function print_result(result)
    println("\nFinite-difference step = ", result.fd_step)
    println("  converged      = ", result.converged)
    println("  iterations     = ", result.iterations)
    println("  accepted steps = ", result.accepted_steps)
    @printf("  final merit    = %.12e\n", result.merit)
    println("  final residual = ", result.residual)
    println("  strengths      = ", result.strengths)

    println("  selected history:")
    for row in result.history
        if row.iteration <= 5 || row.iteration % 10 == 0 || row.iteration == result.iterations
            @printf(
                "    iter %3d  merit = %.6e  step = %.3e  accepted = %s\n",
                row.iteration,
                row.merit,
                row.step_norm,
                string(row.accepted),
            )
        end
    end
end

print_result(result_1e6)
print_result(result_1e4)

println("\nComparison:")
@printf("  merit ratio (1e-4 / 1e-6) = %.6e\n", result_1e4.merit / result_1e6.merit)
@printf(
    "  strength difference norm   = %.6e\n",
    norm(result_1e4.strengths - result_1e6.strengths),
)
@printf(
    "  initial Jacobian relative difference = %.6e\n",
    relative_start_jacobian_difference,
)
@printf(
    "  final Jacobian relative difference   = %.6e\n",
    relative_final_jacobian_difference,
)
