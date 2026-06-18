# Reverse Chapter 3: match the reverse suppressor to the reverse straight FODO.

using SciBmad
using GTPSA
using DifferentiationInterface
import DifferentiationInterface as DI
using LinearAlgebra
using Printf

const tutorial_root = normpath(joinpath(@__DIR__, "..", ".."))
const chapter2_solution = joinpath(tutorial_root, "lattices", "chapter_2", "chapter2_dispsupR_solution.jl")

isfile(chapter2_solution) ||
    error("Run reverse_chapter2_dispersion_suppressor.jl first.")
include(chapter2_solution)

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

D1() = Drift(L=D1_len)
D2() = Drift(L=D2_len)
DB() = Drift(L=DB_len)
BH() = SBend(L=B_len, angle=BH_angle)

QF_arc_R() = Quadrupole(L=L_quad, Kn1=kQF_arc_R)
QD_arc_R() = Quadrupole(L=L_quad, Kn1=kQD_arc_R)
QFR1() = Quadrupole(L=L_quad, Kn1=kQFR1)
QDR1() = Quadrupole(L=L_quad, Kn1=kQDR1)
QFSS() = Quadrupole(L=L_quad, Kn1=K_ss)
QDSS() = Quadrupole(L=L_quad, Kn1=-K_ss)

function build_reverse_suppressor()
    return Beamline(
        [
            QF_arc_R(), D2(), BH(), D1(),
            QD_arc_R(), D2(), BH(), D1(),
            QFR1(), D2(), BH(), D1(),
            QDR1(), D2(), BH(), D1(),
        ];
        species_ref=species_ref,
        E_ref=E_ref,
    )
end

function build_reverse_straight_fodo()
    return Beamline(
        [
            QFSS(), D2(), DB(), D1(),
            QDSS(), D2(), DB(), D1(),
        ];
        species_ref=species_ref,
        E_ref=E_ref,
    )
end

function build_MSSR(k)
    return Beamline(
        [
            Quadrupole(L=L_quad, Kn1=k[1]), D2(), DB(), D1(),
            Quadrupole(L=L_quad, Kn1=k[2]), D2(), DB(), D1(),
            Quadrupole(L=L_quad, Kn1=k[3]), D2(), DB(), D1(),
            Quadrupole(L=L_quad, Kn1=k[4]), D2(), DB(), D1(),
        ];
        species_ref=species_ref,
        E_ref=E_ref,
    )
end

function track_a_particle(v0, beamline)
    v = similar(v0)
    v .= v0
    bunch = Bunch(v; species=beamline.species_ref, p_over_q_ref=beamline.p_over_q_ref)
    track!(bunch, beamline)
    return copy(bunch.coords.v)
end

function transfer_matrix_gtpsa(beamline; x0=zeros(6))
    prep = DI.prepare_jacobian(
        track_a_particle,
        AutoGTPSA(),
        x0,
        DI.Constant(beamline),
    )
    return DI.jacobian(
        track_a_particle,
        prep,
        AutoGTPSA(),
        x0,
        DI.Constant(beamline),
    )
end

function transverse_blocks(M)
    return M[1:2, 1:2], M[3:4, 3:4]
end

struct Twiss
    beta::Float64
    alpha::Float64
end

gamma(t::Twiss) = (1 + t.alpha^2) / t.beta

function sigma_matrix(t::Twiss)
    return [
        t.beta -t.alpha
        -t.alpha gamma(t)
    ]
end

function twiss_from_sigma(S)
    return Twiss(S[1, 1], -S[1, 2])
end

function propagate_twiss(t::Twiss, M)
    return twiss_from_sigma(M * sigma_matrix(t) * transpose(M))
end

function periodic_twiss_from_matrix(M)
    cos_mu = tr(M) / 2
    abs(cos_mu) < 1 || error("Unstable one-cell matrix: |Tr(M)/2| >= 1")
    sin_mu = sqrt(1 - cos_mu^2)
    M[1, 2] / sin_mu < 0 && (sin_mu = -sin_mu)
    return Twiss(M[1, 2] / sin_mu, (M[1, 1] - M[2, 2]) / (2sin_mu))
end

M_suppressor = transfer_matrix_gtpsa(build_reverse_suppressor())
Mx_suppressor, My_suppressor = transverse_blocks(M_suppressor)

# Start with the periodic reverse-arc Twiss at the suppressor entrance.
function build_reverse_arc_fodo()
    return Beamline(
        [
            QF_arc_R(), D2(), SBend(L=B_len, angle=B_angle), D1(),
            QD_arc_R(), D2(), SBend(L=B_len, angle=B_angle), D1(),
        ];
        species_ref=species_ref,
        E_ref=E_ref,
    )
end

Mx_arc, My_arc = transverse_blocks(transfer_matrix_gtpsa(build_reverse_arc_fodo()))
arc_x = periodic_twiss_from_matrix(Mx_arc)
arc_y = periodic_twiss_from_matrix(My_arc)
input_x = propagate_twiss(arc_x, Mx_suppressor)
input_y = propagate_twiss(arc_y, My_suppressor)

# The target is the periodic Twiss solution of the reverse straight FODO.
Mx_ss, My_ss = transverse_blocks(transfer_matrix_gtpsa(build_reverse_straight_fodo()))
target_x = periodic_twiss_from_matrix(Mx_ss)
target_y = periodic_twiss_from_matrix(My_ss)

function matching_residual(k)
    Mx, My = transverse_blocks(transfer_matrix_gtpsa(build_MSSR(k)))
    output_x = propagate_twiss(input_x, Mx)
    output_y = propagate_twiss(input_y, My)
    return [
        (output_x.beta - target_x.beta) / target_x.beta,
        output_x.alpha - target_x.alpha,
        (output_y.beta - target_y.beta) / target_y.beta,
        output_y.alpha - target_y.alpha,
    ]
end

function residual_jacobian(f, x; fd_step=1e-6)
    r0 = f(x)
    J = zeros(length(r0), length(x))
    for j in eachindex(x)
        xp = copy(x)
        xm = copy(x)
        xp[j] += fd_step
        xm[j] -= fd_step
        J[:, j] = (f(xp) - f(xm)) / (2fd_step)
    end
    return J
end

function damped_least_squares(f, x0; maxiter=60, tol=1e-12, lambda0=1e-3)
    x = copy(x0)
    lambda = lambda0

    for iter in 1:maxiter
        r = f(x)
        merit_now = sum(abs2, r)
        J = residual_jacobian(f, x)
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

        if norm(r) < tol || norm(step) < tol
            break
        end
    end
    return x
end

K_start = [K_ss, -K_ss, K_ss, -K_ss]
K_match_R = damped_least_squares(matching_residual, K_start)
K_QFR2, K_QDR2, K_QFR3, K_QDR3 = K_match_R

println("\nReverse straight FODO target Twiss:")
@printf("  beta_x  = %.12f, alpha_x = %.12f\n", target_x.beta, target_x.alpha)
@printf("  beta_y  = %.12f, alpha_y = %.12f\n", target_y.beta, target_y.alpha)

println("\nOptimized reverse matching section:")
@printf("  K_QFR2 = %.15f\n", K_QFR2)
@printf("  K_QDR2 = %.15f\n", K_QDR2)
@printf("  K_QFR3 = %.15f\n", K_QFR3)
@printf("  K_QDR3 = %.15f\n", K_QDR3)
println("  final residual = ", matching_residual(K_match_R))

solution_text = """
# chapter3_mSSR_solution.jl
# Auto-generated by reverse_chapter3_twiss_matching.jl.

K_QFR2 = $(repr(K_QFR2))
K_QDR2 = $(repr(K_QDR2))
K_QFR3 = $(repr(K_QFR3))
K_QDR3 = $(repr(K_QDR3))
"""

solution_path = joinpath(tutorial_root, "lattices", "chapter_3", "chapter3_mSSR_solution.jl")
mkpath(dirname(solution_path))
write(solution_path, solution_text)
println("\nWrote ", solution_path)
