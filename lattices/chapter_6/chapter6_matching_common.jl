# Shared linear-optics and matching helpers for Chapter 6.

using SciBmad
using GTPSA
using DifferentiationInterface
import DifferentiationInterface as DI
using LinearAlgebra
using Printf

const species_ref = Species("electron")
const E_ref = 18e9
const K_ss = 0.351957452649287

const D1_len = 0.609
const D2_len = 1.241
const DB_len = 5.855

# Constructors return fresh elements so each trial beamline owns its elements.
D1() = Drift(name="D1", L=D1_len)
D2() = Drift(name="D2", L=D2_len)
DB() = Drift(name="DB", L=DB_len)
IP6() = Marker(name="IP6")

# The forward and reverse straight FODO cells have the same strengths but
# opposite drift ordering. Their periodic Twiss solutions define IR boundaries.
function build_forward_straight_fodo()
    return Beamline(
        [
            Quadrupole(name="QFSS", L=0.5, Kn1=K_ss), D1(), DB(), D2(),
            Quadrupole(name="QDSS", L=0.5, Kn1=-K_ss), D1(), DB(), D2(),
        ];
        species_ref=species_ref,
        E_ref=E_ref,
    )
end

function build_reverse_straight_fodo()
    return Beamline(
        [
            Quadrupole(name="QFSS", L=0.5, Kn1=K_ss), D2(), DB(), D1(),
            Quadrupole(name="QDSS", L=0.5, Kn1=-K_ss), D2(), DB(), D1(),
        ];
        species_ref=species_ref,
        E_ref=E_ref,
    )
end

function track_a_particle(v0, beamline)
    # SciBmad tracking is differentiated with GTPSA to obtain the linear map.
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

function linear_map_with_descriptor(beamline, d; x0=zeros(6))
    # Keep descriptor-parameter dependence in the first-order map coefficients.
    xvars = vars(d)
    v0 = [x0[i] + xvars[i] for i in 1:6]
    vout = track_a_particle(v0, beamline)
    M = Matrix{Any}(undef, 6, 6)
    for i in 1:6, j in 1:6
        powers = zeros(Int, 6)
        powers[j] = 1
        M[i, j] = par(vout[i], [powers...,:])
    end
    return M
end

parameter_gradient(x) = GTPSA.gradient(x, include_params=true)[7:end]
tps_const(x) = try x[zeros(Int, 6)] catch; x end

transverse_blocks(M) = M[1:2, 1:2], M[3:4, 3:4]

function concrete_matrix(M)
    T = foldl(promote_type, typeof.(M); init=Float64)
    return Matrix{T}(M)
end

# Minimal uncoupled Twiss representation used by the Chapter 6 matches.
struct Twiss
    beta
    alpha
end

gamma(t::Twiss) = (1 + t.alpha^2) / t.beta

function sigma_matrix(t::Twiss)
    return [
        t.beta -t.alpha
        -t.alpha gamma(t)
    ]
end

twiss_from_sigma(S) = Twiss(S[1, 1], -S[1, 2])
propagate_twiss(t::Twiss, M) = twiss_from_sigma(M * sigma_matrix(t) * transpose(M))

function periodic_twiss_from_matrix(M)
    # Choose the phase-advance branch that gives a positive periodic beta.
    cos_mu = tr(M) / 2
    abs(cos_mu) < 1 || error("Unstable one-cell matrix: |Tr(M)/2| >= 1")
    sin_mu = sqrt(1 - cos_mu^2)
    M[1, 2] / sin_mu < 0 && (sin_mu = -sin_mu)
    return Twiss(M[1, 2] / sin_mu, (M[1, 1] - M[2, 2]) / (2sin_mu))
end

function periodic_twiss(beamline)
    Mx, My = transverse_blocks(transfer_matrix_gtpsa(beamline))
    return periodic_twiss_from_matrix(Mx), periodic_twiss_from_matrix(My)
end

function residual_jacobian_fd(f, x; fd_step=1e-6)
    # Centered finite differences are kept for Exercise 6, where the step size
    # itself is the quantity being studied.
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

residual_jacobian(f, x; fd_step=1e-6) = residual_jacobian_fd(f, x; fd_step=fd_step)

function damped_least_squares(
    f,
    x0;
    jacobian=nothing,
    fd_step=1e-6,
    maxiter=100,
    tol=1e-12,
    lambda0=1e-3,
    verbose=true,
)
    # Levenberg-Marquardt-style damping stabilizes this sensitive low-beta
    # matching problem when a full Gauss-Newton step is too aggressive.
    x = copy(x0)
    lambda = lambda0

    for iter in 1:maxiter
        r = f(x)
        merit_now = sum(abs2, r)
        J = jacobian === nothing ? residual_jacobian_fd(f, x; fd_step=fd_step) : jacobian(x)
        step = -(J' * J + lambda * I) \ (J' * r)
        trial = x + step
        merit_trial = sum(abs2, f(trial))

        if verbose
            @printf(
                "iter %3d  merit = %.6e  step = %.3e  lambda = %.3e\n",
                iter, merit_now, norm(step), lambda,
            )
        end

        if merit_trial < merit_now
            x = trial
            lambda /= 10
        else
            lambda *= 10
        end

        norm(r) < tol && break
        norm(step) < tol && break
    end
    return x
end

function damped_least_squares_with_history(
    f,
    x0;
    fd_step=1e-6,
    maxiter=100,
    tol=1e-12,
    lambda0=1e-3,
)
    # This variant records every trial so Exercise 6 can compare how the
    # finite-difference step size changes convergence.
    x = copy(x0)
    lambda = lambda0
    history = NamedTuple[]
    converged = false

    for iter in 1:maxiter
        r = f(x)
        merit_now = sum(abs2, r)
        J = residual_jacobian_fd(f, x; fd_step=fd_step)
        step = -(J' * J + lambda * I) \ (J' * r)
        trial = x + step
        merit_trial = sum(abs2, f(trial))
        accepted = merit_trial < merit_now

        push!(
            history,
            (
                iteration=iter,
                merit=merit_now,
                residual_norm=norm(r),
                step_norm=norm(step),
                lambda=lambda,
                accepted=accepted,
            ),
        )

        if accepted
            x = trial
            lambda /= 10
        else
            lambda *= 10
        end

        if norm(r) < tol || norm(step) < tol
            converged = norm(r) < tol
            break
        end
    end

    return x, history, converged
end

const target_a = Twiss(0.6, 0.0)
const target_b = Twiss(0.06, 0.0)

function normalized_ip_residual(output_a, output_b)
    # Normalize beta errors by their targets so the two planes have comparable
    # weight despite their factor-of-ten difference in target beta.
    return [
        (output_a.beta - target_a.beta) / target_a.beta,
        output_a.alpha - target_a.alpha,
        (output_b.beta - target_b.beta) / target_b.beta,
        output_b.alpha - target_b.alpha,
    ]
end
