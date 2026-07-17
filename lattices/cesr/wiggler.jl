module WigglerModels

using SciBmad

export PlanarWiggler,
       Wiggler,
       planar_wiggler_field,
       planar_wiggler_four_potential,
       planar_wiggler_scales,
       gtpsa_transport_map,
       gtpsa_linear_map

"""
    planar_wiggler_four_potential(x, y, s, t, params)

Return the magnetostatic four-potential and its 4-by-4 Jacobian for the
`kx = 0` specialization of Bmad's periodic planar-wiggler model.

`params` is `(B_max, k_w, phase)`, where `B_max` is in tesla, `k_w` is in
`m^-1`, and `phase` is in radians. The potential is returned in physical
units, so the corresponding element must set
`four_potential_normalized = false`.

The gauge used here is

    phi = Ay = As = 0
    Ax = (B_max / k_w) cosh(k_w y) sin(k_w s + phase).

The derivative tuple follows Beamlines' order:

    (dphi/dx, dphi/dy, dphi/ds, dphi/dt,
     dAx/dx,  dAx/dy,  dAx/ds,  dAx/dt,
     dAy/dx,  dAy/dy,  dAy/ds,  dAy/dt,
     dAs/dx,  dAs/dy,  dAs/ds,  dAs/dt).
"""
@inline function planar_wiggler_four_potential(x, y, s, t, params)
    B_max, k_w, phase = params
    theta = k_w * s + phase

    A_x = (B_max / k_w) * cosh(k_w * y) * sin(theta)
    dA_x_dy = B_max * sinh(k_w * y) * sin(theta)
    dA_x_ds = B_max * cosh(k_w * y) * cos(theta)
    z = zero(A_x)

    potential = (z, A_x, z, z)
    derivatives = (
        z, z, z, z,
        z, dA_x_dy, dA_x_ds, z,
        z, z, z, z,
        z, z, z, z,
    )
    return potential, derivatives
end

"""
    planar_wiggler_field(x, y, s, params)

Return `(B_x, B_y, B_s)` in tesla for the same planar-wiggler field used by
`planar_wiggler_four_potential`. This helper is useful for validation and
plotting; tracking itself uses the four-potential.
"""
@inline function planar_wiggler_field(x, y, s, params)
    B_max, k_w, phase = params
    theta = k_w * s + phase
    z = zero(B_max * cos(theta))
    B_y = B_max * cosh(k_w * y) * cos(theta)
    B_s = -B_max * sinh(k_w * y) * sin(theta)
    return (z, B_y, B_s)
end

"""
    planar_wiggler_scales(; B_max, L_period, p0c)

Return useful small-amplitude scales for a relativistic unit-charge particle.
`p0c` is in eV, `B_max` in tesla, and `L_period` in metres. The returned
named tuple contains magnetic rigidity, peak curvature radius, wavenumber,
horizontal oscillation amplitude and angle, and average vertical focusing.
"""
function planar_wiggler_scales(; B_max, L_period, p0c)
    B_max > 0 || throw(ArgumentError("B_max must be positive"))
    L_period > 0 || throw(ArgumentError("L_period must be positive"))
    p0c > 0 || throw(ArgumentError("p0c must be positive"))

    # B*rho [T m] = p [GeV/c] / 0.299792458 for |q/e| = 1.
    B_rho = (p0c / 1e9) / 0.299792458
    rho_w = B_rho / B_max
    k_w = 2pi / L_period
    orbit_amplitude = inv(rho_w * k_w^2)
    angle_amplitude = inv(rho_w * k_w)
    K_y_average = inv(2rho_w^2)

    return (;
        B_rho,
        rho_w,
        k_w,
        orbit_amplitude,
        angle_amplitude,
        K_y_average,
    )
end

"""
    PlanarWiggler(; B_max, L_period, N_period, L=nothing, ...)

Construct a Beamlines `LineElement` whose field is Bmad's periodic planar
model with `k_x = 0`. The field is integrated by SciBmad's implicit Yoshida
integrator, so ordinary tracking, automatic differentiation, GTPSA maps,
spin, and radiation options all use the same physical field definition.

The default phase, `-k_w*L/2`, is Bmad's convention: the on-axis vertical
field is symmetric about the element centre. `slices_per_period` controls
accuracy; 16 is a practical starting value for CESR optics studies.

This is a continuous-field model. It is not expected to reproduce every term
of Bmad's faster `bmad_standard` planar-wiggler approximation exactly; compare
its map against both Bmad field integration and the CESR baseline.
"""
function PlanarWiggler(;
    B_max,
    L_period,
    N_period::Integer,
    L=nothing,
    phase=nothing,
    slices_per_period::Integer=16,
    order::Integer=6,
    radiation_damping_on::Bool=false,
    radiation_fluctuations_on::Bool=false,
    kwargs...,
)
    B_max >= 0 || throw(ArgumentError("B_max must be nonnegative"))
    L_period > 0 || throw(ArgumentError("L_period must be positive"))
    N_period > 0 || throw(ArgumentError("N_period must be positive"))
    slices_per_period > 0 || throw(ArgumentError("slices_per_period must be positive"))
    order in (2, 4, 6, 8) ||
        throw(ArgumentError("Yoshida order must be one of 2, 4, 6, or 8"))

    expected_L = N_period * L_period
    resolved_L = isnothing(L) ? expected_L : L
    resolved_L > 0 || throw(ArgumentError("L must be positive"))
    isapprox(resolved_L, expected_L; rtol=1e-10, atol=1e-12) ||
        throw(ArgumentError("L=$resolved_L is inconsistent with N_period*L_period=$expected_L"))

    k_w = 2pi / L_period
    resolved_phase = isnothing(phase) ? -k_w * resolved_L / 2 : phase
    n_steps = N_period * slices_per_period

    return LineElement(;
        kind="Wiggler",
        L=resolved_L,
        four_potential=planar_wiggler_four_potential,
        four_potential_params=(B_max, k_w, resolved_phase),
        four_potential_normalized=false,
        tracking_method=Yoshida(;
            order,
            n_steps,
            radiation_damping_on,
            radiation_fluctuations_on,
        ),
        kwargs...,
    )
end

# A concise alias that matches Bmad's element name while retaining explicit,
# reusable physical parameters.
Wiggler(; kwargs...) = PlanarWiggler(; kwargs...)

"""
    gtpsa_transport_map(beamline; v0=zeros(6), order=1)

Track a six-dimensional GTPSA identity map through `beamline` about `v0` and
return the resulting `DAMap`. No finite differencing is used: every Taylor
coefficient is propagated by the same element tracking kernels used for
ordinary particles.

`order=1` is sufficient for the linear transport matrix. Increase `order` to
retain nonlinear map terms. The explicit descriptor is deliberately kept
current while tracking because BeamTracking's implicit four-potential kernel
allocates temporary dynamic TPS objects from `GTPSA.desc_current`.
"""
function gtpsa_transport_map(beamline; v0=zeros(6), order::Integer=1)
    length(v0) == 6 || throw(ArgumentError("v0 must contain six canonical coordinates"))
    order > 0 || throw(ArgumentError("GTPSA order must be positive"))

    descriptor = Descriptor(6, order)
    delta = vars(descriptor)
    reference = collect(v0)
    coordinates = reshape(
        [reference[i] + copy(delta[i]) for i in eachindex(reference)],
        1,
        6,
    )
    bunch = Bunch(v=coordinates)
    SciBmad.BTBL.check_bl_bunch!(bunch, beamline, false)
    track!(bunch, beamline)

    return DAMap(v0=reference, v=vec(bunch.coords.v))
end

"""
    gtpsa_linear_map(beamline; v0=zeros(6))

Return the exact first-order Jacobian propagated by GTPSA through `beamline`.
"""
gtpsa_linear_map(beamline; v0=zeros(6)) =
    jacobian(gtpsa_transport_map(beamline; v0, order=1).v)

end # module WigglerModels
