# Generated from chapter15_orbit_correction_scibmad.ipynb.
# This result uses SciBmad deterministic radiation tracking:
# Yoshida(order=2, radiation_damping_on=true, radiation_fluctuations_on=false).
#
# The taper scale was applied to SBend.Kn0/g_ref, Quadrupole.Kn1, and nonzero
# Sextupole.Kn2. Corrector kickers were deliberately left unchanged.

const CH15_REAL_TAPERING_RESULT = (
    source = "chapter15_b_sawtooth_ring0_beamline.jl",
    tapered_strengths = (
        total = 956,
        sbends = 288,
        quadrupoles = 428,
        sextupoles = 240,
    ),
    pz_before = (
        peak_to_peak = 1.665154e-3,
    ),
    pz_after = (
        peak_to_peak = 1.665157e-3,
    ),
    orbit_x_before = (
        rms_mm = 0.175936,
        peak_to_peak_mm = 0.988046,
    ),
    orbit_x_after = (
        rms_mm = 0.176016,
        peak_to_peak_mm = 0.991835,
    ),
    taper_scale_range = (
        min = 0.999164972,
        max = 1.000830120,
    ),
    largest_relative_strength_change = 8.350e-4,
)
