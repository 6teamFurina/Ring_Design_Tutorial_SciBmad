#!/usr/bin/env julia

using LinearAlgebra
using Printf

include("test_bmad_scibmad.jl")

const RF_ON_REFERENCE = joinpath(HERE, "bmad_reference_rf_on_output.tar.gz")

function main()
  text = read_reference(RF_ON_REFERENCE)
  reference = parse_bmad_reference(text)
  validate_reference(reference)

  R_bmad = Matrix{Float64}(I, 6, 6)
  for i in 1:869
    R_bmad = reference[i].R * R_bmad
  end

  bmad_eigenvalues = eigvals(R_bmad)
  bmad_mode_tunes = sort(
    abs.(angle.(filter(value -> imag(value) > 0, bmad_eigenvalues))) ./ (2pi),
  )

  phi_matches = collect(eachmatch(
    Regex("(?m)^\\s*Phi \\(rad\\)\\s+(" * NUMBER_RE.pattern * ")\\s+(" * NUMBER_RE.pattern * ")"),
    text,
  ))
  isempty(phi_matches) && error("No Bmad Twiss phase rows found")
  final_phi = parse_number.(phi_matches[end].captures)
  bmad_transverse_tunes = mod.(final_phi ./ (2pi), 1)

  ring = load_cesr()
  set_cesr_rf!(ring; on=true)
  closed_orbit = find_closed_orbit(ring)
  descriptor = GTPSA.Descriptor(6, 1)
  optics = twiss(
    ring;
    GTPSA_descriptor=descriptor,
    v0=closed_orbit.v0,
    v0_and_coast=(closed_orbit.v0, closed_orbit.coasting_beam),
  )
  scibmad_tunes = GTPSA.scalar.(optics.tunes)
  bmad_orbit = reference[0].orbit_out
  scibmad_orbit = vec(closed_orbit.v0)

  println("Bmad RF-on closed orbit:    ", bmad_orbit)
  println("SciBmad RF-on closed orbit: ", scibmad_orbit)
  println("SciBmad - Bmad orbit:       ", scibmad_orbit - bmad_orbit)
  @printf("Maximum closed-orbit difference: %.9e\n", maximum(abs, scibmad_orbit - bmad_orbit))
  println("Bmad transverse fractional tunes from Twiss phase: ", bmad_transverse_tunes)
  println("Bmad eigenphase magnitudes:                         ", bmad_mode_tunes)
  println("SciBmad signed tunes:                               ", scibmad_tunes)
  println("SciBmad eigenphase magnitudes:                      ", sort(abs.(scibmad_tunes)))
  println("Bmad one-turn eigenvalue moduli:                    ", sort(abs.(bmad_eigenvalues)))
end

main()
