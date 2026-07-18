#!/usr/bin/env julia

using SciBmad
using GTPSA: Descriptor

include("cesr.jl")

ring = load_cesr()
set_cesr_rf!(ring; on=true)

println("Enabled RF cavities: ", join(sort!(collect(CESR_RF_NAMES)), ", "))
println("Voltage per cavity: ", CESR_RF_VOLTAGE, " V")
println("Computing the six-dimensional closed orbit...")
closed_orbit = find_closed_orbit(ring)
println("Closed orbit: ", closed_orbit)

println("Computing Twiss parameters...")
desc = Descriptor(6, 1)
optics = twiss(
  ring;
  GTPSA_descriptor=desc,
  v0=closed_orbit.v0,
  v0_and_coast=(closed_orbit.v0, closed_orbit.coasting_beam),
)
println("Twiss calculation completed successfully.")
println("Tunes [Qx, Qy, Qz]: ", optics.tunes)
println(optics)
