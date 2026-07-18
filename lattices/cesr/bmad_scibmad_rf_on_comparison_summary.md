# CESR RF-On Bmad–SciBmad Comparison

## Configuration

- Beam energy: 5.289 GeV
- RF cavities: `RF_W1`, `RF_W2`, `RF_E1`, and `RF_E2`
- Voltage: 1.5 MV per cavity (6 MV total)
- RF frequency: 499.7669603784 MHz
- Bmad reference: `bmad_reference_rf_on_output.tar.gz`
- SciBmad lattice: `cesr.jl`, loaded with `load_cesr()` and enabled with
  `set_cesr_rf!(ring; on=true)`

The element comparison linearizes every SciBmad element with GTPSA about the
corresponding Bmad entrance orbit. This assigns discrepancies to the element
that creates them instead of hiding them in a one-turn product.

## Closed Orbit and Tunes

| Coordinate | Bmad | SciBmad | SciBmad - Bmad |
|---|---:|---:|---:|
| x | -1.666982000e-5 | -1.668519860e-5 | -1.537860e-8 |
| px | 2.389568950e-3 | 2.390112503e-3 | 5.435531e-7 |
| y | 1.054880000e-6 | 1.054311223e-6 | -5.687774e-10 |
| py | 1.777040000e-6 | 1.796930342e-6 | 1.989034e-8 |
| z | -3.989291500e-4 | -3.993611670e-4 | -4.320170e-7 |
| pz | -7.163370000e-6 | -7.803721195e-6 | -6.403512e-7 |

The maximum closed-orbit coordinate difference is `6.404e-7` in `pz`.

| Tune | Bmad | SciBmad | SciBmad - Bmad |
|---|---:|---:|---:|
| Qx (fractional) | 0.530014110 | 0.529926212 | -8.790e-5 |
| Qy (fractional) | 0.578460931 | 0.579077486 | 6.166e-4 |
| Qz (eigenphase magnitude) | 0.051738795 | 0.051738677 | -1.182e-7 |

The Bmad transverse tunes are taken from the accumulated Twiss phases printed
by Tao. The longitudinal tune is estimated from the product of Tao's printed
per-element matrices. Those matrices contain only seven printed decimal
places, so their product is less accurate than the Bmad internal one-turn map.
The resulting Bmad eigenvalue moduli differ from one by at most `5.65e-6`.

## Element-by-Element Results

- Maximum local matrix discrepancy: `2.476e-4` at `B06E` (`R44`).
- Maximum cumulative matrix discrepancy: `1.267e-2` at `Q23E`.
- Maximum isolated-element exit-orbit discrepancy: `2.526e-6` at `WIG_E`.
- Maximum element-length discrepancy: `8.338e-7 m`, limited by printed Tao
  longitudinal positions.
- Bmad printed affine-map consistency: `1.629e-9`.

| Element | Maximum local matrix difference | Exit-orbit difference |
|---|---:|---:|
| Q00W/CLEO_SOL overlap | 9.145e-5 | 1.935e-7 |
| WIG_W | 3.951e-5 | 2.526e-6 |
| RF_W1 | 4.857e-8 | 6.999e-12 |
| RF_W2 | 4.857e-8 | 1.571e-11 |
| RF_E1 | 4.245e-8 | 8.176e-12 |
| RF_E2 | 4.191e-8 | 7.410e-12 |
| WIG_E | 3.775e-5 | 2.526e-6 |
| Q00E/CLEO_SOL overlap | 9.141e-5 | 4.189e-7 |

The RF cavity maps agree to approximately `5e-8`, so RF translation is not a
significant source of the present Bmad–SciBmad difference. The dominant local
matrix discrepancies remain in sector bends, followed by the combined
solenoid–quadrupole overlaps and the two wigglers. The largest affine/orbit
difference remains associated with the wiggler model.

## RF-On Versus RF-Off

| Metric | RF off | RF on |
|---|---:|---:|
| Maximum local matrix difference | 2.478e-4 | 2.476e-4 |
| Maximum cumulative matrix difference | 1.258e-2 | 1.267e-2 |
| Maximum isolated-element orbit difference | 2.526e-6 | 2.526e-6 |

Turning on the RF changes the maximum cumulative matrix discrepancy by less
than one percent and does not change the identity of the dominant local error
sources. It does, however, make the closed-orbit problem six-dimensional and
avoids the four-variable ForwardDiff failure encountered by the default
coasting-beam closed-orbit calculation.

## Reproduction

```bash
julia --project=. lattices/cesr/test_bmad_scibmad.jl \
  --rf-on \
  --reference=lattices/cesr/bmad_reference_rf_on_output.tar.gz \
  --csv=lattices/cesr/bmad_scibmad_rf_on_comparison.csv

julia --project=. lattices/cesr/compare_rf_on_optics.jl
```
