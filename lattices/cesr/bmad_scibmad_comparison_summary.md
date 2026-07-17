# CESR Bmad--SciBmad comparison

Run on 2026-07-17 with `bmad_reference_output.tar.gz` and the consolidated
`cesr.jl` model.

## Method

- The Tao fallback export contains tracking indices 0 through 869 and lord
  indices 870 through 872. Only tracking indices are compared.
- For element `i`, the SciBmad element is linearized with first-order GTPSA
  about the Bmad orbit at the exit of element `i-1`.
- The two Q00--CLEO solenoid overlap regions use the Bmad integration setting:
  second-order `SolenoidKick` with 47 steps.
- The local matrix, affine kick, exit orbit, and cumulative matrix product are
  compared separately. No finite differences are used.
- Tao printed the matrix entries with seven digits after the decimal point, so
  differences of approximately `5e-8` represent agreement at export precision.

## Result

- All 869 tracking elements were compared.
- Bmad's printed affine map reproduces its printed exit orbit to `1.474e-9` or
  better.
- The maximum element-length discrepancy is `8.338e-7 m`, consistent with the
  six-decimal longitudinal positions in the Tao output.
- Element 5, `Q00W\CLEO_SOL`, is Bmad's first `Sol_Quad` overlap region. After
  matching Bmad's 47-step integration, its maximum local matrix difference is
  `9.142e-5` and its exit-orbit difference is `1.935e-7`. Before this correction
  the corresponding differences were `1.071e-1` and `3.409e-4`.
- The matching east-side combined element, index 865 `Q00E\CLEO_SOL`, now has
  a maximum local matrix difference of `9.139e-5`.
- The two wiggler discrepancies are `3.944e-5` (`WIG_W`) and `3.789e-5`
  (`WIG_E`), with exit-orbit discrepancies of about `2.526e-6`.
- The largest cumulative matrix discrepancy is now `1.258e-2` at index 674
  `Q23E`, down from `8.076` before correcting the overlap integration. This is
  an accumulated consequence and should not be interpreted as a local error in
  `Q23E`.

The complete numerical record is in `bmad_scibmad_comparison.csv`. The largest
remaining local matrix discrepancies are now in the sector bends, followed by
the two overlap elements and the wigglers.
