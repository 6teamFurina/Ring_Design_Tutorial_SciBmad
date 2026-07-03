# SciBmad Ring Design Tutorial

A ring design tutorial in SciBmad/Julia.

This tutorial introduces SciBmad for the simulation and optimization of particle
accelerators. The numbered Jupyter notebooks should be read in order. Together
they build an example storage ring similar to the Electron Storage Ring of the
Electron-Ion Collider. The main tutorial text is followed by exercises, with
example solutions provided for comparison. Readers are encouraged to try the
exercises before opening the solutions.

## Chapters

| # | Notebook | Topic |
|---|---|---|
| 0 | `chapter00_power_series_and_optimization_scibmad.ipynb` | Power series, differentiation, and optimizers |
| 1 | `chapter01_fodo_scibmad.ipynb` | FODO cells |
| 2 | `chapter02_dispersion_suppressor_scibmad.ipynb` | Dispersion suppressor |
| 3 | `chapter03_twiss_matching_scibmad.ipynb` | Matching the dispersion suppressor to the straight section |
| 4 | `chapter04_machine_coordinates_scibmad.ipynb` | Machine coordinates in SciBmad |
| 5 | `chapter05_constructing_the_ring_scibmad.ipynb` | Constructing the ring |
| 6 | `chapter06_low_beta_ir_scibmad.ipynb` | Low-beta interaction region insertion |
| 7 | `chapter07_tune_cell_scibmad.ipynb` | Tune cell |
| 8 | `chapter08_phase_space_scibmad.ipynb` | Particle phase-space coordinates |
| 9 | `chapter09_rf_cavities_scibmad.ipynb` | RF cavities |
| 10 | `chapter10_long_term_tracking_scibmad.ipynb` | Long-term tracking |
| 11 | `chapter11_control_elements_scibmad.ipynb` | Control elements |
| 12 | `chapter12_dynamic_aperture_scibmad.ipynb` | Dynamic aperture |
| 13 | `chapter13_nonlinear_twiss_scibmad.ipynb` | Nonlinear Twiss |
| 14 | `chapter14_model_design_base_lattices_scibmad.ipynb` | Model, design, and base lattices |
| 15 | `chapter15_orbit_correction_scibmad.ipynb` | Orbit correction |
| 16 | `chapter16_error_fitting_scibmad.ipynb` | Error fitting |
| 17 | `chapter17_spin_tracking_with_ramping_scibmad.ipynb` | Spin tracking with ramping |

## Which Lattice Each Chapter Uses

The chapters fall into three groups. Knowing which group a chapter belongs to
explains why some chapters build on each other while others stand alone.

- **Building the example ring (chapters 1–3, 5–7, 9).** These chapters
  progressively construct one EIC-ESR-like storage ring. Each stage optimizes a
  piece and saves its strengths to a small solution file under `lattices/`, which
  the next stage loads. This is the main through-line of the tutorial.
- **Standalone teaching examples (chapters 0, 4, 8, 11, 14).** These use a small
  purpose-built example (a toy objective, a single misaligned element, one
  bend-and-quadrupole line, etc.) to isolate one concept. They do not use the
  example ring.
- **Separate pre-built lattices (chapters 10, 12, 13, 15, 16, 17).** These
  demonstrate a capability on an independent lattice supplied with the chapter,
  rather than on the ring built in chapters 1–9.

| Chapter | Lattice / model | Group |
|---|---|---|
| 0 | toy two-variable objective | standalone example |
| 1 | forward arc FODO cell, built from scratch → `chapter1_fodoF_solution.jl` | builds the ring |
| 2 | dispersion suppressor; loads chapter 1 → `chapter2_dispsupF_solution.jl` | builds the ring |
| 3 | matching section; loads chapter 2 → `chapter3_mSSF_solution.jl` | builds the ring |
| 4 | single misaligned quadrupole / patch element | standalone example |
| 5 | full ring assembled from the chapter 1–3 solutions → `chapter5_ring_definition.jl` | builds the ring |
| 6 | low-beta interaction region inserted into the chapter 5 ring | builds the ring |
| 7 | tune cell; loads `chapter5_ring_definition.jl` + `chapter6_IR_solution.jl` | builds the ring |
| 8 | one bend-and-quadrupole line, single-particle tracking | standalone example |
| 9 | RF cavities inserted into `chapter5_ring_definition.jl` | builds the ring |
| 10 | compact 16-cell electron storage ring, built inline | separate lattice |
| 11 | quadrupole-and-bend control example | standalone example |
| 12 | compact demonstration ring + ESR-style lattice (`esr-da-opt.jl`) | separate lattice |
| 13 | ESR v6.3.1 lattice (`esr-v6.3.1.jl`) | separate lattice |
| 14 | drift–bend–quadrupole example | standalone example |
| 15 | sawtooth `ring0` beamline (`chapter15_b_sawtooth_ring0_beamline.jl`) | separate lattice |
| 16 | RCS lattice (`RCSV5S0.jl`) | separate lattice |
| 17 | AGS-like spin lattice (`spin_lat.bmad` / `ags.jl`) | separate lattice |

## Reading the Tutorial

To read a chapter, open the corresponding numbered `.ipynb` notebook. The
notebooks already include the tutorial text, source code, and saved output from
the code cells, so they can be read without running anything.

The `assets/` directory contains images loaded by the tutorial notebooks. Keep
it next to the notebooks when viewing them locally so the figures render
correctly.

## Running a Notebook Yourself

The notebooks share a single Julia environment defined by the `Project.toml` at
the project root.

Some SciBmad dependencies (for example `SimUtils`) live in the `bmad-sim`
package registry rather than the Julia General registry, so add that registry
once before instantiating:

```julia
import Pkg
Pkg.Registry.add("General")  # normally already present
Pkg.Registry.add(Pkg.RegistrySpec(url = "https://github.com/bmad-sim/BmadRegistry.jl"))
```

Then, rather than installing packages one at a time, activate and instantiate
the shared environment once from the project root:

```julia
import Pkg
Pkg.activate(".")
Pkg.instantiate()
```

The `[compat]` bounds record the package versions the notebooks were last run
against. (Individual notebooks still contain `Pkg.add` cells from an earlier
per-chapter setup; with the shared environment activated these are redundant and
can be skipped.)

Then start JupyterLab from the project root so the relative paths used in the
notebooks stay valid.

For chapter `N`, also make sure you have:

- `lattices/common/`
- the non-exercise lattice files in `lattices/chapter_1/` through
  `lattices/chapter_N/` (not every chapter has a lattice directory)
- `assets/`, if the notebook displays tutorial figures or writes updated plots

Exercise files are kept with their chapter under `lattices/chapter_*`. Download
the exercise solution codes/notebooks and exercise output results if you want to
work through those exercises.

## Repository Layout

- numbered `chapter*.ipynb` files: main tutorial chapters
- `Project.toml`: shared Julia environment for all chapters
- `lattices/common/`: shared lattices and helper scripts used by later chapters
- `lattices/chapter_*/`: chapter-specific lattices and results generated by tutorial/solution codes
- `lattices/chapter_*/exercises`: solution codes to exercises
- `assets/`: figures loaded by the tutorial notebooks

## License

This tutorial is licensed under the
[Creative Commons Attribution-ShareAlike 4.0 International License](https://creativecommons.org/licenses/by-sa/4.0/)
(CC BY-SA 4.0). See the [`LICENSE`](LICENSE) file for the full text.
