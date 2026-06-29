# Ring_Design_Tutorial_SciBmad

A ring design tutorial in SciBmad/Julia.

This tutorial introduces SciBmad for the simulation and optimization of particle
accelerators. The numbered Jupyter notebooks should be read in order. Together
they build an example storage ring similar to the Electron Storage Ring of the
Electron-Ion Collider. The main tutorial text is followed by exercises, with
example solutions provided for comparison. Readers are encouraged to try the
exercises before opening the solutions.

## Reading the Tutorial

To read a chapter, open the corresponding numbered `.ipynb` notebook. The
notebooks already include the tutorial text, source code, and saved output from
the code cells, so they can be read without running anything.

The `assets/` directory contains images loaded by the tutorial notebooks. Keep
it next to the notebooks when viewing them locally so the figures render
correctly.

## Running a Notebook Yourself

To re-run a chapter, download the notebook and open it in JupyterLab from the
root of this project directory. Running from the project root keeps the relative
paths used in the notebooks valid.

Downloading commands for necessary packages are shown when they are used in the tutorial first time.

For chapter `N`, also download:

- `lattices/common/`
- the non-exercise lattice files in `lattices/chapter_1/` through
  `lattices/chapter_N/`
- `assets/`, if the notebook displays tutorial figures or writes updated plots

Exercise files are kept with their chapter under `lattices/chapter_*`. Download
the exercise solution codes/notebooks and exercise output results if you want to work
through those exercises.

## Repository Layout

- numbered `chapter*.ipynb` files: main tutorial chapters
- `lattices/common/`: shared lattices and helper scripts used by later chapters
- `lattices/chapter_*/`: chapter-specific lattices, generated solutions, and
  exercise files
- `assets/`: figures loaded by the tutorial notebooks
