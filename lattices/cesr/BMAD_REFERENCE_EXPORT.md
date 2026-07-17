# Exporting the CESR Bmad reference on Linux

The Bmad/Tao side of the comparison must run in the Linux environment where
Tao is installed. SciBmad can remain on the Windows machine. The exporter
records maps in the Bmad coordinate order `(x, px, y, py, z, pz)`.

## Files to copy to the Linux VM

- `run_bmad_reference.sh`
- `export_bmad_reference.py`

The existing Linux directory must also contain `cesr.bmad`.

From a Windows PowerShell prompt, not from inside the SSH session, one possible
copy command is:

    scp run_bmad_reference.sh export_bmad_reference.py \
      jn577@lnx201.classe.cornell.edu:/home/jn577/cesr_scibmad/

## Run on Linux

    cd /home/jn577/cesr_scibmad
    chmod +x run_bmad_reference.sh
    ./run_bmad_reference.sh /home/jn577/cesr_scibmad/cesr.bmad

The preferred path uses PyTao and may take several minutes because it requests
both a local and a cumulative 6-by-6 map at every tracking position. It writes:

- `bmad_reference_output/bmad_reference.json`
- `bmad_reference_output/element_index.csv`
- `bmad_reference_output/run.log`
- `bmad_reference_output.tar.gz`

If PyTao is not installed, the wrapper automatically asks the Tao command-line
interface for all element transfer matrices and writes `tao_fallback.txt`.

## Copy the result back to Windows

Run this from Windows PowerShell:

    scp jn577@lnx201.classe.cornell.edu:/home/jn577/cesr_scibmad/bmad_reference_output.tar.gz .

Return the archive without modifying its contents. The JSON contains enough
information to compare Bmad and SciBmad element maps and to identify the first
location where their cumulative maps diverge.

## Run the local SciBmad comparison

From the repository root on the machine with Julia and SciBmad installed:

    julia --project=. lattices/cesr/test_bmad_scibmad.jl

The test validates the archive, aligns Bmad tracking indices 1 through 869 with
the SciBmad line, and uses GTPSA to linearize every SciBmad element about the
same Bmad entrance trajectory used for the corresponding Tao matrix. It writes
the per-element diagnostics to `bmad_scibmad_comparison.csv`.

For a quick archive check or a short smoke test, use:

    julia --project=. lattices/cesr/test_bmad_scibmad.jl --parse-only
    julia --project=. lattices/cesr/test_bmad_scibmad.jl --max-elements=10

The normal run is diagnostic and succeeds after producing the report even when
physics discrepancies exist. Add `--strict` when the models are expected to
agree and the script should act as a numerical regression gate. The strict
matrix tolerance is `1e-6`, reflecting the seven-decimal matrix precision in
the Tao CLI fallback output.
