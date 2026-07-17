#!/usr/bin/env bash
# Run in the Linux VM that contains Bmad/Tao. The preferred structured export
# uses PyTao. If PyTao is unavailable, a raw Tao CLI fallback is produced.

set -u
set -o pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
lattice="${1:-${script_dir}/cesr.bmad}"
output_dir="${2:-${script_dir}/bmad_reference_output}"
archive="${output_dir}.tar.gz"

if [[ ! -f "${lattice}" ]]; then
  echo "ERROR: lattice file not found: ${lattice}" >&2
  echo "Usage: $0 [/absolute/path/to/cesr.bmad] [output_directory]" >&2
  exit 1
fi

mkdir -p "${output_dir}"
exec > >(tee "${output_dir}/run.log") 2>&1

echo "CESR Bmad reference export"
echo "lattice=${lattice}"
echo "output_dir=${output_dir}"
echo "date=$(date --iso-8601=seconds)"
echo "host=$(hostname)"
echo "user=$(id -un)"
echo "tao=$(command -v tao || true)"
echo "python3=$(command -v python3 || true)"
tao -version 2>&1 || true
python3 --version 2>&1 || true

if ! command -v tao >/dev/null 2>&1; then
  echo "ERROR: tao is not on PATH." >&2
  exit 1
fi

structured_ok=false
if command -v python3 >/dev/null 2>&1 && python3 -c 'import pytao' >/dev/null 2>&1; then
  echo "PyTao found. Starting structured per-element export."
  if python3 "${script_dir}/export_bmad_reference.py" \
      --lattice "${lattice}" \
      --output "${output_dir}"; then
    structured_ok=true
  else
    echo "WARNING: structured PyTao export failed. Running Tao CLI fallback."
  fi
else
  echo "WARNING: PyTao is not installed. Running Tao CLI fallback."
fi

if [[ "${structured_ok}" != true ]]; then
  command_file="${output_dir}/tao_fallback.commands"
  raw_output="${output_dir}/tao_fallback.txt"
  {
    echo "show version"
    echo "show global"
    echo "show lattice"
    echo "show ele WIG_W -xfer_mat"
    echo "show ele WIG_E -xfer_mat"
    echo "show ele Q00W -xfer_mat"
    echo "show ele Q00E -xfer_mat"
    echo "show ele CLEO_SOL -xfer_mat"
    echo "show ele B03W -xfer_mat"
    echo "show ele B03E -xfer_mat"
    for index in $(seq 0 869); do
      echo "show ele ${index} -xfer_mat"
    done
    echo "quit"
  } > "${command_file}"

  tao -noinit -noplot -lat "$(realpath "${lattice}")" \
    < "${command_file}" > "${raw_output}" 2>&1 || true
  echo "Raw Tao fallback written to ${raw_output}"
fi

tar -czf "${archive}" -C "$(dirname "${output_dir}")" "$(basename "${output_dir}")"
echo "Archive ready: ${archive}"
echo "Copy this .tar.gz file back to the Windows CESR folder."
