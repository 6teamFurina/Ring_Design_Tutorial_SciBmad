#!/usr/bin/env python3
"""Export a CESR Bmad/Tao reference data set for SciBmad comparison.

Run this script in the Linux environment that has Bmad, Tao, and PyTao.
The output JSON contains the local map of every tracking element and the
cumulative map from BEGINNING to every element exit.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional


def as_jsonable(value: Any) -> Any:
    """Convert numpy/PyTao values into objects accepted by json.dump."""
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        return value if math.isfinite(value) else str(value)
    if hasattr(value, "item"):
        try:
            return as_jsonable(value.item())
        except (TypeError, ValueError):
            pass
    if hasattr(value, "tolist"):
        return as_jsonable(value.tolist())
    if isinstance(value, dict):
        return {str(key): as_jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [as_jsonable(item) for item in value]
    if hasattr(value, "model_dump"):
        return as_jsonable(value.model_dump())
    if hasattr(value, "dict"):
        return as_jsonable(value.dict())
    return str(value)


def safe_lat_list(tao: Any, who: str, count: int) -> List[Any]:
    try:
        values = tao.lat_list(
            "*",
            who,
            which="model",
            flags="-array_out -track_only",
        )
        result = list(as_jsonable(values))
        if len(result) == count:
            return result
        print(
            "WARNING: query {!r} returned {} values instead of {}".format(
                who, len(result), count
            ),
            file=sys.stderr,
        )
    except Exception as exc:
        print("WARNING: query {!r} failed: {}".format(who, exc), file=sys.stderr)
    return [None] * count


def normalize_map_result(result: Any) -> Dict[str, Any]:
    """Normalize PyTao's dict and model representations of a linear map."""
    if result is None:
        raise RuntimeError("Tao returned no map")
    if isinstance(result, dict):
        matrix = result.get("mat6", result.get("matrix"))
        vector = result.get("vec0", result.get("vector"))
        error = result.get("symplectic_error")
    else:
        matrix = getattr(result, "mat6", None)
        vector = getattr(result, "vec0", None)
        error = getattr(result, "symplectic_error", None)
    if matrix is None:
        raise RuntimeError("Tao map result has no mat6 field: {!r}".format(result))
    return {
        "vec0": as_jsonable(vector),
        "mat6": as_jsonable(matrix),
        "symplectic_error": as_jsonable(error),
    }


def local_element_map(tao: Any, index: int) -> Dict[str, Any]:
    """Read the map across one element, supporting old and new PyTao APIs."""
    if hasattr(tao, "ele"):
        element = tao.ele(
            index,
            which="model",
            defaults=False,
            mat6=True,
            warn=False,
        )
        return normalize_map_result(element.mat6)
    return normalize_map_result(tao.ele_mat6(str(index), which="model"))


def cumulative_map(tao: Any, index: int) -> Dict[str, Any]:
    """Read the map from the lattice beginning to an element's downstream end."""
    if index == 0:
        return {
            "vec0": [0.0] * 6,
            "mat6": [
                [1.0 if row == column else 0.0 for column in range(6)]
                for row in range(6)
            ],
            "symplectic_error": 0.0,
        }
    return normalize_map_result(tao.matrix("beginning", str(index)))


def make_tao(lattice: Path) -> Any:
    from pytao import Tao

    try:
        return Tao(lattice_file=str(lattice), noplot=True)
    except TypeError:
        # Compatibility with older PyTao releases.
        return Tao("-lat {} -noplot".format(lattice))


def export_reference(lattice: Path, output_dir: Path) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    tao = make_tao(lattice)

    try:
        version = as_jsonable(tao.version())
    except Exception:
        version = as_jsonable(tao.cmd("show version", raises=False))

    indices = list(
        as_jsonable(
            tao.lat_list(
                "*",
                "ele.ix_ele",
                which="model",
                flags="-array_out -track_only",
            )
        )
    )
    indices = [int(index) for index in indices]
    count = len(indices)
    if count < 2:
        raise RuntimeError("Tao returned too few tracking elements: {}".format(count))

    columns = {
        "name": safe_lat_list(tao, "ele.name", count),
        "key": safe_lat_list(tao, "ele.key", count),
        "s": safe_lat_list(tao, "ele.s", count),
        "length": safe_lat_list(tao, "ele.l", count),
        "p0c": safe_lat_list(tao, "ele.p0c", count),
        "x": safe_lat_list(tao, "orbit.vec.1", count),
        "px": safe_lat_list(tao, "orbit.vec.2", count),
        "y": safe_lat_list(tao, "orbit.vec.3", count),
        "py": safe_lat_list(tao, "orbit.vec.4", count),
        "z": safe_lat_list(tao, "orbit.vec.5", count),
        "pz": safe_lat_list(tao, "orbit.vec.6", count),
        "beta_a": safe_lat_list(tao, "twiss.beta_a", count),
        "alpha_a": safe_lat_list(tao, "twiss.alpha_a", count),
        "phi_a": safe_lat_list(tao, "twiss.phi_a", count),
        "beta_b": safe_lat_list(tao, "twiss.beta_b", count),
        "alpha_b": safe_lat_list(tao, "twiss.alpha_b", count),
        "phi_b": safe_lat_list(tao, "twiss.phi_b", count),
        "eta_x": safe_lat_list(tao, "twiss.eta_x", count),
        "etap_x": safe_lat_list(tao, "twiss.etap_x", count),
        "eta_y": safe_lat_list(tao, "twiss.eta_y", count),
        "etap_y": safe_lat_list(tao, "twiss.etap_y", count),
    }

    elements: List[Dict[str, Any]] = []
    errors: List[Dict[str, Any]] = []
    for position, index in enumerate(indices):
        record = {"index": index}
        for name, values in columns.items():
            record[name] = values[position]

        try:
            record["local_map"] = local_element_map(tao, index)
        except Exception as exc:
            message = "local map failed: {}".format(exc)
            record["local_map_error"] = message
            errors.append({"index": index, "stage": "local_map", "error": str(exc)})

        try:
            record["cumulative_map"] = cumulative_map(tao, index)
        except Exception as exc:
            message = "cumulative map failed: {}".format(exc)
            record["cumulative_map_error"] = message
            errors.append(
                {"index": index, "stage": "cumulative_map", "error": str(exc)}
            )

        elements.append(record)
        if position % 25 == 0 or position + 1 == count:
            print("Exported {}/{} tracking positions".format(position + 1, count))

    try:
        one_turn = normalize_map_result(tao.matrix("beginning", "end"))
    except Exception as exc:
        one_turn = {"error": str(exc)}
        errors.append({"index": None, "stage": "one_turn", "error": str(exc)})

    try:
        ring_general = as_jsonable(tao.ring_general())
    except Exception as exc:
        ring_general = {"error": str(exc)}

    payload = {
        "format": "cesr-bmad-reference-v1",
        "lattice": str(lattice),
        "tao_version": version,
        "tracking_position_count": count,
        "coordinate_order": ["x", "px", "y", "py", "z", "pz"],
        "map_location": "downstream end of each tracking element",
        "ring_general": ring_general,
        "one_turn_map": one_turn,
        "elements": elements,
        "errors": errors,
    }

    json_path = output_dir / "bmad_reference.json"
    with json_path.open("w", encoding="utf-8") as stream:
        json.dump(payload, stream, indent=2, sort_keys=False, allow_nan=False)
        stream.write("\n")

    csv_path = output_dir / "element_index.csv"
    csv_fields = [
        "index",
        "name",
        "key",
        "s",
        "length",
        "p0c",
        "x",
        "px",
        "y",
        "py",
        "z",
        "pz",
        "beta_a",
        "alpha_a",
        "phi_a",
        "beta_b",
        "alpha_b",
        "phi_b",
        "eta_x",
        "etap_x",
        "eta_y",
        "etap_y",
    ]
    with csv_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=csv_fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(elements)

    print("Wrote {}".format(json_path))
    print("Wrote {}".format(csv_path))
    if errors:
        print("WARNING: {} map queries failed; see the errors array".format(len(errors)))
    return json_path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lattice", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    lattice = args.lattice.expanduser().resolve()
    output = args.output.expanduser().resolve()
    if not lattice.is_file():
        parser.error("lattice file does not exist: {}".format(lattice))

    export_reference(lattice, output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
