"""
utils/inspect_file.py
---------------------
Tools for inspecting unknown MEA data files.

Designed for new users who want to:
  1. See the full internal structure of any HDF5 or BXR file
     before writing a loader for a new manufacturer.
  2. Compare two files side by side (e.g. two manufacturers)
     to find what keys are shared or different.

Usage (command line)
--------------------
# Dump the structure of one or more files:
    python -m utils.inspect_file dump recording.bxr
    python -m utils.inspect_file dump data.raw.h5
    python -m utils.inspect_file dump file1.h5 file2.h5

# Compare two files:
    python -m utils.inspect_file compare 3brain_file.bxr maxwell_file.h5

Output is saved as <filename>.structure.txt next to the data file,
and also printed to the terminal.
"""

import os
import sys
import h5py


# ── single file dump ──────────────────────────────────────────────────────────
def dump_structure(path: str, out_txt: str = None, print_to_screen: bool = True):
    """
    Traverse every group and dataset in an HDF5 / BXR file and write
    a human-readable structure report.

    Parameters
    ----------
    path          : path to the .h5 or .bxr file
    out_txt       : where to save the report (default: <filename>.structure.txt
                    in the same folder as the data file)
    print_to_screen: also print to terminal (default True)

    Returns
    -------
    Path to the saved .txt file.
    """
    if out_txt is None:
        base    = os.path.splitext(os.path.basename(path))[0]
        folder  = os.path.dirname(os.path.abspath(path))
        out_txt = os.path.join(folder, base + ".structure.txt")

    lines = []

    def add(text=""):
        lines.append(text)
        if print_to_screen:
            print(text)

    add(f"FILE : {path}")
    add("=" * 70)

    with h5py.File(path, "r") as f:

        # root-level attributes (e.g. SamplingRate in 3Brain)
        if f.attrs:
            add("\nROOT ATTRIBUTES")
            add("-" * 40)
            for k, v in f.attrs.items():
                add(f"  {k} = {v}")

        add("\nSTRUCTURE")
        add("-" * 40)

        def visit(name, obj):
            depth  = name.count("/")
            indent = "  " * depth
            short  = name.split("/")[-1]

            if isinstance(obj, h5py.Group):
                add(f"{indent}[GROUP]  {short}/")
                if obj.attrs:
                    for k, v in obj.attrs.items():
                        add(f"{indent}  attr: {k} = {v}")

            elif isinstance(obj, h5py.Dataset):
                # structured dtypes (e.g. spike tables with named fields)
                if obj.dtype.names:
                    fields = ", ".join(
                        f"{n}:{obj.dtype[n]}" for n in obj.dtype.names
                    )
                    add(f"{indent}[DATA]   {short}  "
                        f"shape={obj.shape}  fields=({fields})")
                else:
                    add(f"{indent}[DATA]   {short}  "
                        f"shape={obj.shape}  dtype={obj.dtype}")
                if obj.attrs:
                    for k, v in obj.attrs.items():
                        add(f"{indent}  attr: {k} = {v}")

        f.visititems(visit)

    add("")
    add(f"Total lines: {len(lines)}")

    with open(out_txt, "w", encoding="utf-8") as w:
        w.write("\n".join(lines) + "\n")

    if print_to_screen:
        print(f"\n[OK] Saved to: {out_txt}")

    return out_txt


# ── two-file comparison ───────────────────────────────────────────────────────
def compare_structures(path_a: str, path_b: str, print_to_screen: bool = True):
    """
    Compare the dataset keys of two HDF5/BXR files.
    Useful for mapping a new manufacturer's format to an existing loader.

    Shows:
      - Datasets only in file A
      - Datasets only in file B
      - Datasets present in both (with their dtypes)

    Parameters
    ----------
    path_a / path_b : paths to the two files
    """

    def collect(path):
        datasets = {}
        with h5py.File(path, "r") as f:
            def visit(name, obj):
                if isinstance(obj, h5py.Dataset):
                    if obj.dtype.names:
                        dtype_str = "structured(" + ",".join(obj.dtype.names) + ")"
                    else:
                        dtype_str = str(obj.dtype)
                    datasets[name] = dtype_str
            f.visititems(visit)
        return datasets

    d1 = collect(path_a)
    d2 = collect(path_b)

    only_a  = sorted(set(d1) - set(d2))
    only_b  = sorted(set(d2) - set(d1))
    in_both = sorted(set(d1) & set(d2))

    lines = []

    def add(text=""):
        lines.append(text)
        if print_to_screen:
            print(text)

    add("=" * 70)
    add(f"FILE A: {path_a}")
    add(f"FILE B: {path_b}")
    add("=" * 70)

    add(f"\nDATASETS only in A  ({len(only_a)} total)")
    add("-" * 40)
    for k in only_a:
        add(f"  {k}  [{d1[k]}]")

    add(f"\nDATASETS only in B  ({len(only_b)} total)")
    add("-" * 40)
    for k in only_b:
        add(f"  {k}  [{d2[k]}]")

    add(f"\nDATASETS in both  ({len(in_both)} total)")
    add("-" * 40)
    for k in in_both:
        same = "same dtype" if d1[k] == d2[k] else f"A={d1[k]}  B={d2[k]}"
        add(f"  {k}  [{same}]")

    # save report
    out_txt = "compare_" + \
        os.path.splitext(os.path.basename(path_a))[0] + "_vs_" + \
        os.path.splitext(os.path.basename(path_b))[0] + ".txt"

    with open(out_txt, "w", encoding="utf-8") as w:
        w.write("\n".join(lines) + "\n")

    if print_to_screen:
        print(f"\n[OK] Saved to: {out_txt}")

    return out_txt


# ── command-line interface ────────────────────────────────────────────────────
if __name__ == "__main__":
    usage = (
        "Usage:\n"
        "  python -m utils.inspect_file dump <file1> [file2] ...\n"
        "  python -m utils.inspect_file compare <fileA> <fileB>\n"
    )

    if len(sys.argv) < 3:
        print(usage)
        sys.exit(1)

    command = sys.argv[1].lower()

    if command == "dump":
        for path in sys.argv[2:]:
            if not os.path.exists(path):
                print(f"[FAIL] File not found: {path}")
                continue
            dump_structure(path)

    elif command == "compare":
        if len(sys.argv) != 4:
            print("compare needs exactly two files.\n" + usage)
            sys.exit(1)
        compare_structures(sys.argv[2], sys.argv[3])

    else:
        print(f"Unknown command '{command}'\n" + usage)
        sys.exit(1)
