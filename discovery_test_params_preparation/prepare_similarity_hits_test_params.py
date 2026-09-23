#!/usr/bin/env python3

"""
prepare_similarity_hits_test_params.py

Prepare Rosetta test_params directories directly from the output CSV
of the Enamine similarity-search pipeline.

Expected input CSV format:

    SMILES,ligand_name,source_bz2,tanimoto

The CSV may either be headerless or contain:

    smiles,ligand_name,source_bz2,tanimoto

Output organization:

    OUTPUT_ROOT/
        0/
            ligand_1/
                test_params/
            ligand_2/
                test_params/
            ...
        1/
            ...
        2/
            ...

By default, each numbered organizational directory contains up to
100 ligands.

This script is designed to run as an LSF array worker.

Example:

    python prepare_similarity_hits_test_params.py \
        --input-csv top_1m_pipeline.csv \
        --output-root /path/to/refined_conformers \
        --task-index 1 \
        --ligands-per-task 1

If --task-index is not supplied, LSB_JOBINDEX is used.
"""

import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
import time

from itertools import islice
from pathlib import Path


# ============================================================
# Utility functions
# ============================================================

def run_command(cmd, cwd=None):
    """
    Run a command and fail immediately if it returns non-zero.
    """

    print(
        "[CMD] " + " ".join(str(x) for x in cmd),
        flush=True,
    )

    subprocess.run(
        [str(x) for x in cmd],
        cwd=cwd,
        check=True,
    )


def looks_like_header(row):
    """
    Determine whether the first CSV row is a header.
    """

    if not row:
        return False

    first = row[0].strip().lower()

    return first in {
        "smiles",
        "canonical_smiles",
        "canonicalsmiles",
    }


def count_records(csv_path):
    """
    Count ligand records in the CSV, excluding a header if present.

    Intended mainly for standalone validation/debugging.
    The LSF worker itself does not need this.
    """

    count = 0

    with open(csv_path, "r", newline="") as handle:

        reader = csv.reader(handle)

        first = next(reader, None)

        if first is None:
            return 0

        if not looks_like_header(first):
            count += 1

        for row in reader:
            if row:
                count += 1

    return count


def iter_assigned_records(
    csv_path,
    task_index,
    ligands_per_task,
):
    """
    Yield only records assigned to this LSF array element.

    rank is zero-based relative to ligand records, NOT including
    the optional header.

    Assignment:

        task 1 -> records 0 ... N-1
        task 2 -> records N ... 2N-1
        ...

    where N = ligands_per_task.
    """

    start_record = (
        (task_index - 1)
        * ligands_per_task
    )

    end_record = (
        start_record
        + ligands_per_task
    )

    with open(
        csv_path,
        "r",
        newline="",
    ) as handle:

        reader = csv.reader(handle)

        first = next(reader, None)

        if first is None:
            return

        # ----------------------------------------------------
        # Convert reader into a clean stream of ligand rows.
        # ----------------------------------------------------

        def rows():

            if not looks_like_header(first):
                yield first

            for row in reader:

                if row:
                    yield row

        # ----------------------------------------------------
        # Skip directly to this task's assigned section.
        # ----------------------------------------------------

        selected_rows = islice(
            rows(),
            start_record,
            end_record,
        )

        for offset, row in enumerate(
            selected_rows
        ):

            rank = (
                start_record
                + offset
            )

            if len(row) < 4:

                raise RuntimeError(
                    f"CSV record {rank + 1} has "
                    f"{len(row)} fields; expected at least 4.\n"
                    f"Record: {row}"
                )

            smiles = row[0].strip()
            ligand_name = row[1].strip()
            source_bz2 = row[2].strip()

            try:
                tanimoto = float(
                    row[3].strip()
                )

            except ValueError:

                raise RuntimeError(
                    f"Invalid Tanimoto value at "
                    f"record {rank + 1}: "
                    f"{row[3]}"
                )

            if not smiles:
                raise RuntimeError(
                    f"Missing SMILES at record "
                    f"{rank + 1}"
                )

            if not ligand_name:
                raise RuntimeError(
                    f"Missing ligand name at record "
                    f"{rank + 1}"
                )

            yield {
                "rank": rank,
                "smiles": smiles,
                "ligand_name": ligand_name,
                "source_bz2": source_bz2,
                "tanimoto": tanimoto,
            }


# ============================================================
# Validate finished test_params
# ============================================================

def test_params_complete(path):
    """
    Basic validation of an existing test_params directory.

    Require:
        residue_types.txt
        exclude_pdb_component_list.txt
        patches.txt
        at least one .params
    """

    path = Path(path)

    if not path.is_dir():
        return False

    required = [
        "residue_types.txt",
        "exclude_pdb_component_list.txt",
        "patches.txt",
    ]

    for filename in required:

        if not (
            path / filename
        ).is_file():
            return False

    params_files = list(
        path.glob("*.params")
    )

    if not params_files:
        return False

    return True


# ============================================================
# Process one ligand
# ============================================================

def prepare_ligand(
    record,
    output_root,
    ligands_per_directory,
    conformator_container,
    conformator_executable,
    molfile_to_params,
    obabel_executable,
    max_conformers,
    conformator_license,
    overwrite,
):
    """
    Generate Rosetta test_params for one ligand.
    """

    rank = record["rank"]

    smiles = record["smiles"]
    ligand_name = record["ligand_name"]
    source_bz2 = record["source_bz2"]
    tanimoto = record["tanimoto"]

    # --------------------------------------------------------
    # Organizational directory.
    #
    # rank 0-99   -> 0
    # rank 100-199 -> 1
    # etc.
    # --------------------------------------------------------

    group_index = (
        rank
        // ligands_per_directory
    )

    group_dir = (
        output_root
        / str(group_index)
    )

    ligand_dir = (
        group_dir
        / ligand_name
    )

    final_test_params = (
        ligand_dir
        / "test_params"
    )

    staging_test_params = (
        ligand_dir
        / "test_params.inprogress"
    )

    group_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    ligand_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    print(
        "============================================================",
        flush=True,
    )
    print(
        f"Rank:         {rank + 1}",
        flush=True,
    )
    print(
        f"Ligand:       {ligand_name}",
        flush=True,
    )
    print(
        f"Tanimoto:     {tanimoto:.10f}",
        flush=True,
    )
    print(
        f"Group:        {group_index}",
        flush=True,
    )
    print(
        f"Output:       {ligand_dir}",
        flush=True,
    )
    print(
        "============================================================",
        flush=True,
    )

    # --------------------------------------------------------
    # Resume behavior.
    # --------------------------------------------------------

    if (
        final_test_params.exists()
        and test_params_complete(
            final_test_params
        )
        and not overwrite
    ):

        print(
            f"[SKIP] Complete test_params already exists "
            f"for {ligand_name}",
            flush=True,
        )

        return {
            "ligand_name": ligand_name,
            "status": "skipped",
            "params_count": len(
                list(
                    final_test_params.glob(
                        "*.params"
                    )
                )
            ),
        }

    if final_test_params.exists():

        if overwrite:

            print(
                f"[INFO] Removing existing "
                f"{final_test_params}",
                flush=True,
            )

            shutil.rmtree(
                final_test_params
            )

        else:

            raise RuntimeError(
                f"Existing test_params directory for "
                f"{ligand_name} appears incomplete:\n"
                f"  {final_test_params}\n"
                f"Use --overwrite to rebuild it."
            )

    # --------------------------------------------------------
    # Remove abandoned staging directory from prior failure.
    # --------------------------------------------------------

    if staging_test_params.exists():

        print(
            f"[INFO] Removing previous incomplete staging "
            f"directory for {ligand_name}",
            flush=True,
        )

        shutil.rmtree(
            staging_test_params
        )

    staging_test_params.mkdir(
        parents=True
    )

    start_time = time.time()

    # --------------------------------------------------------
    # Save provenance.
    # --------------------------------------------------------

    metadata = {
        "rank": rank + 1,
        "ligand_name": ligand_name,
        "smiles": smiles,
        "source_bz2": source_bz2,
        "tanimoto": tanimoto,
    }

    with open(
        ligand_dir
        / "similarity_hit.json",
        "w",
    ) as handle:

        json.dump(
            metadata,
            handle,
            indent=2,
        )

    # --------------------------------------------------------
    # Write input SMILES.
    #
    # Include the ligand name as a second field. This is
    # generally safe/useful for molecule tools and preserves ID.
    # --------------------------------------------------------

    smiles_file = (
        staging_test_params
        / f"{ligand_name}.smi"
    )

    with open(
        smiles_file,
        "w",
    ) as handle:

        handle.write(
            f"{smiles}\t{ligand_name}\n"
        )

    # --------------------------------------------------------
    # Conformator output.
    # --------------------------------------------------------

    conformer_sdf = (
        staging_test_params
        / f"{ligand_name}_confs.sdf"
    )

    conformator_cmd_inside = [
        conformator_executable,
        "-i",
        smiles_file.name,
        "-o",
        conformer_sdf.name,
        "--keep3d",
        "--hydrogens",
        "-v",
        "0",
    ]

    # If your Conformator supports an explicit conformer-count
    # argument, add it here.
    #
    # The older pipeline relied on Conformator's default
    # (~250) rather than explicitly passing a count.
    #
    # max_conformers is currently metadata/configuration only
    # unless the appropriate command-line flag is supplied.
    #
    # See note after the script.

    if conformator_license:

        activation_cmd = [
            "singularity",
            "exec",
            str(conformator_container),
            conformator_executable,
            "--license",
            conformator_license,
        ]

        run_command(
            activation_cmd,
            cwd=staging_test_params,
        )

    conformator_cmd = [
        "singularity",
        "exec",
        str(conformator_container),
    ] + conformator_cmd_inside

    run_command(
        conformator_cmd,
        cwd=staging_test_params,
    )

    if not conformer_sdf.is_file():

        raise RuntimeError(
            f"Conformator did not generate expected file:\n"
            f"  {conformer_sdf}"
        )

    # --------------------------------------------------------
    # Split multi-conformer SDF.
    # --------------------------------------------------------

    split_prefix = (
        f"{ligand_name}_.sdf"
    )

    run_command(
        [
            obabel_executable,
            "-isdf",
            conformer_sdf.name,
            "-O",
            split_prefix,
            "-m",
        ],
        cwd=staging_test_params,
    )

    # --------------------------------------------------------
    # Find individual conformer SDFs.
    # --------------------------------------------------------

    individual_sdfs = sorted(
        path
        for path
        in staging_test_params.glob(
            f"{ligand_name}_*.sdf"
        )
        if path.name
        != conformer_sdf.name
    )

    if not individual_sdfs:

        raise RuntimeError(
            f"No individual conformer SDF files "
            f"were produced for {ligand_name}."
        )

    print(
        f"[INFO] {len(individual_sdfs):,} "
        f"conformers generated.",
        flush=True,
    )

    # --------------------------------------------------------
    # Convert every conformer to Rosetta params.
    # --------------------------------------------------------

    params_names = []

    for sdf_path in individual_sdfs:

        residue_name = (
            sdf_path.stem
        )

        run_command(
            [
                "singularity",
                "exec",
                str(conformator_container),
                "python",
                molfile_to_params,
                sdf_path.name,
                "-n",
                residue_name,
                "--keep-names",
                "--long-names",
                "--clobber",
                "--no-pdb",
            ],
            cwd=staging_test_params,
        )

        expected_params = (
            staging_test_params
            / f"{residue_name}.params"
        )

        if not expected_params.is_file():

            raise RuntimeError(
                f"Expected params file was not generated:\n"
                f"  {expected_params}"
            )

        params_names.append(
            expected_params.name
        )

    # --------------------------------------------------------
    # Create Rosetta database control files.
    # --------------------------------------------------------

    (
        staging_test_params
        / "exclude_pdb_component_list.txt"
    ).touch()

    (
        staging_test_params
        / "patches.txt"
    ).touch()

    residue_types_file = (
        staging_test_params
        / "residue_types.txt"
    )

    with open(
        residue_types_file,
        "w",
    ) as handle:

        handle.write(
            "## the atom_type_set and mm-atom_type_set "
            "to be used for the subsequent parameter\n"
        )

        handle.write(
            "TYPE_SET_MODE full_atom\n"
        )

        handle.write(
            "ATOM_TYPE_SET fa_standard\n"
        )

        handle.write(
            "ELEMENT_SET default\n"
        )

        handle.write(
            "MM_ATOM_TYPE_SET fa_standard\n"
        )

        handle.write(
            "ORBITAL_TYPE_SET fa_standard\n"
        )

        handle.write(
            "## Params files\n"
        )

        for params_name in params_names:

            handle.write(
                params_name + "\n"
            )

    # --------------------------------------------------------
    # Remove temporary chemistry files only AFTER all params
    # have been successfully generated.
    # --------------------------------------------------------

    for path in staging_test_params.glob(
        "*.sdf"
    ):
        path.unlink()

    for path in staging_test_params.glob(
        "*.smi"
    ):
        path.unlink()

    # --------------------------------------------------------
    # Validate staging result.
    # --------------------------------------------------------

    if not test_params_complete(
        staging_test_params
    ):

        raise RuntimeError(
            f"Generated test_params failed validation "
            f"for {ligand_name}."
        )

    # --------------------------------------------------------
    # Promote staging directory to completed result.
    # --------------------------------------------------------

    staging_test_params.rename(
        final_test_params
    )

    elapsed = (
        time.time()
        - start_time
    )

    # --------------------------------------------------------
    # Completion marker.
    # --------------------------------------------------------

    completion = {
        **metadata,
        "params_count": len(
            params_names
        ),
        "elapsed_seconds": elapsed,
        "status": "success",
    }

    with open(
        ligand_dir
        / "test_params_complete.json",
        "w",
    ) as handle:

        json.dump(
            completion,
            handle,
            indent=2,
        )

    print(
        f"[SUCCESS] {ligand_name}: "
        f"{len(params_names):,} params files "
        f"in {elapsed:.1f} sec",
        flush=True,
    )

    return {
        "ligand_name": ligand_name,
        "status": "success",
        "params_count": len(
            params_names
        ),
    }


# ============================================================
# Main
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description=(
            "Prepare Rosetta test_params directories from "
            "similarity-search pipeline CSV hits."
        )
    )

    parser.add_argument(
        "--input-csv",
        required=True,
        help=(
            "Similarity pipeline CSV containing "
            "SMILES, ligand_name, source_bz2, Tanimoto."
        ),
    )

    parser.add_argument(
        "--output-root",
        required=True,
        help=(
            "Root directory where numbered ligand directories "
            "will be created."
        ),
    )

    parser.add_argument(
        "--task-index",
        type=int,
        default=None,
        help=(
            "1-based task index. If omitted, use "
            "LSB_JOBINDEX."
        ),
    )

    parser.add_argument(
        "--ligands-per-task",
        type=int,
        default=1,
        help=(
            "Number of CSV hits processed sequentially by "
            "each LSF array element."
        ),
    )

    parser.add_argument(
        "--ligands-per-directory",
        type=int,
        default=100,
        help=(
            "Number of ligand directories per numbered "
            "organizational directory. Default: 100."
        ),
    )

    parser.add_argument(
        "--conformator-container",
        default=(
            "/pi/summer.thyme-umw/"
            "enamine-REAL-2.6billion/"
            "conformator_container.sif"
        ),
        help="Singularity container containing Conformator.",
    )

    parser.add_argument(
        "--conformator-executable",
        default=(
            "/conformator_for_container/"
            "conformator_1.2.1/"
            "conformator"
        ),
        help=(
            "Conformator executable path INSIDE container."
        ),
    )

    parser.add_argument(
        "--molfile-to-params",
        default=(
            "/conformator_for_container/"
            "molfile_to_params.py"
        ),
        help=(
            "molfile_to_params.py path INSIDE container."
        ),
    )

    parser.add_argument(
        "--obabel-executable",
        default="obabel",
        help="Open Babel executable. Default: obabel.",
    )

    parser.add_argument(
        "--max-conformers",
        type=int,
        default=250,
        help=(
            "Requested/default conformer target. Currently "
            "informational unless a Conformator count flag "
            "is added for your installed version."
        ),
    )

    parser.add_argument(
        "--license-env",
        default="CONFORMATOR_LICENSE",
        help=(
            "Environment-variable name containing optional "
            "Conformator license string. The secret itself "
            "is never supplied as an ordinary script argument."
        ),
    )

    parser.add_argument(
        "--overwrite",
        action="store_true",
        help=(
            "Rebuild existing test_params directories."
        ),
    )

    args = parser.parse_args()

    if args.ligands_per_task < 1:
        raise ValueError(
            "--ligands-per-task must be >= 1"
        )

    if args.ligands_per_directory < 1:
        raise ValueError(
            "--ligands-per-directory must be >= 1"
        )

    # --------------------------------------------------------
    # Resolve task index.
    # --------------------------------------------------------

    if args.task_index is None:

        value = os.environ.get(
            "LSB_JOBINDEX"
        )

        if value is None:

            raise RuntimeError(
                "--task-index was not provided and "
                "LSB_JOBINDEX is undefined."
            )

        task_index = int(
            value
        )

    else:

        task_index = (
            args.task_index
        )

    if task_index < 1:

        raise ValueError(
            "task index must be >= 1"
        )

    input_csv = Path(
        args.input_csv
    ).resolve()

    output_root = Path(
        args.output_root
    ).resolve()

    if not input_csv.is_file():

        raise FileNotFoundError(
            f"Input CSV not found: "
            f"{input_csv}"
        )

    output_root.mkdir(
        parents=True,
        exist_ok=True,
    )

    conformator_container = Path(
        args.conformator_container
    )

    if not conformator_container.is_file():

        raise FileNotFoundError(
            f"Conformator container not found: "
            f"{conformator_container}"
        )

    # --------------------------------------------------------
    # Safely retrieve optional license from environment.
    # --------------------------------------------------------

    conformator_license = (
        os.environ.get(
            args.license_env,
            "",
        )
    )

    # --------------------------------------------------------
    # Get this task's records.
    # --------------------------------------------------------

    records = list(
        iter_assigned_records(
            input_csv,
            task_index,
            args.ligands_per_task,
        )
    )

    if not records:

        print(
            f"[INFO] Task {task_index} has no assigned "
            f"ligands. Exiting successfully.",
            flush=True,
        )

        return

    print(
        "============================================================",
        flush=True,
    )

    print(
        f"Test-params preparation task {task_index}",
        flush=True,
    )

    print(
        f"Ligands assigned: {len(records)}",
        flush=True,
    )

    print(
        f"Input CSV:        {input_csv}",
        flush=True,
    )

    print(
        f"Output root:      {output_root}",
        flush=True,
    )

    print(
        "============================================================",
        flush=True,
    )

    # --------------------------------------------------------
    # Process assigned records sequentially.
    #
    # We intentionally do NOT multiprocessing-parallelize
    # ligands inside this job. LSF provides the outer
    # parallelism.
    # --------------------------------------------------------

    successes = 0
    skipped = 0

    for record in records:

        result = prepare_ligand(
            record=record,
            output_root=output_root,
            ligands_per_directory=(
                args.ligands_per_directory
            ),
            conformator_container=(
                conformator_container
            ),
            conformator_executable=(
                args.conformator_executable
            ),
            molfile_to_params=(
                args.molfile_to_params
            ),
            obabel_executable=(
                args.obabel_executable
            ),
            max_conformers=(
                args.max_conformers
            ),
            conformator_license=(
                conformator_license
            ),
            overwrite=(
                args.overwrite
            ),
        )

        if result["status"] == "success":
            successes += 1

        elif result["status"] == "skipped":
            skipped += 1

    print(
        "============================================================",
        flush=True,
    )

    print(
        f"[DONE] Task {task_index}: "
        f"{successes} generated, "
        f"{skipped} skipped",
        flush=True,
    )

    print(
        "============================================================",
        flush=True,
    )


if __name__ == "__main__":
    main()
