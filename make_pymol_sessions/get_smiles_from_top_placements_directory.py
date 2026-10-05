#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path
from collections import defaultdict

from rdkit import Chem


# Common HETATM residues that usually are not the ligand of interest.
IGNORE_RESNAMES = {
    "HOH", "WAT", "DOD",
    "NA", "K", "CL", "CA", "MG", "ZN", "MN", "FE", "CU", "CO",
    "SO4", "PO4",
    "GOL", "EDO", "PEG",
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract ligand HETATM records from PDB files and convert them to SMILES."
    )

    parser.add_argument(
        "input_dir",
        help="Directory containing PDB files."
    )

    parser.add_argument(
        "-o",
        "--output",
        default="pdb_ligand_smiles.csv",
        help="Output CSV file (default: pdb_ligand_smiles.csv)."
    )

    parser.add_argument(
        "--recursive",
        action="store_true",
        help="Search recursively for PDB files."
    )

    return parser.parse_args()


def collect_hetatm_residues(pdb_file):
    """
    Return HETATM records grouped by:
        (resname, chain, residue_number, insertion_code)
    """

    residues = defaultdict(list)

    with open(pdb_file, "r") as f:
        for line in f:
            if not line.startswith("HETATM"):
                continue

            # PDB fixed-width fields
            resname = line[17:20].strip()
            chain = line[21].strip()
            resid = line[22:26].strip()
            icode = line[26].strip()

            if resname.upper() in IGNORE_RESNAMES:
                continue

            key = (resname, chain, resid, icode)
            residues[key].append(line.rstrip("\n"))

    return residues


def ligand_block_from_lines(lines):
    """
    Construct a small standalone PDB block for RDKit.
    """
    return "\n".join(lines) + "\nEND\n"


def pdb_lines_to_mol(lines):
    """
    Convert HETATM lines to an RDKit molecule.

    proximityBonding=True lets RDKit infer bonds from coordinates.
    """

    block = ligand_block_from_lines(lines)

    try:
        mol = Chem.MolFromPDBBlock(
            block,
            sanitize=True,
            removeHs=False,
            proximityBonding=True,
        )

        if mol is not None:
            return mol

    except Exception:
        pass

    # Second attempt without initial sanitization
    try:
        mol = Chem.MolFromPDBBlock(
            block,
            sanitize=False,
            removeHs=False,
            proximityBonding=True,
        )

        if mol is None:
            return None

        Chem.SanitizeMol(mol)
        return mol

    except Exception:
        return None


def heavy_atom_count(mol):
    return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1)


def organic_atom_count(mol):
    """
    Count typical organic atoms.
    Useful for distinguishing a ligand from metal ions, etc.
    """
    organic_atomic_numbers = {6, 7, 8, 9, 15, 16, 17, 35, 53}

    return sum(
        1
        for atom in mol.GetAtoms()
        if atom.GetAtomicNum() in organic_atomic_numbers
    )


def find_ligand(pdb_file):
    """
    Find the most likely ligand in a PDB.

    If multiple HETATM residues exist, choose the residue with
    the largest number of organic heavy atoms.
    """

    residues = collect_hetatm_residues(pdb_file)

    candidates = []

    for key, lines in residues.items():

        mol = pdb_lines_to_mol(lines)

        if mol is None:
            continue

        heavy = heavy_atom_count(mol)
        organic = organic_atom_count(mol)

        # Require at least one carbon to avoid selecting ions,
        # sulfates, etc.
        has_carbon = any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())

        if not has_carbon:
            continue

        candidates.append(
            {
                "key": key,
                "mol": mol,
                "lines": lines,
                "heavy_atoms": heavy,
                "organic_atoms": organic,
            }
        )

    if not candidates:
        return None

    # Largest organic molecule is assumed to be the ligand.
    candidates.sort(
        key=lambda x: (x["organic_atoms"], x["heavy_atoms"]),
        reverse=True
    )

    return candidates[0]


def main():

    args = parse_args()

    input_dir = Path(args.input_dir)

    if not input_dir.is_dir():
        raise SystemExit(f"ERROR: directory does not exist: {input_dir}")

    if args.recursive:
        pdb_files = sorted(input_dir.rglob("*.pdb"))
    else:
        pdb_files = sorted(input_dir.glob("*.pdb"))

    print(f"Found {len(pdb_files):,} PDB files")

    rows = []

    for i, pdb_file in enumerate(pdb_files, start=1):

        ligand = find_ligand(pdb_file)

        if ligand is None:
            print(f"[WARNING] No ligand found: {pdb_file.name}")

            rows.append(
                {
                    "source_pdb": pdb_file.name,
                    "resname": "",
                    "chain": "",
                    "resid": "",
                    "heavy_atoms": "",
                    "smiles": "",
                }
            )

            continue

        resname, chain, resid, icode = ligand["key"]

        # Remove explicit hydrogens before generating SMILES
        mol = Chem.RemoveHs(ligand["mol"])

        try:
            smiles = Chem.MolToSmiles(
                mol,
                canonical=True,
                isomericSmiles=True,
            )
        except Exception:
            smiles = ""

        rows.append(
            {
                "source_pdb": pdb_file.name,
                "resname": resname,
                "chain": chain,
                "resid": resid + icode,
                "heavy_atoms": ligand["heavy_atoms"],
                "smiles": smiles,
            }
        )

        print(
            f"[{i:,}/{len(pdb_files):,}] "
            f"{pdb_file.name} -> "
            f"{resname} {chain}:{resid}{icode} -> "
            f"{smiles}"
        )

    with open(args.output, "w", newline="") as f:

        fieldnames = [
            "source_pdb",
            "resname",
            "chain",
            "resid",
            "heavy_atoms",
            "smiles",
        ]

        writer = csv.DictWriter(f, fieldnames=fieldnames)

        writer.writeheader()
        writer.writerows(rows)

    print()
    print(f"Wrote {len(rows):,} entries to:")
    print(args.output)


if __name__ == "__main__":
    main()
