# this script works with a placements directory that was dehydrated with
# shrink_placement_pdbs_to_placement_and_surrounding_residues.py,
# using a skeleton pdb to rebuild full pdbs

import os, sys

working_location = os.getcwd()

if len(sys.argv) == 2:
    working_location = sys.argv[1]

os.chdir(working_location)

os.system("tar -xzf placements.tar.gz")
os.chdir("placements")
placements_location = os.getcwd()


def parse_protein_residue_blocks(pdb_path):
    """
    Return:
      res_blocks: dict[int_resseq] -> list[str ATOM lines]
      first_nonatom_prefix: list[str]  # any lines before first ATOM (rare, but safe)
      nonatom_lines: list[str]         # all non-ATOM lines (HETATM, REMARK, footer, etc.)
      saw_ter: bool                    # whether TER was present
    Notes:
      - We treat any line starting with "ATOM" as protein.
      - Residue key is resseq only (single-chain assumption).
    """
    res_blocks = {}
    first_nonatom_prefix = []
    nonatom_lines = []
    saw_atom = False
    saw_ter = False

    cur_resseq = None
    cur_lines = []

    with open(pdb_path, "r") as fh:
        for line in fh:
            if line.startswith("ATOM"):
                saw_atom = True
                resseq = int(line.split()[5])

                if cur_resseq is None:
                    cur_resseq = resseq
                    cur_lines = [line]
                elif resseq == cur_resseq:
                    cur_lines.append(line)
                else:
                    # finalize previous residue
                    res_blocks[cur_resseq] = cur_lines
                    # start new residue
                    cur_resseq = resseq
                    cur_lines = [line]
            else:
                # non-ATOM
                if not saw_atom:
                    first_nonatom_prefix.append(line)
                else:
                    nonatom_lines.append(line)

                if line.startswith("TER"):
                    saw_ter = True

        # finalize last residue if file ended after ATOMs
        if cur_resseq is not None and len(cur_lines) > 0:
            res_blocks[cur_resseq] = cur_lines

    return res_blocks, first_nonatom_prefix, nonatom_lines, saw_ter


# --- parse skeleton protein once ---
skeleton_blocks, skeleton_prefix, skeleton_nonatom, skeleton_saw_ter = parse_protein_residue_blocks("skeleton.pdb")

# For rehydration we typically want the "universe" of residues to be the skeleton's residues
skeleton_resseqs_sorted = sorted(skeleton_blocks.keys())


# --- rehydrate each dehydrated placement ---
for r, d, f in os.walk(placements_location):
    for file in f:
        if r == placements_location and file.endswith(".pdb") and file != "skeleton.pdb":
            print(file)

            placement_blocks, placement_prefix, placement_nonatom, placement_saw_ter = parse_protein_residue_blocks(file)

            out_path = "temp.pdb"
            with open(out_path, "w") as out:

                # 1) Write any header/prefix lines from the placement first (keeps remarks, etc.)
                #    (If you prefer skeleton header, swap these.)
                for line in placement_prefix:
                    out.write(line)

                # 2) Write full protein in sorted residue order once:
                #    placement overrides skeleton where present
                for resseq in skeleton_resseqs_sorted:
                    if resseq in placement_blocks:
                        for line in placement_blocks[resseq]:
                            out.write(line)
                    else:
                        for line in skeleton_blocks[resseq]:
                            out.write(line)

                # 3) Ensure there is a TER after protein (optional but usually correct)
                #    If your placement_nonatom already includes TER/END/footer, you can skip this.
                #    We'll only add TER if the placement didn't already have one in its non-ATOM.
                if not any(l.startswith("TER") for l in placement_nonatom):
                    out.write("TER\n")

                # 4) Write all non-ATOM lines from the placement (ligand HETATM + footer comments)
                #    BUT: avoid duplicating TER/END if you injected one above.
                for line in placement_nonatom:
                    # If we wrote TER already, skip placement TER lines to avoid duplicates
                    if line.startswith("TER") and not any(l.startswith("TER") for l in placement_nonatom):
                        continue
                    out.write(line)

            os.system(f"mv {out_path} {file}")


# recompress and clean up
os.chdir(working_location)
os.system("tar -czf placements.tar.gz placements")
os.system("rm -drf placements")
