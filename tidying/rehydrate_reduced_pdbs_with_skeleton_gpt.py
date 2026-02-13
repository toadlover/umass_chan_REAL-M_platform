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
    Returns:
      res_blocks: dict[int] -> list[str]  # ATOM lines grouped by resseq
      prefix_lines: list[str]             # non-ATOM lines before first ATOM
      nonatom_lines: list[str]            # non-ATOM lines after first ATOM (HETATM, REMARK, footer, TER, END, etc.)
    Assumes single chain and uses resseq = int(line.split()[5]).
    """
    res_blocks = {}
    prefix_lines = []
    nonatom_lines = []

    saw_atom = False
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
                if not saw_atom:
                    prefix_lines.append(line)
                else:
                    nonatom_lines.append(line)

        # finalize last residue at EOF
        if cur_resseq is not None and cur_lines:
            res_blocks[cur_resseq] = cur_lines

    return res_blocks, prefix_lines, nonatom_lines


# Parse skeleton once
skeleton_blocks, skeleton_prefix, skeleton_nonatom = parse_protein_residue_blocks("skeleton.pdb")


for r, d, f in os.walk(placements_location):
    for file in f:
        if r == placements_location and file.endswith(".pdb") and file != "skeleton.pdb":
            print(file)

            placement_blocks, placement_prefix, placement_nonatom = parse_protein_residue_blocks(file)

            # --- FIX: union of residue indices (handles terminal residues present only in placement) ---
            all_resseqs_sorted = sorted(set(skeleton_blocks.keys()) | set(placement_blocks.keys()))

            out_path = "temp.pdb"
            with open(out_path, "w") as out:

                # Write header/prefix from placement (keeps remarks at top if any)
                for line in placement_prefix:
                    out.write(line)

                # Write full protein in order: placement overrides skeleton
                for resseq in all_resseqs_sorted:
                    if resseq in placement_blocks:
                        for line in placement_blocks[resseq]:
                            out.write(line)
                    elif resseq in skeleton_blocks:
                        for line in skeleton_blocks[resseq]:
                            out.write(line)
                    else:
                        # Shouldn't happen because we're iterating the union
                        pass

                # Write the non-ATOM lines from placement (ligand + footer + END, etc.)
                # (This preserves your ligand and any comment block.)
                for line in placement_nonatom:
                    out.write(line)

            os.system(f"mv {out_path} {file}")


# recompress and clean up
os.chdir(working_location)
os.system("tar -czf placements.tar.gz placements")
os.system("rm -drf placements")
