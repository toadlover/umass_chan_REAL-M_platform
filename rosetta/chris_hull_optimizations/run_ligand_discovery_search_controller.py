#!/usr/bin/env python3
"""
HPC‑optimized controller for launching Rosetta ligand discovery jobs on LSF.

Features:
- Clean, deterministic job list generation
- LSF‑safe array submission
- Automatic tmp-space reservation via rusage[tmp=X]
- Scratch-size sanity checking
- Drop‑in replacement for the user's original controller
"""

import os
import sys
from pathlib import Path
import subprocess

starting_location = os.getcwd()

def run(cmd: str):
    """Run a shell command safely."""
    subprocess.run(cmd, shell=True, check=True)


# ---------------------------------------------------------------------
# Parse controller arguments
# ---------------------------------------------------------------------
if len(sys.argv) < 8:
    print(
        "Usage: run_ligand_discovery_search_controller.py "
        "<target_pdb> <anchor_list> <motifs_file> <discovery_root> "
        "<atr> <rep> <ddg> [extra_args_file]"
    )
    sys.exit(1)

target_pdb = Path(sys.argv[1])
anchor_list = sys.argv[2]        # e.g. "63,87,96,179"
motifs_file = Path(sys.argv[3])
discovery_root = Path(sys.argv[4])
atr, rep, ddg = sys.argv[5:8]
extra_args_file = Path(sys.argv[8]) if len(sys.argv) >= 9 else None


#look down the discovery root for all test_params directories

test_params_directories = []

for r,d,f in os.walk(discovery_root):
    for dire in d:
        if dire == "test_params":
            test_params_directories.append(r + "/" + dire + "/")

# ---------------------------------------------------------------------
# Build job list (one job per anchor residue)
# ---------------------------------------------------------------------
anchors = [a.strip() for a in anchor_list.split(",") if a.strip()]

#make a list that holds all job lists by name and their corresponding lengths
#this way, we can handle larger directory clusters where we would exceed the umass hpc array job size limit of 10,000 jobs (safely cap individual arrays at 8,000)
all_joblist_list = []
#count the number of jobs in the joblist
joblist_job_counter = 0
#count the number of joblist files
joblist_file_counter = 0

#working joblist file name
working_joblist_file = "joblist_" + str(joblist_file_counter) + ".txt"

joblist_path = Path(starting_location + "/" + working_joblist_file)
#with joblist_path.open("w") as f:

#do an implicit check to see if there is a compiled raw scores file from round 1 discovery that is in the discovery root, which would have been used to create discovery directories fro round 2 discovery
#the purpose of this is to be smartr about expanded conformer discovery
#if there are multiple anchor residues being used, we only want to run expanded conformer discovery on the anchor(s) that led to placements in the expanded conformer set
#this first step is to perform a possible read-in of the corresponding compiled raw scores csv if it is applicable
#if we do want to do this, we will read in all ligands from the csv by name and pair them with the anchor that they are listed with
secondary_discovery_ligands_with_anchors = {}
#only even consider if there is more than one anchor, otherwise we implicitly will have all ligands only run expanded conformers on the single anchor anyway
if len(anchors > 1):
    #runthrough and see if there is a raw scores file in the directory
    #makign this choice instead of adding another argument to not have to deal with changing the entire pipeline architecture
    #if you are runnign secondary discovery, there should eb a raw scores file in the discovery space anyway, since it had to be used to make the test params directories
    #if no usable csv is found, all anchor residues will just implicitly be used instead
    for r,d,f in os.walk(discovery_root):
        for file in f:
            if r == "discovery_root" and file.endswith(".csv"):
                #read the csv and extract all ligands with their corresponding anchor
                placements_file = open(r + "/" + file, "r")

                for line in placements_file.readlines():
                    #primary filter to filter out potential bad lines, lines must have "res" and ".pdb" in the line, line.split(",") must also be at leat 7 entries
                    if len(line.split(",")) < 7:
                        continue
                    if "/res" not in line and ".pdb," not in line:
                        continue
                    #derive the anchor and ligand name
                    #anchor is derived between "/res" and the following udnerscore
                    anchor = line.split("/res")[1].split("_")[0]
                    #confirm that this is in the anchors list
                    if anchor not in anchors:
                        continue
                    #next, derive the ligand, which should be the 3rd to last entry when splitting by underscores when first splitting by split(".pdb")[0]
                    this_ligand = line.split(".pdb")[0].split("_")[-3]

                    #pair the ligand and anchor
                    #if the ligands is already listed, add the anchor
                    if this_ligand not in secondary_discovery_ligands_with_anchors.keys():
                        secondary_discovery_ligands_with_anchors[this_ligand] = [anchor]
                    else:
                        if anchor in secondary_discovery_ligands_with_anchors[this_ligand]:
                            continue
                        else:
                            #append the anchor to the list in the dictionary if it was not already present
                            secondary_discovery_ligands_with_anchors[this_ligand].append(anchor)

f = joblist_path.open("w")

for tp_dire in test_params_directories:
    for anchor in anchors:
        line = (
            f"{target_pdb} "
            f"{anchor} "
            f"{motifs_file} "
            f"{tp_dire} "
            f"{atr} {rep} {ddg}"
        )

        
        #skip the anchor if it is not in the secondary list, if there is info on the ligand anchor (for secondary discovery)
        if len(secondary_discovery_ligands_with_anchors) > 0:
            #determine what the ligand here should be from the tp_dire
            #ligand should be the 2nd to last entry when splitting by slashes
            #to account for where there may or may not be a trailing slash, split by /test_params first and then it is the first thing after teh last slash
            cur_lig = tp_dire.split("/test_params")[0].split("/")[-1]

            if cur_lig in secondary_discovery_ligands_with_anchors.keys():
                if anchor not in secondary_discovery_ligands_with_anchors[cur_lig]:
                    #skip entry if the ligand is in the ligands list by we do not have a placement from it for the anchor
                    continue


        if extra_args_file:
            line += f" {extra_args_file}"
        f.write(line + "\n")

        joblist_job_counter = joblist_job_counter + 1

        #if capping the limit of 8,000
        if joblist_job_counter == 8000:

            #append the file and count to the list
            all_joblist_list.append([joblist_path,joblist_job_counter])

            #change the working joblist file name and full path
            joblist_file_counter = joblist_file_counter + 1
            working_joblist_file = "joblist_" + str(joblist_file_counter) + ".txt"
            joblist_path = Path(starting_location + "/" + working_joblist_file)

            #reset the job counter
            joblist_job_counter = 0

            #open a new write stream to keep going in a new file
            f.close()
            f = joblist_path.open("w")

#handle the final dangling job if there are any remaining jobs (which there likely are and this handles if we never hit the cap)
if joblist_job_counter > 0:
    f.close()
    all_joblist_list.append([joblist_path,joblist_job_counter])

num_jobs = len(anchors) * len(test_params_directories)
print(f"Prepared job list with {num_jobs} jobs.")


# ---------------------------------------------------------------------
# Ensure logs directory exists
# ---------------------------------------------------------------------
logs = Path("logs")
logs.mkdir(exist_ok=True)


# ---------------------------------------------------------------------
# SCRATCH SIZE COMMENTARY + SANITY CHECK
#
# The worker script performs all heavy I/O inside a loopback-mounted
# filesystem located in /tmp. The wrapper script creates this filesystem
# with a size determined by SCRATCH_SIZE_GB.
#
# To prevent oversubscription of /tmp on compute nodes, we request the
# same amount of space from LSF using:
#
#     -R "rusage[tmp=X]"
#
# where X is in megabytes.
#
# This ensures:
#   • LSF will only dispatch the job to nodes with enough free /tmp space.
#   • Multiple jobs cannot overrun a node's local disk.
#   • The wrapper and scheduler always agree on the required scratch size.
#
# Users may override the default (10 GB) by setting:
#
#     export SCRATCH_SIZE_GB=20
#
# before running this controller.
#
# SANITY CHECK:
#   We ensure the scratch size is large enough to hold:
#     • The motif file
#     • Temporary Rosetta outputs
#     • Placement PDBs
#   Minimum recommended scratch = max(5 GB, 2 × motif file size)
# ---------------------------------------------------------------------

scratch_gb = int(os.environ.get("SCRATCH_SIZE_GB", 20))
tmp_mb = scratch_gb * 1024  # LSF expects MB

motif_size_gb = motifs_file.stat().st_size / (1024**3)
min_required_gb = max(5, motif_size_gb * 2)

if scratch_gb < min_required_gb:
    print(
        f"\nERROR: Requested SCRATCH_SIZE_GB={scratch_gb} GB is too small.\n"
        f"Motif file alone is {motif_size_gb:.2f} GB.\n"
        f"Minimum recommended scratch is {min_required_gb:.1f} GB.\n"
        f"Please increase SCRATCH_SIZE_GB and resubmit.\n"
    )
    sys.exit(1)

if motif_size_gb > scratch_gb * 0.25:
    print(
        f"WARNING: Motif file ({motif_size_gb:.2f} GB) is more than 25% "
        f"of requested scratch ({scratch_gb} GB)."
    )


# ---------------------------------------------------------------------
# Submit LSF job array
# ---------------------------------------------------------------------
wrapper = Path("/pi/summer.thyme-umw/enamine-REAL-2.6billion/umass_chan_REAL-M_platform/rosetta/chris_hull_optimizations/lsf_scratch_wrapper.sh").resolve()

if not wrapper.exists():
    print("ERROR: lsf_scratch_wrapper.sh not found.")
    sys.exit(1)

#iterate over every job list
for joblist in all_joblist_list:
    bsub_cmd = (
        f"bsub "
        f"-J rosetta_ld[1-{joblist[1]}] "
        f"-R \"rusage[mem=10000,tmp={tmp_mb}]\" "
        f"-q long "
        f"-W 96:00 "
        f"-o logs/%J_%I.out "
        f"-e logs/%J_%I.err "
    #    f"bash {wrapper} $(sed -n \"\\$LSB_JOBINDEX\"p " + starting_location + "/joblist.txt)"
        f"bash {wrapper} " + str(joblist[0])
    )

    print("\nSubmitting LSF job array:")
    print(bsub_cmd, "\n")

    run(bsub_cmd)

