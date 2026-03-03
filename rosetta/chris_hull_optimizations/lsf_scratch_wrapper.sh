#!/bin/bash
set -euo pipefail

###############################################################################
# LSF SCRATCH WRAPPER
#
# This script:
#   • Determines a unique scratch directory for each job (array‑safe)
#   • Creates a loopback filesystem of size SCRATCH_SIZE_GB
#   • Mounts it at $SCRATCH_DIR
#   • Exports SCRATCH_DIR for the worker script
#   • Runs the worker script inside the scratch filesystem
#   • Cleans up automatically on exit
#
# It is invoked automatically by the controller:
#   bash lsf_scratch_wrapper.sh <args from joblist.txt>
###############################################################################


###############################################################################
# CONFIGURATION
###############################################################################

# Default scratch size (GB) if not provided by user
SCRATCH_SIZE_GB="${SCRATCH_SIZE_GB:-10}"


###############################################################################
# DETERMINE SCRATCH DIRECTORY
#
# If SCRATCH_DIR is not set, create a unique path under /tmp:
#   /tmp/$USER/$LSB_JOBID.scratch
#   /tmp/$USER/$LSB_JOBID_$LSB_JOBINDEX.scratch   (for array jobs)
###############################################################################

if [ -z "${SCRATCH_DIR:-}" ]; then
    if [ -n "${LSB_JOBINDEX:-}" ]; then
        SCRATCH_DIR="/tmp/${USER}/${LSB_JOBID}_${LSB_JOBINDEX}.scratch"
    else
        SCRATCH_DIR="/tmp/${USER}/${LSB_JOBID}.scratch"
    fi
fi

SCRATCH_FILE="${SCRATCH_DIR}.img"
SCRATCH_MNT="${SCRATCH_DIR}"


###############################################################################
# CHECK /tmp SPACE
#
# Ensure the node has enough free space before attempting to create the
# loopback filesystem. LSF also enforces this via rusage[tmp=X], but this
# check protects against unexpected conditions.
###############################################################################

TMP_FREE=$(df --output=avail -k /tmp | tail -1)
REQUIRED=$((SCRATCH_SIZE_GB * 1024 * 1024))   # GB → KB

if [ "$TMP_FREE" -lt "$REQUIRED" ]; then
    echo "ERROR: Not enough free space in /tmp for ${SCRATCH_SIZE_GB} GB scratch."
    echo "       Required: ${REQUIRED} KB"
    echo "       Available: ${TMP_FREE} KB"
    exit 1
fi


###############################################################################
# CREATE LOOPBACK SCRATCH FILESYSTEM
###############################################################################

# Ensure parent directory exists
mkdir -p "$(dirname "$SCRATCH_FILE")"

# Create sparse file
truncate -s "${SCRATCH_SIZE_GB}G" "$SCRATCH_FILE"

# Format as ext4
mkfs.ext4 -F "$SCRATCH_FILE" >/dev/null

# Create mount point
mkdir -p "$SCRATCH_MNT"

# Mount loopback
mount -o loop "$SCRATCH_FILE" "$SCRATCH_MNT"


###############################################################################
# CLEANUP HANDLER
#
# Ensures scratch is unmounted and removed even if the job crashes.
###############################################################################

cleanup() {
    umount "$SCRATCH_MNT" 2>/dev/null || true
    rm -f "$SCRATCH_FILE" 2>/dev/null || true
    rmdir "$SCRATCH_MNT" 2>/dev/null || true
}
trap cleanup EXIT


###############################################################################
# RUN WORKER SCRIPT INSIDE SCRATCH
###############################################################################

export SCRATCH_DIR="$SCRATCH_MNT"

# Pass all arguments through to the worker script
python run_ligand_discovery_search.py "$@"

