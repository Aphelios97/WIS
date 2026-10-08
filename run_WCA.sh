#!/bin/bash

# ============================================================
# Run WCA
# Usage:
#   ./run_WCA.sh dataset delta
#
# Example:
#   ./run_WCA.sh email-Enron-full 86400000
# ============================================================

if [ $# -ne 2 ]; then
    echo "Usage: $0 <dataset> <delta>"
    exit 1
fi

DATASET=$1
DELTA=$2

echo "Running WCA..."
echo "Dataset: $DATASET"
echo "delta: $DELTA"

./WCA "$DATASET" "$DELTA"
