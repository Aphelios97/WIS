#!/bin/bash

# ============================================================
# Run WIS
# Usage:
#   ./run_WIS.sh dataset S b c delta
#
# Example:
#   ./run_WIS.sh email-Eu-full 100 2 2.0 86400
# ============================================================

if [ $# -ne 5 ]; then
    echo "Usage: $0 <dataset> <S> <b> <c> <delta>"
    exit 1
fi

DATASET=$1
S=$2
B=$3
C=$4
DELTA=$5

echo "Running WIS..."
echo "Dataset: $DATASET"
echo "Sample Times: $S"
echo "b: $B"
echo "c: $C"
echo "delta: $DELTA"

./WIS_modified "$DATASET" "$S" "$B" "$C" "$DELTA"
