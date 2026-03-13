#!/bin/bash

set -euo pipefail

exec > "${snakemake_log[0]}" 2>&1

BC="${snakemake_input[bc]}"      # R1 = barcode read (bc)
CDNA="${snakemake_input[cdna]}"  # R2 = cDNA read (cdna)

TSOseq="${snakemake_params[tso_seq]}"

cutadapt \
    -G "${TSOseq}" \
    -m :40 \
    -j "${snakemake[threads]}" \
    -o "${snakemake_output[bc_trimmed]}" \
    -p "${snakemake_output[cdna_trimmed]}" \
    "${BC}" \
    "${CDNA}" \
    > "${snakemake_output[report]}"
