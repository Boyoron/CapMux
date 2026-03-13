#!/usr/bin/env bash
set -euo pipefail


run_dir="${snakemake_input[run_dir]}"
out_dir="${snakemake_output[out_dir]}"
sample_sheet="${snakemake_input[sample_sheet]}"
threads="${snakemake[threads]}"

use_bases_mask="${snakemake_params[use_bases_mask]}"
demux_index="${snakemake_params[demux_index]}"

log="${snakemake_log[0]}"

# --use-bases-mask selection
if [[ -z "${use_bases_mask}" ]]; then
    case "${demux_index}" in
        i7|i5)
            bases_mask="y*,I*,y*"
            ;;
        i7_i5)
            bases_mask="y*,I*,I*,y*"
            ;;
        *)
            echo "ERROR: demux_index must be one of: i5, i7, i7_i5 (got: '${demux_index}')" >&2
            exit 1
            ;;
    esac
else
    bases_mask="${use_bases_mask}"
fi

bcl2fastq \
    --runfolder-dir "${run_dir}" \
    --output-dir "${out_dir}" \
    --sample-sheet "${sample_sheet}" \
    --mask-short-adapter-reads 0 \
    --minimum-trimmed-read-length 0 \
    --use-bases-mask "${bases_mask}" \
    --no-lane-splitting \
    --create-fastq-for-index-reads \
    --barcode-mismatches 1 \
    --processing-threads "${threads}" \
    &> "${log}"
