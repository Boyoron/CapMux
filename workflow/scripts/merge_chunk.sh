#!/usr/bin/env bash
set -euo pipefail

# merge_chunk.sh: Merge chunk-level FASTQs produced by capseq_demux.py into per-sample FASTQs
# Usage: merge_chunk.sh INPUT_DIR OUTPUT_DIR [threads]
#
# Expected inputs:
#   1) <chunk>_<sample>_cdna.fastq
#      <chunk>_<sample>_bc.fastq
#      where <chunk> itself may contain underscores, typically ending with _S<number>
#   2) <sample>_cdna.fastq / <sample>_bc.fastq  (no chunk prefix)
#
# After successful merge, INPUT_DIR is deleted.

if [[ $# -lt 2 || $# -gt 3 ]]; then
    echo "Usage: $0 INPUT_DIR OUTPUT_DIR [threads]" >&2
    exit 1
fi

INPUT_DIR="$1"
OUTPUT_DIR="$2"
THREADS="${3:-4}"

if ! command -v pigz >/dev/null 2>&1; then
    echo "ERROR: pigz not found in PATH" >&2
    exit 1
fi

if [[ ! -d "$INPUT_DIR" ]]; then
    echo "ERROR: input directory '$INPUT_DIR' does not exist" >&2
    exit 1
fi

# safety check
case "$INPUT_DIR" in
    ""|"/")
        echo "Refusing to rm -rf '$INPUT_DIR'" >&2
        exit 1
        ;;
esac

mkdir -p "$OUTPUT_DIR"
# use absolute OUTPUT_DIR
OUTPUT_DIR="$(cd "$OUTPUT_DIR" && pwd)"

pushd "$INPUT_DIR" >/dev/null
shopt -s nullglob

# collect candidate inputs
fastq_files=( *.fastq )

if (( ${#fastq_files[@]} == 0 )); then
    echo "No .fastq files found in $INPUT_DIR" >&2
    popd >/dev/null
    exit 1
fi

# groups[key] = newline separated list of files to merge
declare -A groups

# parse sample/read
for f in "${fastq_files[@]}"; do
    base="${f%.fastq}"

    if [[ "$base" != *_* ]]; then
        echo "Skipping $f: does not look like <...>_<read>.fastq" >&2
        continue
    fi

    read_part="${base##*_}"   # bc or cdna
    rest="${base%_*}"         # everything before read

    if [[ "$read_part" != "bc" && "$read_part" != "cdna" ]]; then
        continue
    fi

    sample_part=""

    # case a): chunk ends with _S<number> (bcl2fastq-style)
    if [[ "$rest" =~ ^(.*_S[0-9]+)_(.*)$ ]]; then
        sample_part="${BASH_REMATCH[2]}"
    else
        # case b): "<chunk>_<sample>"
        if [[ "$rest" == *_* ]]; then
            sample_part="${rest#*_}"
        else
            # case c): no chunk
            sample_part="$rest"
        fi
    fi

    key="${sample_part}_${read_part}"

    # appending file to the group's list
    if [[ -z "${groups[$key]+x}" ]]; then
        groups["$key"]="$f"
    else
        groups["$key"]+=$'\n'"$f"
    fi
done

if (( ${#groups[@]} == 0 )); then
    echo "No mergeable files found (expected *_bc.fastq or *_cdna.fastq) in $INPUT_DIR" >&2
    popd >/dev/null
    exit 1
fi

# merge each group
for key in "${!groups[@]}"; do
    out="${OUTPUT_DIR%/}/${key}.fastq.gz"

    # read newline separated file list into array and sort
    IFS=$'\n' read -r -d '' -a files < <(printf '%s\0' "${groups[$key]}")
    # sort
    mapfile -t files < <(printf "%s\n" "${files[@]}" | LC_ALL=C sort)

    echo "Merging + compressing ${#files[@]} file(s) into $out" >&2
    pigz -p "$THREADS" -c "${files[@]}" > "$out"
done

popd >/dev/null

echo "Removing input directory: $INPUT_DIR" >&2
rm -rf -- "$INPUT_DIR"

echo "Done."
