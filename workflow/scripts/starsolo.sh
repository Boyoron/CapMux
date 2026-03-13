#!/usr/bin/env bash
set -euo pipefail

mkdir -p "$(dirname "${snakemake_log[0]}")"
exec > "${snakemake_log[0]}" 2>&1

index="${snakemake_input[index]}"
cdna="${snakemake_input[cdna]}"
bc="${snakemake_input[bc]}"

threads="${snakemake[threads]}"
out_prefix="${snakemake_params[out_prefix]}"

mkdir -p "$(dirname "${out_prefix}")"
features="${snakemake_params[features]}"
demux_by="${snakemake_params[demux_by]}"   # "bc1" or "index"
cb_match="${snakemake_params[cb_match]}"

bc1_start="${snakemake_params[bc1_start]}"; bc1_len="${snakemake_params[bc1_len]}"
bc2_start="${snakemake_params[bc2_start]}"; bc2_len="${snakemake_params[bc2_len]}"
bc3_start="${snakemake_params[bc3_start]}"; bc3_len="${snakemake_params[bc3_len]}"
umi1_start="${snakemake_params[umi1_start]}"; umi1_len="${snakemake_params[umi1_len]}"
umi2_start="${snakemake_params[umi2_start]}"; umi2_len="${snakemake_params[umi2_len]}"

wl_bc1="${snakemake_params[wl_bc1]}"
wl_bc2="${snakemake_params[wl_bc2]}"
wl_bc3="${snakemake_params[wl_bc3]}"
wl_index="${snakemake_params[wl_index]}"

is_uint() { [[ "$1" =~ ^[0-9]+$ ]]; }

validate_segment() {
  local name="$1" s="$2" l="$3"

  if ! is_uint "$s" || ! is_uint "$l"; then
    echo "ERROR: segment '$name' has non-integer coords: start='$s' len='$l'" >&2
    exit 2
  fi
  if [[ "$s" -eq 0 && "$l" -eq 0 ]]; then
    return 1
  fi
  if [[ "$s" -eq 0 || "$l" -eq 0 ]]; then
    echo "ERROR: segment '$name' must be disabled as start=0,len=0 OR enabled with start>=1,len>=1 (got start=$s len=$l)" >&2
    exit 2
  fi
  if [[ "$s" -lt 1 || "$l" -lt 1 ]]; then
    echo "ERROR: segment '$name' enabled requires start>=1,len>=1 (got start=$s len=$l)" >&2
    exit 2
  fi
  return 0
}

first_seq_len() {
  awk 'NF && $0 !~ /^#/ {gsub(/\r/, ""); print length($0); exit}' "$1"
}

require_file() {
  local label="$1" path="$2"
  if [[ -z "$path" ]]; then
    echo "ERROR: missing path for ${label} whitelist in config.yaml" >&2
    exit 2
  fi
  if [[ ! -f "$path" ]]; then
    echo "ERROR: ${label} whitelist file not found: $path" >&2
    exit 2
  fi
}

USED_BC=0
USED_UMI=0
UMI_LEN_TOTAL=0

declare -a CB_WL=()
declare -a CB_LEN=()
declare -a CB_POS=()

add_cb_segment() {
  local name="$1" seg_len="$2" wl="$3"
  require_file "$name" "$wl"
  got_len="$(first_seq_len "$wl" || true)"
  if [[ -z "$got_len" ]]; then
    echo "ERROR: whitelist '$wl' is empty" >&2
    exit 2
  fi
  if [[ "$got_len" -ne "$seg_len" ]]; then
    echo "ERROR: whitelist '$wl' has sequence length $got_len, but config expects ${seg_len} for ${name}" >&2
    exit 2
  fi
  CB_WL+=("$wl")
  CB_LEN+=("$seg_len")
  USED_BC=$((USED_BC+1))
}

if validate_segment bc1 "$bc1_start" "$bc1_len"; then add_cb_segment "bc1" "$bc1_len" "$wl_bc1"; fi
if validate_segment bc2 "$bc2_start" "$bc2_len"; then add_cb_segment "bc2" "$bc2_len" "$wl_bc2"; fi
if validate_segment bc3 "$bc3_start" "$bc3_len"; then add_cb_segment "bc3" "$bc3_len" "$wl_bc3"; fi

if validate_segment umi1 "$umi1_start" "$umi1_len"; then UMI_LEN_TOTAL=$((UMI_LEN_TOTAL + umi1_len)); USED_UMI=$((USED_UMI+1)); fi
if validate_segment umi2 "$umi2_start" "$umi2_len"; then UMI_LEN_TOTAL=$((UMI_LEN_TOTAL + umi2_len)); USED_UMI=$((USED_UMI+1)); fi

if [[ "$USED_BC" -lt 1 ]]; then
  echo "ERROR: At least one of bc1/bc2/bc3 must be enabled." >&2
  exit 2
fi
if [[ "$USED_UMI" -lt 1 ]]; then
  echo "ERROR: At least one of umi1/umi2 must be enabled." >&2
  exit 2
fi

if [[ "$demux_by" != "bc1" && "$demux_by" != "index" ]]; then
  echo "ERROR: demux_by must be 'bc1' or 'index' (got '$demux_by')" >&2
  exit 2
fi

index_len=0
if [[ "$demux_by" == "bc1" ]]; then
  if ! validate_segment bc1 "$bc1_start" "$bc1_len"; then
    echo "ERROR: demux_by=bc1 requires bc1 to be enabled (it is used for demultiplexing)." >&2
    exit 2
  fi
  require_file "index" "$wl_index"
  index_len="$(first_seq_len "$wl_index" || true)"
  if [[ -z "$index_len" ]]; then
    echo "ERROR: index whitelist '$wl_index' is empty" >&2
    exit 2
  fi
  if ! is_uint "$index_len" || [[ "$index_len" -lt 1 ]]; then
    echo "ERROR: could not determine index length from '$wl_index' (got '$index_len')" >&2
    exit 2
  fi
  CB_WL+=("$wl_index")
  CB_LEN+=("$index_len")
fi

pos=0
for l in "${CB_LEN[@]}"; do
  start0=$pos
  end0=$((pos + l - 1))
  CB_POS+=("0_${start0}_0_${end0}")
  pos=$((pos + l))
done

UMI_START0=$pos
UMI_END0=$((pos + UMI_LEN_TOTAL - 1))
UMI_POS="0_${UMI_START0}_0_${UMI_END0}"
expected_total=$((pos + UMI_LEN_TOTAL))

set +o pipefail
actual_len="$(zcat "$bc" | awk 'NR==2 {gsub(/\r/, ""); print length($0); exit}')"
set -o pipefail
if [[ -z "$actual_len" ]]; then
  echo "ERROR: could not read barcode FASTQ: $bc" >&2
  exit 2
fi
if [[ "$actual_len" -ne "$expected_total" ]]; then
  echo "ERROR: barcode read length mismatch. Expected $expected_total (CB+UMI stitched), got $actual_len from $bc" >&2
  echo "       demux_by=$demux_by, CB lens=${CB_LEN[*]}, UMI total len=$UMI_LEN_TOTAL" >&2
  exit 2
fi

echo "INFO: demux_by=$demux_by"
echo "INFO: soloCBposition: ${CB_POS[*]}"
echo "INFO: soloUMIposition: ${UMI_POS}"
echo "INFO: soloCBwhitelist: ${CB_WL[*]}"
echo "INFO: soloCBmatchWLtype: ${cb_match}"

STAR \
  --genomeDir "${index}" \
  --readFilesIn "${cdna}" "${bc}" \
  --soloCBwhitelist "${CB_WL[@]}" \
  --runThreadN "${threads}" \
  --outFileNamePrefix "${out_prefix}" \
  --readFilesCommand zcat \
  --runDirPerm All_RWX \
  --outReadsUnmapped None \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMattributes NH HI nM AS CR UR CB UB sS sQ sM GX GN \
  --outSAMunmapped Within \
  --soloType CB_UMI_Complex \
  --soloMultiMappers EM \
  --soloFeatures "${features}" \
  --soloUMIdedup Exact \
  --soloCBposition "${CB_POS[@]}" \
  --soloUMIposition "${UMI_POS}" \
  --soloBarcodeReadLength 1 \
  --soloCBmatchWLtype "${cb_match}"
