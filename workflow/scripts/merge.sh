#!/usr/bin/env bash
set -euo pipefail

exec > "${snakemake_log[0]}" 2>&1

# inputs
R1="${snakemake_input[r1]}"
R2="${snakemake_input[r2]}"

# outputs
BC="${snakemake_output[bc]}"
CDNA="${snakemake_output[cdna]}"

# segment coords
bc1_start="${snakemake_params[bc1_start]}"; bc1_len="${snakemake_params[bc1_len]}"
bc2_start="${snakemake_params[bc2_start]}"; bc2_len="${snakemake_params[bc2_len]}"
bc3_start="${snakemake_params[bc3_start]}"; bc3_len="${snakemake_params[bc3_len]}"
umi1_start="${snakemake_params[umi1_start]}"; umi1_len="${snakemake_params[umi1_len]}"
umi2_start="${snakemake_params[umi2_start]}"; umi2_len="${snakemake_params[umi2_len]}"

is_uint() { [[ "$1" =~ ^[0-9]+$ ]]; }

# return 0 if enabled, 1 if disabled.
validate_segment() {
  local name="$1" s="$2" l="$3"

  if ! is_uint "$s" || ! is_uint "$l"; then
    echo "ERROR: segment '$name' has non-integer coords: start='$s' len='$l'" >&2
    exit 2
  fi

  if [[ "$s" -eq 0 && "$l" -eq 0 ]]; then
    echo "INFO: segment '$name' disabled (start=0,len=0) -> skipping" >&2
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

# enforce at least one barcode and at least one UMI enabled
USED_BC=0
USED_UMI=0
if validate_segment bc1 "$bc1_start" "$bc1_len"; then USED_BC=$((USED_BC+1)); fi
if validate_segment bc2 "$bc2_start" "$bc2_len"; then USED_BC=$((USED_BC+1)); fi
if validate_segment bc3 "$bc3_start" "$bc3_len"; then USED_BC=$((USED_BC+1)); fi
if validate_segment umi1 "$umi1_start" "$umi1_len"; then USED_UMI=$((USED_UMI+1)); fi
if validate_segment umi2 "$umi2_start" "$umi2_len"; then USED_UMI=$((USED_UMI+1)); fi

if [[ "$USED_BC" -lt 1 ]]; then
  echo "ERROR: At least one of bc1/bc2/bc3 must be enabled." >&2
  exit 2
fi
if [[ "$USED_UMI" -lt 1 ]]; then
  echo "ERROR: At least one of umi1/umi2 must be enabled." >&2
  exit 2
fi

# building awk concatenation expression in fixed order: bc1, bc2, bc3, umi1, umi2
expr=""
STRUCT=""

add_expr() {
  local name="$1" s="$2" l="$3"
  expr+="substr(\$0,${s},${l})"
  if [[ -z "$STRUCT" ]]; then
    STRUCT="${name}(${s},${l})"
  else
    STRUCT="${STRUCT} + ${name}(${s},${l})"
  fi
}

if validate_segment bc1 "$bc1_start" "$bc1_len"; then add_expr "bc1" "$bc1_start" "$bc1_len"; fi
if validate_segment bc2 "$bc2_start" "$bc2_len"; then add_expr "bc2" "$bc2_start" "$bc2_len"; fi
if validate_segment bc3 "$bc3_start" "$bc3_len"; then add_expr "bc3" "$bc3_start" "$bc3_len"; fi
if validate_segment umi1 "$umi1_start" "$umi1_len"; then add_expr "umi1" "$umi1_start" "$umi1_len"; fi
if validate_segment umi2 "$umi2_start" "$umi2_len"; then add_expr "umi2" "$umi2_start" "$umi2_len"; fi

echo "INFO: Concatenating segments (R1 only): ${STRUCT}"
echo "INFO: AWK expr: ${expr}"

zcat "${R1}" \
  | awk "NR%4==2 || NR%4==0 { \$0=${expr} } 1" \
  | gzip > "${BC}"

cp "${R2}" "${CDNA}"
