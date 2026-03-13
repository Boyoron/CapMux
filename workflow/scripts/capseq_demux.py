#!/usr/bin/env python3
import csv
import gzip
import os
import argparse
from collections import defaultdict

BASES = "ACGT"
AMBIG = object()


def safe_slice(s, start, end):
    """Slice s[start:end], returning '' if indices out of range."""
    if s is None:
        return ""
    if start < 0 or end < 0:
        return ""
    if start >= len(s) or start >= end:
        return ""
    return s[start:end]


def open_maybe_gzip(path, mode="rt"):
    """
    Open a file that may be gzipped (.gz) or plain text.
    """
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    else:
        return open(path, mode)


def normalize_segment(name, start, length):
    """
    Rules:
      - start=0 and len=0 => disabled
      - start=0 xor len=0 => error
      - enabled => start>=1 and len>=1 (1-based start)
    Returns: (enabled: bool, start0: int, end0: int)
    """
    if start is None or length is None:
        raise ValueError(f"{name}: start/len must be provided")

    if start == 0 and length == 0:
        return False, 0, 0

    if (start == 0) != (length == 0):
        raise ValueError(f"{name}: must be disabled as start=0,len=0 OR enabled with start>=1,len>=1 (got start={start}, len={length})")

    if start < 1 or length < 1:
        raise ValueError(f"{name}: enabled segments require start>=1,len>=1 (got start={start}, len={length})")

    start0 = start - 1
    end0 = start0 + length
    return True, start0, end0


def load_bc1_map(csv_path, allow_mismatches, bc1_len):
    """
    Load bc1_seq -> sample_id from CSV with columns: bc1_seq,sample_id.
    Build 1-mismatch neighbors if allow_mismatches>=1.

    bc1_len is used to validate / generate neighbors.

    Returns:
      exact: dict of exact BC1 -> sample
      neighbors: dict of 1-mismatch BC1 -> sample or AMBIG
    """
    exact = {}
    neighbors = {}

    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            bc = row["bc1_seq"].strip()
            sample = row["sample_id"].strip()
            exact[bc] = sample

    if allow_mismatches >= 1:
        for bc, sample in exact.items():
            if len(bc) != bc1_len:
                continue
            for i in range(len(bc)):
                for b in BASES:
                    if b == bc[i]:
                        continue
                    nb = bc[:i] + b + bc[i + 1 :]
                    prev = neighbors.get(nb)
                    if prev is None:
                        neighbors[nb] = sample
                    elif prev is not sample:
                        neighbors[nb] = AMBIG

    return exact, neighbors


def assign_sample_from_bc1(seq, exact, neighbors, allow_mismatches):
    """
    seq: observed BC1 sequence

    Returns:
      (sample_id or None, status_code)

    status_code:
        0 = exact BC1 match
        1 = 1-mismatch BC1 match (unique)
        2 = ambiguous 1-mismatch
        3 = no BC1 match
    """
    # exact match
    if seq in exact:
        return exact[seq], 0

    if allow_mismatches < 1:
        return None, 3
        
    # 1-mm neighbor map
    candidate = neighbors.get(seq)
    if candidate is None:
        return None, 3
    if candidate is AMBIG:
        return None, 2
    return candidate, 1


def extract_segments(seq, qual, seg_defs):
    """
    seg_defs: dict name -> (enabled, start0, end0)
    Returns: dict name -> (seg_seq, seg_qual) for enabled segments,
             disabled segments map to ("","").
    """
    out = {}
    for name, (enabled, start0, end0) in seg_defs.items():
        if not enabled:
            out[name] = ("", "")
        else:
            out[name] = (safe_slice(seq, start0, end0), safe_slice(qual, start0, end0))
    return out


def demux_all(
    r1_path,
    r2_path,
    i1_path,
    output_dir,
    bc1_csv,
    allow_mismatches,
    bc1_start, bc1_len,
    bc2_start, bc2_len,
    bc3_start, bc3_len,
    umi1_start, umi1_len,
    umi2_start, umi2_len,
):
    os.makedirs(output_dir, exist_ok=True)

    # normalize segment definitions
    seg_defs = {}
    seg_defs["bc1"]  = normalize_segment("bc1",  bc1_start, bc1_len)
    seg_defs["bc2"]  = normalize_segment("bc2",  bc2_start, bc2_len)
    seg_defs["bc3"]  = normalize_segment("bc3",  bc3_start, bc3_len)
    seg_defs["umi1"] = normalize_segment("umi1", umi1_start, umi1_len)
    seg_defs["umi2"] = normalize_segment("umi2", umi2_start, umi2_len)

    used_bc = sum(1 for k in ("bc1", "bc2", "bc3") if seg_defs[k][0])
    used_umi = sum(1 for k in ("umi1", "umi2") if seg_defs[k][0])

    if used_bc < 1:
        raise ValueError("Config error: at least one of bc1, bc2, bc3 must be enabled (not start=0,len=0).")
    if used_umi < 1:
        raise ValueError("Config error: at least one of umi1, umi2 must be enabled (not start=0,len=0).")

    if not seg_defs["bc1"][0]:
        raise ValueError("Config error: bc1 is disabled but this demux step requires bc1 (demultiplexing is by BC1).")

    bc1_len_effective = bc1_len
    exact, neighbors = load_bc1_map(bc1_csv, allow_mismatches, bc1_len_effective)

    common_prefix = os.path.commonprefix([
        os.path.basename(r1_path),
        os.path.basename(r2_path),
        os.path.basename(i1_path),
    ])

    out_R2 = {}
    out_BC = {}

    unknown_R2_path = os.path.join(output_dir, f"{common_prefix}Undetermined_cdna.fastq")
    unknown_BC_path = os.path.join(output_dir, f"{common_prefix}Undetermined_bc.fastq")

    unknown_R2 = open(unknown_R2_path, "w")
    unknown_BC = open(unknown_BC_path, "w")

    counts = defaultdict(int)

    for p in (r1_path, r2_path, i1_path):
        if not os.path.exists(p):
            raise FileNotFoundError(f"Input file not found: {p}")

    print(f"Processing R1 = {r1_path}")
    print(f"          R2 = {r2_path}")
    print(f"          I1 = {i1_path}")

    with open_maybe_gzip(r1_path, "rt") as r1, \
         open_maybe_gzip(r2_path, "rt") as r2, \
         open_maybe_gzip(i1_path, "rt") as i1:

        while True:
            h1 = r1.readline()
            if not h1:
                break

            s1 = r1.readline().rstrip()
            _p1 = r1.readline()
            q1 = r1.readline().rstrip()

            h2 = r2.readline()
            s2 = r2.readline().rstrip()
            p2 = r2.readline()
            q2 = r2.readline().rstrip()

            hi1 = i1.readline()
            if not hi1:
                break
            si1 = i1.readline().rstrip()
            pi1 = i1.readline()
            qi1 = i1.readline().rstrip()

            segs = extract_segments(s1, q1, seg_defs)

            bc1_seq, bc1_qual = segs["bc1"]

            # use BC1 for demuxing
            if len(bc1_seq) != bc1_len_effective:
                sample = None
                status = 3
            else:
                sample, status = assign_sample_from_bc1(bc1_seq, exact, neighbors, allow_mismatches)

            counts[status] += 1

            if sample is None:
                oR2 = unknown_R2
                oBC = unknown_BC
            else:
                if sample not in out_R2:
                    r2_out_path = os.path.join(output_dir, f"{common_prefix}{sample}_cdna.fastq")
                    bc_out_path = os.path.join(output_dir, f"{common_prefix}{sample}_bc.fastq")
                    out_R2[sample] = open(r2_out_path, "w")
                    out_BC[sample] = open(bc_out_path, "w")
                oR2 = out_R2[sample]
                oBC = out_BC[sample]

            # write R2 unchanged
            oR2.write(h2)
            oR2.write(s2 + "\n")
            oR2.write(p2)
            oR2.write(q2 + "\n")

            # index from I1
            index_seq = si1
            index_qual = qi1

            # build BC
            bc_parts = []
            q_parts = []

            for name in ("bc1", "bc2", "bc3"):
                if seg_defs[name][0]:
                    bc_parts.append(segs[name][0])
                    q_parts.append(segs[name][1])

            bc_parts.append(index_seq)
            q_parts.append(index_qual)

            for name in ("umi1", "umi2"):
                if seg_defs[name][0]:
                    bc_parts.append(segs[name][0])
                    q_parts.append(segs[name][1])

            bc_full_seq = "".join(bc_parts)
            bc_full_qual = "".join(q_parts)

            oBC.write(h1)
            oBC.write(bc_full_seq + "\n")
            oBC.write("+\n")
            oBC.write(bc_full_qual + "\n")

    for f in list(out_R2.values()) + list(out_BC.values()):
        f.close()
    unknown_R2.close()
    unknown_BC.close()

    print("Demux finished.")
    print("Status counts:")
    print("  0 = exact BC1 match       :", counts[0])
    print("  1 = 1-mismatch BC1 match  :", counts[1])
    print("  2 = ambiguous 1-mm        :", counts[2])
    print("  3 = no BC1 match          :", counts[3])


def parse_args():
    ap = argparse.ArgumentParser(
        description="Demultiplex by BC1 and build BC FASTQs from a single R1/R2/I1 triplet."
    )
    ap.add_argument("--r1", required=True)
    ap.add_argument("--r2", required=True)
    ap.add_argument("--i1", required=True)
    ap.add_argument("--output-dir", "-o", required=True)
    ap.add_argument("--bc1-csv", default="bc1_to_sample.csv")
    ap.add_argument("--allow-mismatches", type=int, default=1)

    ap.add_argument("--bc1-start", type=int, default=5)
    ap.add_argument("--bc1-len",   type=int, default=8)
    ap.add_argument("--bc2-start", type=int, default=18)
    ap.add_argument("--bc2-len",   type=int, default=10)
    ap.add_argument("--bc3-start", type=int, default=32)
    ap.add_argument("--bc3-len",   type=int, default=8)
    ap.add_argument("--umi1-start", type=int, default=1)
    ap.add_argument("--umi1-len",   type=int, default=4)
    ap.add_argument("--umi2-start", type=int, default=40)
    ap.add_argument("--umi2-len",   type=int, default=4)

    return ap.parse_args()


if __name__ == "__main__":
    args = parse_args()
    demux_all(
        r1_path=args.r1,
        r2_path=args.r2,
        i1_path=args.i1,
        output_dir=args.output_dir,
        bc1_csv=args.bc1_csv,
        allow_mismatches=args.allow_mismatches,
        bc1_start=args.bc1_start, bc1_len=args.bc1_len,
        bc2_start=args.bc2_start, bc2_len=args.bc2_len,
        bc3_start=args.bc3_start, bc3_len=args.bc3_len,
        umi1_start=args.umi1_start, umi1_len=args.umi1_len,
        umi2_start=args.umi2_start, umi2_len=args.umi2_len,
    )
