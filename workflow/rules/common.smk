import os
import re
import glob
import pandas as pd
from snakemake.utils import validate


def get_sample_ids(wildcards):
    if config.get("demux_by") == "index": 
        ckpt = checkpoints.parse_demux.get()
        sample_file = ckpt.output.sample_ids
        with open(sample_file) as f:
            samples = f.read().strip().splitlines()
        return samples
    if config.get("demux_by") == "bc1": 
        """
        Sample IDs after BC1 demux.
        Uses checkpoint parse_bc1_demux which:
        """
        ckpt = checkpoints.parse_bc1_demux.get()
        sample_file = ckpt.output.sample_ids
        with open(sample_file) as f:
            samples = [line.strip() for line in f if line.strip()]
        return samples


def get_demux_root(wildcards):
    ckpt = checkpoints.parse_demux.get()
    with open(ckpt.output.demux_root) as f:
        demux_root = f.read().strip()
        return demux_root


def get_fastqs_for_sample(wildcards):
    """
    Detect the demultiplexed FASTQ files for a given sample and return R1, R2 and I1.
    """
    pattern = os.path.join(
        config["out_dir"], "demuxed", "**", f"{wildcards.sample}_S*_*_001.fastq.gz"
    )
    matches = sorted(glob.glob(pattern, recursive=True))

    if len(matches) < 3:
        raise FileNotFoundError(
            f"Did not find at least 3 FASTQ files (R1, R2, I1) for id={wildcards.sample} matching {pattern}"
        )

    R1 = next((f for f in matches if "_R1_" in f), None)
    R2 = next((f for f in matches if "_R2_" in f), None)
    I1 = next((f for f in matches if "_I1_" in f), None)
    
    if not (R1 and R2 and I1):
        raise FileNotFoundError(
            f"Could not identify R1, R2 and I1 among these matches for id={wildcards.sample}: {matches}"
        )

    return {"r1": R1, "r2": R2, "i1": I1}


if config.get("demux_by") == "bc1":

    def get_bc1_chunks(wildcards):
        """
        Return a list of 'chunk' prefixes for each R1/R2/I1 triplet in the demuxed run.

        Example:
            demuxed/<demux_root>/XYZ_S1_L001_R1_001.fastq.gz
        ->  chunk = "XYZ_S1_L001"
        """
        demux_root = get_demux_root(wildcards)
        demux_dir = os.path.join(config["out_dir"], "demuxed", demux_root)

        pattern = os.path.join(demux_dir, "*_R1_*.fastq.gz")
        r1_files = glob.glob(pattern)

        chunks = []
        for path in r1_files:
            base = os.path.basename(path)
            
            chunk = re.sub(r"_R1_001\.fastq\.gz$", "", base)
            chunks.append(chunk)

        if not chunks:
            raise ValueError(f"No R1 FASTQs found matching {pattern}")

        return sorted(set(chunks))

    def bc1_done_files(wildcards):
        """
        Expand to all .done flags produced by the per-chunk demux rule.
        Used as an input to the merge_bc1_demux rule.
        """
        chunks = get_bc1_chunks(wildcards)
        return expand(
            os.path.join(config["out_dir"], "bc1_demux_tmp", "{chunk}.done"),
            chunk=chunks,
        )
