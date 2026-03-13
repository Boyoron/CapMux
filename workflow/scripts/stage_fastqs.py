"""Stage user-provided FASTQ files for the pipeline.

This script is used by rule `stage_fastqs` in `demux.smk` when config["run_mode"] == "fastq".

It validates Illumina style FASTQ names and symlinks them into the staged directory i.e.: {out_dir}/demuxed/<basename(fastq_dir)>/

It also ensures directory {out_dir}/demuxed/Stats/ exists.
"""

import os
import re
import glob

def _die(msg: str) -> None:
    raise ValueError(msg)

fastq_dir = str(snakemake.input.fastq_dir).strip()

if not fastq_dir:
    _die("run_mode=fastq but config['fastq_dir'] is empty")

if not os.path.isdir(fastq_dir):
    _die(f"fastq_dir does not exist or is not a directory: {fastq_dir!r}")

dest_root = str(snakemake.output.stage_dir)
stats_dir = str(snakemake.output.stats_dir)
marker_path = str(snakemake.output.marker)

os.makedirs(dest_root, exist_ok=True)
os.makedirs(stats_dir, exist_ok=True)

# searching for fastq.gz files with Illumina naming
name_re = re.compile(r"^.+?_S\d+(?:_L\d{3})?_(?:R1|R2|I1)_001\.fastq\.gz$")

candidates = sorted(glob.glob(os.path.join(fastq_dir, "*.fastq.gz")))
if not candidates:
    _die(f"No *.fastq.gz files found in fastq_dir={fastq_dir!r}")

bad = [os.path.basename(p) for p in candidates if not name_re.match(os.path.basename(p))]
if bad:
    _die(
        "Found FASTQs that do not match expected Illumina naming. "
        "Please rename them first. Examples: SAMPLE_S1_R1_001.fastq.gz or SAMPLE_S1_L001_R1_001.fastq.gz. "
        f"Offenders (first 10): {bad[:10]}"
    )

# symlink into results
for src in candidates:
    dst = os.path.join(dest_root, os.path.basename(src))
    if os.path.lexists(dst):
        os.remove(dst)
    os.symlink(os.path.abspath(src), dst)

# writing a marker after successful staging
with open(marker_path, "w") as fh:
    fh.write(f"source_fastq_dir\t{os.path.abspath(fastq_dir)}\n")
    fh.write(f"staged_dir\t{os.path.abspath(dest_root)}\n")
    fh.write(f"n_fastqs\t{len(candidates)}\n")
    fh.write("files\t" + ",".join(os.path.basename(p) for p in candidates) + "\n")

# writing log
log_path = snakemake.log[0] if getattr(snakemake, "log", None) else None
if log_path:
    with open(log_path, "w") as fh:
        fh.write(f"Staged {len(candidates)} FASTQ(s) from: {fastq_dir}\n")
        fh.write(f"Into: {dest_root}\n")
