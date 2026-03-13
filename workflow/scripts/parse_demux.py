#!/usr/bin/env python3
import os
import re

# pylint: disable=undefined-variable

demux_dir = snakemake.input.demux_dir
sample_ids_file = snakemake.output.sample_ids
demux_root_file = snakemake.output.demux_root

sample_name_set = set()
sample_dirs = set()

for root, dirs, files in os.walk(demux_dir):
    for f in files:
        m = re.match(r"^(.+?)_S\d+(?:_L\d{3})?_(?:R[123]|I1)_001\.fastq\.gz$", f)
        if m and m.group(1) != "Undetermined":
            sample_name_set.add(m.group(1))
            sample_dirs.add(root)

if not sample_name_set:
    raise ValueError(f"No demuxed FASTQ files found under {demux_dir!r}")

if len(sample_dirs) > 1:
    raise ValueError(
        "Found demuxed FASTQ files in multiple directories under "
        f"{demux_dir!r}: {sorted(sample_dirs)}; expected a single subfolder."
    )

# full path
full_dir = next(iter(sample_dirs))
# only folder name
folder_name = os.path.basename(full_dir.rstrip("/"))

# creating sample_ids.txt with sample ID
samples = sorted(sample_name_set)
with open(sample_ids_file, "w") as outF:
    for s in samples:
        outF.write(s + "\n")

# creating demux_root.txt with folder name
with open(demux_root_file, "w") as f:
    f.write(folder_name + "\n")
