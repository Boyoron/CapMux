#!/usr/bin/env python3
import glob
import os
import re

merged_dir = snakemake.input.merged_dir
output_file = snakemake.output.sample_ids

pattern = os.path.join(merged_dir, "*_bc.fastq.gz")

samples = []
for path in glob.glob(pattern):
    base = os.path.basename(path)
    if base.startswith("Undetermined"):
        continue
    sample = re.sub(r"_bc\.fastq\.gz$", "", base)
    samples.append(sample)

samples = sorted(set(samples))

os.makedirs(os.path.dirname(output_file), exist_ok=True)
with open(output_file, "w") as out:
    out.write("\n".join(samples) + "\n")