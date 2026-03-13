#!/usr/bin/env python3
import os
import glob
import re
import pandas as pd
from io import StringIO

# pylint: disable=undefined-variable

sample = snakemake.params.sample
demux_dir = snakemake.input.solo_dir
output_files = snakemake.output

def modify_summary(content):
    df = pd.read_csv(StringIO(content), header=None)
    df = df.transpose()
    new_header = ['Sample Name'] + df.iloc[0, 1:].tolist()
    new_values = [sample] + df.iloc[1, 1:].tolist()
    df.columns = new_header
    df = pd.DataFrame([new_values], columns=new_header)
    return df.to_csv(index=False)


def modify_umi(content):
    lines = content.strip().splitlines()
    numbered_lines = [f"{i+1}\t{line}" for i, line in enumerate(lines)]
    return "\n".join(numbered_lines) + "\n"


def modify_stats(content):
    lines = content.strip().splitlines()
    processed = []
    for line in lines:
        line = re.sub(r'^\s+', '', line)
        line = re.sub(r'\s+', '\t', line)
        processed.append(line)
    return "\n".join(processed) + "\n"


file_map = {
    "Summary.csv": (modify_summary, output_files.summary),
    "UMIperCellSorted.txt": (modify_umi, output_files.umi),
    "Barcodes.stats": (modify_stats, output_files.barcodes),
    "Features.stats": (modify_stats, output_files.features),
}


for suffix, (modifier, outpath) in file_map.items():
    pattern = os.path.join(demux_dir, f"**/*{suffix}")
    matches = glob.glob(pattern, recursive=True)
    if not matches:
        raise FileNotFoundError(f"No file matching {pattern}")
    if len(matches) > 1:
        raise RuntimeError(f"Multiple matches for {pattern}: {matches}")
    
    inpath = matches[0]
    with open(inpath, "r") as infile:
        content = infile.read()

    modified = modifier(content)

    with open(outpath, "w") as outfile:
        outfile.write(modified)
