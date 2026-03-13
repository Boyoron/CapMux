#!/usr/bin/env python3
import os
import pandas as pd



sheet_path = str(snakemake.input[0])
out_path = str(snakemake.output[0])

df = pd.read_excel(sheet_path)

idx = df["index_seq"].dropna().astype(str).str.strip()
idx = idx[idx != ""]

seen = set()
unique = []
for s in idx.tolist():
    if s not in seen:
        seen.add(s)
        unique.append(s)

os.makedirs(os.path.dirname(out_path), exist_ok=True)
with open(out_path, "w") as f:
    for s in unique:
        f.write(s + "\n")
