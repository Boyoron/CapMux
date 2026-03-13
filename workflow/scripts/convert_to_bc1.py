#!/usr/bin/env python3
import sys
import argparse
import pandas as pd
import numpy as np

# pylint: disable=undefined-variable

sys.stderr = open(snakemake.log[0], "w")

bc1_list_path = snakemake.params.bc1_list_path

def convert_to_bc1(inp, out):
    
    xls = pd.ExcelFile(inp)

    if "bc1" in xls.sheet_names:
        bc1_list = pd.read_csv(bc1_list_path, names=["bc1_seq"])

        extended_sample_sheet = (
            pd.read_excel(io=inp, sheet_name="bc1", header=None)
              .iloc[1:, 1:]
              .astype("string")
              .fillna("UNASSIGNED")
        )

        sample_list = (
            extended_sample_sheet
              .stack()
              .reset_index(drop=True)
              .to_frame(name="sample_id")
        )

        # ---- if n_bc1 != n_samples ----
        n_bc1 = len(bc1_list)
        n_samples = len(sample_list)
        if n_bc1 != n_samples:
            sys.exit(
                "ERROR: the count of barcodes in bc1_list does not match bc1 sheet sample count "
                f"(bc1_list={n_bc1}, samples_in_sheet={n_samples}). "
                "Check your bc1_list.txt and the 'bc1' sheet dimensions."
            )
        # -------------------------------

        bc1_to_sample = pd.concat([bc1_list, sample_list], axis=1)
        bc1_to_sample.to_csv(out, index=False)

    # ---- if mandatory bc1 sheet is missing ----
    else:
        sys.exit("ERROR: check if there is bc1 sheet in the extended_sample_sheet.xlsx")
    # -------------------------------------------

convert_to_bc1(
    inp=snakemake.input[0],
    out=snakemake.output[0]
)
