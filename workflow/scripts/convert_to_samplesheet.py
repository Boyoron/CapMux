#!/usr/bin/env python3
import sys
import pandas as pd
import numpy as np

# pylint: disable=undefined-variable

sys.stderr = open(snakemake.log[0], "w")

def _make_header_df(demux_index: str) -> pd.DataFrame:
    """
    Create the fixed Illumina sample sheet header + [Data] column header row.
    Uses 10 columns (0..9), 19 rows.
    """

    demux_index = str(demux_index).strip()

    if demux_index == "i7":
        data_headers = [
            "Sample_ID", "Sample_Name", "Sample_Plate", "Sample_Well",
            "I7_Index_ID", "index", "Sample_Project", "Description", "", ""
        ]
    elif demux_index == "i7_i5":
        # index2 is inserted right after index
        data_headers = [
            "Sample_ID", "Sample_Name", "Sample_Plate", "Sample_Well",
            "I7_Index_ID", "index", "index2", "Sample_Project", "Description", ""
        ]
    elif demux_index == "i5":
        data_headers = [
            "Sample_ID", "Sample_Name", "Sample_Plate", "Sample_Well",
            "I7_Index_ID", "index", "Sample_Project", "Description", "", ""
        ]
    else:
        sys.exit(
            f"error: invalid demux_index='{demux_index}'. "
            "Allowed: i7, i5, i7_i5."
        )

    col0 = [
        "[Header]", "IEMFileVersion", "Investigator Name", "Experiment Name",
        "Date", "Workflow", "Application", "Assay", "Description", "Chemistry",
        "", "[Reads]", "", "", "", "[Settings]", "", "[Data]", data_headers[0]
    ]

    header_dict = {"0": col0}

    # columns 1..9: blank for first 18 rows, then the data header name
    for i in range(1, 10):
        header_dict[str(i)] = ([""] * 18) + [data_headers[i]]

    return pd.DataFrame(header_dict)


def convert_to_samplesheet(inp: str, out: str, demux_index: str) -> None:
    demux_index = str(demux_index).strip()

    ex_sheet = pd.read_excel(io=inp)

    required_cols = {"project_id", "sample_id"}
    if demux_index == "i7":
        required_cols |= {"index_seq"}
    elif demux_index == "i7_i5":
        required_cols |= {"index_seq", "index2_seq"}
    elif demux_index == "i5":
        required_cols |= {"index2_seq"}
    else:
        sys.exit(
            f"error: invalid demux_index='{demux_index}'. "
            "Allowed: i7, i5, i7_i5."
        )

    missing = sorted(required_cols - set(ex_sheet.columns))
    if missing:
        sys.exit(
            "error: extended_sample_sheet.xlsx is missing required column(s): "
            + ", ".join(missing)
        )

    # validate required fields are not empty
    if ex_sheet.loc[:, list(required_cols)].isnull().any().any():
        sys.exit(
            "error: check if there are no empty entries in required columns: "
            + ", ".join(sorted(required_cols))
        )

    header_df = _make_header_df(demux_index)

    # keep everything as strings
    sample_id = ex_sheet["sample_id"].astype("string").fillna("")
    project_id = ex_sheet["project_id"].astype("string").fillna("")

    ex_sheet_dict = {"0": sample_id}

    if demux_index == "i7":
        ex_sheet_dict["5"] = ex_sheet["index_seq"].astype("string").fillna("")
        ex_sheet_dict["6"] = project_id

    elif demux_index == "i7_i5":
        ex_sheet_dict["5"] = ex_sheet["index_seq"].astype("string").fillna("")
        ex_sheet_dict["6"] = ex_sheet["index2_seq"].astype("string").fillna("")
        ex_sheet_dict["7"] = project_id

    elif demux_index == "i5":
        # use index2_seq from .xlsx, but keep column name as "index"
        ex_sheet_dict["5"] = ex_sheet["index2_seq"].astype("string").fillna("")
        ex_sheet_dict["6"] = project_id

    sheet_df = pd.DataFrame(ex_sheet_dict)

    sheet = pd.concat([header_df, sheet_df], ignore_index=True).fillna("")
    sheet.to_csv(out, header=None, index=None)


convert_to_samplesheet(
    inp=snakemake.input[0],
    out=snakemake.output[0],
    demux_index=snakemake.params["demux_index"],
)
