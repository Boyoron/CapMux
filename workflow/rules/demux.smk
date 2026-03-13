import os
import glob

RUN_MODE = config.get('run_mode', 'bcl')

if RUN_MODE == "fastq":
    FASTQ_STAGE_NAME = os.path.basename(str(config["fastq_dir"]).rstrip("/"))

INDEX_WL_AUTO = os.path.join("assets", "barcodes", "index_list.txt")
INDEX_WL_USER = config.get("barcode", {}).get("whitelists", {}).get("index", "")

if str(INDEX_WL_USER).strip() == "":
    rule generate_index_whitelist:
        input:
            config["sample_sheet"],
        output:
            INDEX_WL_AUTO,
        threads: 1
        log:
            os.path.join(config["out_dir"], "logs/generate_index_whitelist.log"),
        conda:
            "../envs/pandas.yaml"
        script:
            "../scripts/generate_index_whitelist.py"


rule convert_sheet:
    """
    Converts extended_sample_sheet_template.xlsx to standard sample_sheet.csv.
    """
    input:
        inp=config["sample_sheet"],
    output:
        out=os.path.join(config["out_dir"], "sample_sheet.csv"),
    threads: 1
    log:
        os.path.join(config["out_dir"], "logs/convert_sheet.log"),
    conda:
        "../envs/pandas.yaml"
    params:
        demux_index = config["demux_index"]
    script:
        "../scripts/convert_to_samplesheet.py"


if config.get("demux_by") == "bc1": 
    rule convert_bc1:
        """
        Converts extended_sample_sheet_template.xlsx sheet "bc1" to bc1_to_sample.csv.
        """
        input:
            inp=config["sample_sheet"],
        output:
            out=os.path.join(config["out_dir"], "bc1_to_sample.csv"),
        threads: 1
        log:
            os.path.join(config["out_dir"], "logs/convert_bc1.log"),
        conda:
            "../envs/pandas.yaml"
        params:
            bc1_list_path = config.get("barcode", {}).get("whitelists", {}).get("bc1", "")
        script:
            "../scripts/convert_to_bc1.py"


if RUN_MODE == "fastq":
    rule stage_fastqs:
        """
        Stage user-provided FASTQ files into {out_dir}/demuxed/<basename(fastq_dir)>/.

        Requirements:
          - config["fastq_dir"] points to a folder containing Illumina-style FASTQ named files, e.g.
              SAMPLE_S1_R1_001.fastq.gz
              SAMPLE_S1_R2_001.fastq.gz
              SAMPLE_S1_I1_001.fastq.gz
        """
        input:
            fastq_dir=config["fastq_dir"],
        output:
            demux_dir=directory(os.path.join(config["out_dir"], "demuxed")),
            stage_dir=directory(os.path.join(config["out_dir"], "demuxed", FASTQ_STAGE_NAME)),
            stats_dir=directory(os.path.join(config["out_dir"], "demuxed", "Stats")),
            marker=os.path.join(config["out_dir"], "demuxed", FASTQ_STAGE_NAME, ".staged.ok"),
        threads: 1
        log:
            os.path.join(config["out_dir"], "logs/stage_fastqs.log"),
        conda:
            "../envs/coreutils.yaml"
        script:
            "../scripts/stage_fastqs.py"

else:
    rule demux:
        """
        Runs bcl2fastq to demultiplex Illumina BCL files.
        """
        input:
            run_dir=config["run_dir"],
            sample_sheet=os.path.join(config["out_dir"], "sample_sheet.csv"),
        output:
            out_dir=directory(os.path.join(config["out_dir"], "demuxed")),
            sta_dir=directory(os.path.join(config["out_dir"], "demuxed/Stats")),
        threads: config.get("threads", 8)
        log:
            os.path.join(config["out_dir"], "logs/demux.log"),
        conda:
            "../envs/bcl2fastq.yaml"
        params:
            use_bases_mask = config["use_bases_mask"],
            demux_index = config["demux_index"]
        script:
            "../scripts/bcl2fastq.sh"


checkpoint parse_demux:
    input:
        demux_dir = os.path.join(config["out_dir"], "demuxed"),
        staged_ok = os.path.join(config["out_dir"], "demuxed", FASTQ_STAGE_NAME, ".staged.ok") if RUN_MODE == "fastq" else [],
    output:
        sample_ids = os.path.join(config["out_dir"], "sample_ids.txt"),
        demux_root = os.path.join(config["out_dir"], "demux_root.txt"),
    threads: 1
    log:
        os.path.join(config["out_dir"], "logs/parse_demux.log"),
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/parse_demux.py"


if config.get("demux_by") == "bc1": 
    rule bc1_demux_chunk:
        """
        Run capseq_demux.py on R1/R2/I1 FASTQ files.
        Outputs FASTQs into bc1_demux_tmp/ and a per-chunk .done flag.
        """
        input:
            r1 = lambda w: os.path.join(
                config["out_dir"], "demuxed", get_demux_root(w), f"{w.chunk}_R1_001.fastq.gz"
            ),
            r2 = lambda w: os.path.join(
                config["out_dir"], "demuxed", get_demux_root(w), f"{w.chunk}_R2_001.fastq.gz"
            ),
            i1 = lambda w: os.path.join(
                config["out_dir"], "demuxed", get_demux_root(w), f"{w.chunk}_I1_001.fastq.gz"
            ),
            bc1_csv = os.path.join(config["out_dir"], "bc1_to_sample.csv"),
        output:
            done = touch(os.path.join(config["out_dir"], "bc1_demux_tmp", "{chunk}.done")),
        threads: 1
        log:
            os.path.join(config["out_dir"], "logs", "capseq_demux_{chunk}.log"),
        conda:
            "../envs/coreutils.yaml"
        params:
            outdir = os.path.join(config["out_dir"], "bc1_demux_tmp"),
            allow_mismatches = config.get("barcode", {}).get("allow_mismatches", 1),
            bc1_start  = config.get("barcode", {}).get("segments", {}).get("bc1",  {}).get("start", 32),
            bc1_len    = config.get("barcode", {}).get("segments", {}).get("bc1",  {}).get("len",   8),
            bc2_start  = config.get("barcode", {}).get("segments", {}).get("bc2",  {}).get("start", 18),
            bc2_len    = config.get("barcode", {}).get("segments", {}).get("bc2",  {}).get("len",   10),
            bc3_start  = config.get("barcode", {}).get("segments", {}).get("bc3",  {}).get("start", 5),
            bc3_len    = config.get("barcode", {}).get("segments", {}).get("bc3",  {}).get("len",   8),
            umi1_start = config.get("barcode", {}).get("segments", {}).get("umi1", {}).get("start", 1),
            umi1_len   = config.get("barcode", {}).get("segments", {}).get("umi1", {}).get("len",   4),
            umi2_start = config.get("barcode", {}).get("segments", {}).get("umi2", {}).get("start", 40),
            umi2_len   = config.get("barcode", {}).get("segments", {}).get("umi2", {}).get("len",   4),
        shell:
            r"""
            set -euo pipefail
            mkdir -p "{params.outdir}"

            python workflow/scripts/capseq_demux.py \
                --r1 {input.r1} \
                --r2 {input.r2} \
                --i1 {input.i1} \
                --output-dir "{params.outdir}" \
                --bc1-csv {input.bc1_csv} \
                --allow-mismatches {params.allow_mismatches} \
                --bc1-start {params.bc1_start} --bc1-len {params.bc1_len} \
                --bc2-start {params.bc2_start} --bc2-len {params.bc2_len} \
                --bc3-start {params.bc3_start} --bc3-len {params.bc3_len} \
                --umi1-start {params.umi1_start} --umi1-len {params.umi1_len} \
                --umi2-start {params.umi2_start} --umi2-len {params.umi2_len} \
                &> {log}

            touch {output.done}
            """


    rule merge_bc1_demux:
        """
        Merge per-chunk bc1 demux outputs from bc1_demux_tmp into final merged_fastq/ using merge_chunk.sh.
        bc1_demux_tmp folder is deleted after merge.
        """
        input:
            done = bc1_done_files,
        output:
            merged_dir = directory(os.path.join(config["out_dir"], "merged_fastq")),
        threads: config["threads"]
        log:
            os.path.join(config["out_dir"], "logs", "merge_bc1_demux.log"),
        conda:
            "../envs/coreutils.yaml"
        params:
            tmpdir = os.path.join(config["out_dir"], "bc1_demux_tmp"),
        shell:
            """
            set -euo pipefail

            bash workflow/scripts/merge_chunk.sh \
                "{params.tmpdir}" \
                "{output.merged_dir}" \
                {threads} \
                &> {log}
            """


if config.get("demux_by") == "bc1": 
    checkpoint parse_bc1_demux:
        """
        Inspect merged_fastq/ and write a sample list based on bc1 demux outputs, excluding Undetermined_* files.
        """
        input:
            merged_dir = os.path.join(config["out_dir"], "merged_fastq"),
        output:
            sample_ids = os.path.join(config["out_dir"], "sample_ids_bc1.txt"),
        threads: 1
        log:
            os.path.join(config["out_dir"], "logs/parse_bc1_demux.log"),
        conda:
            "../envs/pandas.yaml"
        script:
            "../scripts/parse_bc1_demux.py"