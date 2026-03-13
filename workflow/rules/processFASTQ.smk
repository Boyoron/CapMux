import os

def _seg(name, field, default):
    return config.get("barcode", {}).get("segments", {}).get(name, {}).get(field, default)

if config.get("demux_by") == "index":
    rule merge:
        """
        Reorder barcode segments in R1 to build cell barcode FASTQs (bc), copy R2 as cDNA (cdna).
        Output files go to merged_fastq/.
        """
        input:
            r1 = lambda w: get_fastqs_for_sample(w)["r1"],
            r2 = lambda w: get_fastqs_for_sample(w)["r2"],
        output:
            bc = os.path.join(config["out_dir"], "merged_fastq", "{sample}_bc.fastq.gz"),
            cdna = os.path.join(config["out_dir"], "merged_fastq", "{sample}_cdna.fastq.gz"),
        threads: 1
        log:
            os.path.join(config["out_dir"], "logs", "merge_{sample}.log"),
        conda:
            "../envs/coreutils.yaml"
        params:
            bc1_start  = _seg("bc1", "start", 5),
            bc1_len    = _seg("bc1", "len", 8),
            bc2_start  = _seg("bc2", "start", 18),
            bc2_len    = _seg("bc2", "len", 10),
            bc3_start  = _seg("bc3", "start", 32),
            bc3_len    = _seg("bc3", "len", 8),
            umi1_start = _seg("umi1", "start", 1),
            umi1_len   = _seg("umi1", "len", 4),
            umi2_start = _seg("umi2", "start", 40),
            umi2_len   = _seg("umi2", "len", 4),
        script:
            "../scripts/merge.sh"


rule trim_reads:
    """
    Trim cDNA FASTQs (cdna) with cutadapt and keep trimmed files in trimmed_fastq/.
    Uses merged_fastq outputs from the 'merge' rule as inputs.
    """
    input:
        bc = os.path.join(
            config["out_dir"], "merged_fastq", "{sample}_bc.fastq.gz"
        ),
        cdna = os.path.join(
            config["out_dir"], "merged_fastq", "{sample}_cdna.fastq.gz"
        ),
    output:
        bc_trimmed = os.path.join(
            config["out_dir"], "trimmed_fastq", "{sample}_bc.fastq.gz"
        ),
        cdna_trimmed = os.path.join(
            config["out_dir"], "trimmed_fastq", "{sample}_cdna.fastq.gz"
        ),
        report = os.path.join(
            config["out_dir"], "trimmed_fastq", "{sample}_cutadapt_report.txt"
        ),
    threads: 2
    log:
        os.path.join(config["out_dir"], "logs", "trim_fastq_{sample}.log")
    conda:
        "../envs/trim.yaml"
    params:
        tso_seq = config.get("tso_seq")
    script:
        "../scripts/trim.sh"
