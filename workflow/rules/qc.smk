rule fastqc:
    input:
        os.path.join(
            config["out_dir"],
            "merged_fastq/{sample}_{read}.fastq.gz",
        ),
    output:
        html=os.path.join(config["out_dir"], "qc/fastqc/{sample}_{read}.html"),
        zip=os.path.join(config["out_dir"], "qc/fastqc/{sample}_{read}_fastqc.zip"),
    log:
        os.path.join(config["out_dir"], "logs/fastqc_{sample}_{read}.log"),
    threads: 2
    resources:
        mem_mb=8000,
    wrapper:
        "v3.10.0/bio/fastqc"


rule multiqc:
    input:
        demux=os.path.join(config["out_dir"], "demuxed/Stats"),
        fastqc=expand(
            os.path.join(
                config["out_dir"], "qc/fastqc/{sample}_{read}_fastqc.zip"
            ),
            sample=get_sample_ids,
            read=["cdna", "bc"],
        ),

        cutadapt=expand(
            os.path.join(
                config["out_dir"], "trimmed_fastq/{sample}_cutadapt_report.txt"),
            sample=get_sample_ids,
        ),

        star_logs=expand(
            os.path.join(config["out_dir"], "mapped/{sample}/{sample}_Log.final.out"),
            sample=get_sample_ids,
        ),
        star=expand(
            os.path.join(config["out_dir"], "mapped/{sample}/{sample}_Solo.out"),
            sample=get_sample_ids,
        ),
        star_stats=expand(
            os.path.join(
                config["out_dir"], "mapped/{sample}/{sample}_Solo.out/{sample}_{file}"
            ),
            sample=get_sample_ids,
            file=[
                "Summary.csv",
                "UMIperCellSorted.txt",
                "Features.stats",
                "Barcodes.stats",
            ],
        ),
        config_file="config/multiqc_config.yaml",
    output:
        html=os.path.join(config["out_dir"], "qc/multiqc_report.html"),
    threads: 1
    params:
        extra="--verbose",
    log:
        os.path.join(config["out_dir"], "logs/multiqc.log"),
    conda:
        "../envs/multiqc.yaml"
    shell:
        """
        set -euo pipefail

        multiqc -f \
            {input.demux} \
            {input.fastqc} \
            {input.cutadapt} \
            {input.star_logs} \
            {input.star} \
            -c {input.config_file} \
            --outdir $(dirname {output.html}) \
            --filename $(basename {output.html}) \
            {params.extra} &> {log}
        """
