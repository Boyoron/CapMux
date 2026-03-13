import os

def _seg(name, field, default):
    return config.get("barcode", {}).get("segments", {}).get(name, {}).get(field, default)

def _wl(name, default=""):
    return config.get("barcode", {}).get("whitelists", {}).get(name, default)

INDEX_WL_AUTO = os.path.join("assets", "barcodes", "index_list.txt")

def _wl_index_effective():
    p = _wl("index", "")
    return INDEX_WL_AUTO if str(p).strip() == "" else p


rule starsoloUMI:
    """Aligns reads with STARsolo (dynamic CB/UMI layout from config.yaml)."""
    input:
        cdna = os.path.join(config["out_dir"], "trimmed_fastq", "{sample}_cdna.fastq.gz"),
        bc   = os.path.join(config["out_dir"], "trimmed_fastq", "{sample}_bc.fastq.gz"),
        index = config["star"]["index"],
        index_whitelist = _wl_index_effective(),
    output:
        bam = os.path.join(config["out_dir"], "mapped", "{sample}", "{sample}_Aligned.sortedByCoord.out.bam"),
        solo_dir = directory(os.path.join(config["out_dir"], "mapped", "{sample}", "{sample}_Solo.out")),
        star_logs = os.path.join(config["out_dir"], "mapped", "{sample}", "{sample}_Log.final.out"),
    params:
        out_prefix = lambda wc, output: os.path.join(config["out_dir"], "mapped", wc.sample, f"{wc.sample}_"),
        features = config.get("star", {}).get("features", "GeneFull"),
        cb_match = config.get("star", {}).get("cb_match", "EditDist_2"),
        demux_by = config.get("demux_by", "index"),

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

        wl_bc1 = _wl("bc1", ""),
        wl_bc2 = _wl("bc2", ""),
        wl_bc3 = _wl("bc3", ""),
        wl_index  = _wl_index_effective(),
    threads: config.get("threads", 8)
    log:
        os.path.join(config["out_dir"], "logs", "starsoloUMI_{sample}.log"),
    conda:
        "../envs/star.yaml"
    script:
        "../scripts/starsolo.sh"



rule format_starsolo:
    input:
        solo_dir=os.path.join(config["out_dir"], "mapped/{sample}/{sample}_Solo.out"),
    output:
        summary=os.path.join(
            config["out_dir"], "mapped/{sample}/{sample}_Solo.out/{sample}_Summary.csv"
        ),
        umi=os.path.join(
            config["out_dir"],
            "mapped/{sample}/{sample}_Solo.out/{sample}_UMIperCellSorted.txt",
        ),
        barcodes=os.path.join(
            config["out_dir"],
            "mapped/{sample}/{sample}_Solo.out/{sample}_Barcodes.stats",
        ),
        features=os.path.join(
            config["out_dir"],
            "mapped/{sample}/{sample}_Solo.out/{sample}_Features.stats",
        ),
    log:
        os.path.join(config["out_dir"], "logs/format_starsolo_{sample}.log"),
    threads: 1
    conda:
        "../envs/pandas.yaml"
    params:
        sample="{sample}",
    script:
        "../scripts/starsolo_to_multiqc.py"