# =============================================================================
# 01_qc.smk: Quality Control (FastQC and pairsqc)
# Uses official 4DN v43 scripts inside official container
# =============================================================================

def get_fastq1(wildcards):
    return SAMPLE_RUN_MAP[(wildcards.sample, wildcards.run)][0]

def get_fastq2(wildcards):
    return SAMPLE_RUN_MAP[(wildcards.sample, wildcards.run)][1]

# -----------------------------------------------------------------------------
# 1. FastQC Read Quality Control
# -----------------------------------------------------------------------------
rule fastqc_r1:
    input:
        fastq = get_fastq1
    output:
        zip = f"{OUTDIR}/qc/fastqc/{{sample}}_{{run}}_R1_fastqc.zip",
        html = f"{OUTDIR}/qc/fastqc/{{sample}}_{{run}}_R1_fastqc.html"
    params:
        outdir = f"{OUTDIR}/qc/fastqc/tmp_{{sample}}_{{run}}_R1"
    threads:
        4
    container:
        CONTAINER_4DN
    shell:
        """
        mkdir -p {params.outdir}
        /usr/local/bin/run-fastqc.sh {input.fastq} {threads} {params.outdir}
        
        # Standardize output filenames
        mv {params.outdir}/*_fastqc.zip {output.zip}
        mv {params.outdir}/*_fastqc.html {output.html}
        rm -rf {params.outdir}
        """

rule fastqc_r2:
    input:
        fastq = get_fastq2
    output:
        zip = f"{OUTDIR}/qc/fastqc/{{sample}}_{{run}}_R2_fastqc.zip",
        html = f"{OUTDIR}/qc/fastqc/{{sample}}_{{run}}_R2_fastqc.html"
    params:
        outdir = f"{OUTDIR}/qc/fastqc/tmp_{{sample}}_{{run}}_R2"
    threads:
        4
    container:
        CONTAINER_4DN
    shell:
        """
        mkdir -p {params.outdir}
        /usr/local/bin/run-fastqc.sh {input.fastq} {threads} {params.outdir}
        
        # Standardize output filenames
        mv {params.outdir}/*_fastqc.zip {output.zip}
        mv {params.outdir}/*_fastqc.html {output.html}
        rm -rf {params.outdir}
        """

# -----------------------------------------------------------------------------
# 2. 4DN pairsqc Contact Distribution QC
# -----------------------------------------------------------------------------
rule pairsqc:
    input:
        pairs = f"{OUTDIR}/pairs/{{sample}}.ff.pairs.gz",
        px2 = f"{OUTDIR}/pairs/{{sample}}.ff.pairs.gz.px2",
        chrom_sizes = get_chrom_sizes()
    output:
        report_zip = f"{OUTDIR}/qc/pairsqc/{{sample}}_report.zip"
    params:
        sample_name = "{sample}",
        enzyme = config.get("qc", {}).get("pairsqc_enzyme", "4"),
        max_dist = config.get("qc", {}).get("pairsqc_max_distance", 8.4),
        outdir = f"{OUTDIR}/qc/pairsqc"
    threads:
        4
    container:
        CONTAINER_4DN
    shell:
        """
        mkdir -p {params.outdir}
        /usr/local/bin/run-pairsqc-single.sh \
            {input.pairs} \
            {input.chrom_sizes} \
            {params.sample_name} \
            {params.enzyme} \
            {params.outdir} \
            {params.max_dist}
        """
