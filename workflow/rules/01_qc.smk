# =============================================================================
# 01_qc.smk: Quality Control (FastQC and pairsqc)
# Uses official 4DN v43 scripts inside official container
# =============================================================================

# -----------------------------------------------------------------------------
# 1. FastQC Read Quality Control
# -----------------------------------------------------------------------------
rule fastqc_r1:
    input:
        fastq = get_fastq1
    output:
        zip = f"{OUTDIR}/qc/fastqc/{{sample}}_{{run}}_R1_fastqc.zip",
        html = f"{OUTDIR}/qc/fastqc/{{sample}}_{{run}}_R1_fastqc.html"
    log:
        f"{OUTDIR}/logs/fastqc/{{sample}}_{{run}}_R1.log"
    params:
        outdir = lambda wildcards, output: f"{os.path.dirname(output.zip)}/tmp_{wildcards.sample}_{wildcards.run}_R1"
    threads:
        4
    container:
        CONTAINER_4DN
    shell:
        """
        mkdir -p {params.outdir} $(dirname {log})
        /usr/local/bin/run-fastqc.sh {input.fastq} {threads} {params.outdir} > {log} 2>&1
        
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
    log:
        f"{OUTDIR}/logs/fastqc/{{sample}}_{{run}}_R2.log"
    params:
        outdir = lambda wildcards, output: f"{os.path.dirname(output.zip)}/tmp_{wildcards.sample}_{wildcards.run}_R2"
    threads:
        4
    container:
        CONTAINER_4DN
    shell:
        """
        mkdir -p {params.outdir} $(dirname {log})
        /usr/local/bin/run-fastqc.sh {input.fastq} {threads} {params.outdir} > {log} 2>&1
        
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
    log:
        f"{OUTDIR}/logs/pairsqc/{{sample}}.log"
    params:
        sample_name = "{sample}",
        enzyme = config.get("qc", {}).get("pairsqc_enzyme", "4"),
        max_dist = config.get("qc", {}).get("pairsqc_max_distance", 8.4),
        outdir = lambda wildcards, output: os.path.dirname(output.report_zip)
    threads:
        4
    container:
        CONTAINER_4DN
    shell:
        """
        mkdir -p {params.outdir} $(dirname {log})
        /usr/local/bin/run-pairsqc-single.sh \
            {input.pairs} \
            {input.chrom_sizes} \
            {params.sample_name} \
            {params.enzyme} \
            {params.outdir} \
            {params.max_dist} > {log} 2>&1
        """
