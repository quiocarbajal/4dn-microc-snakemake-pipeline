# =============================================================================
# 00_reference.smk: Reference File Preparation Rules
# =============================================================================

rule generate_chrom_sizes:
    input:
        src = get_chrom_sizes_source()
    output:
        chrom_sizes = f"{OUTDIR}/reference/chrom.sizes"
    log:
        f"{OUTDIR}/logs/reference/generate_chrom_sizes.log"
    container:
        CONTAINER_4DN
    shell:
        """
        mkdir -p $(dirname {output.chrom_sizes}) $(dirname {log})
        if [[ "{input.src}" == *.ann ]]; then
            awk 'NR>1 {{ if (NR%2==0) {{ name=$2 }} else {{ print name "\\t" $2 }} }}' {input.src} > {output.chrom_sizes} 2> {log}
        elif [[ "{input.src}" == *.fai ]]; then
            cut -f1,2 {input.src} > {output.chrom_sizes} 2> {log}
        else
            samtools faidx {input.src} > {log} 2>&1
            cut -f1,2 {input.src}.fai > {output.chrom_sizes} 2>> {log}
        fi
        """

rule build_bwa_index:
    input:
        fasta = lambda wildcards: config.get("reference", {}).get("fasta", "") or "reference/genome.fasta"
    output:
        bwa_index_dir = directory(f"{OUTDIR}/reference/bwa_index")
    log:
        f"{OUTDIR}/logs/reference/build_bwa_index.log"
    params:
        prefix = lambda wildcards, output: f"{output.bwa_index_dir}/genome"
    threads:
        8
    container:
        CONTAINER_4DN
    shell:
        """
        mkdir -p {output.bwa_index_dir} $(dirname {log})
        bwa index -p {params.prefix} {input.fasta} > {log} 2>&1
        """
