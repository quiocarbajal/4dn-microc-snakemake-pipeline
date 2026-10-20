# =============================================================================
# 02_alignment.smk: 4DN BWA-MEM Alignment (-SP5M)
# Runs bwa mem -SP5M directly inside the official 4DN v43 container. The 4DN
# run-bwa-mem.sh wrapper (used for .tar.gz indexes) is not used: it unpacks the
# index and writes fixed-name fastq1/fastq2 copies into the shared output
# directory, so concurrent samples could read each other's reads.
# =============================================================================

rule bwa_mem:
    input:
        fastq1 = get_fastq1,
        fastq2 = get_fastq2,
        bwa_index = get_bwa_index()
    output:
        bam = temp(f"{OUTDIR}/bam/{{sample}}_{{run}}.bam") if not SAVE_BAM else f"{OUTDIR}/bam/{{sample}}_{{run}}.bam"
    log:
        f"{OUTDIR}/logs/bwa_mem/{{sample}}_{{run}}.log"
    params:
        outdir = lambda wildcards, output: os.path.dirname(output.bam)
    threads:
        8
    container:
        CONTAINER_4DN
    shell:
        """
        mkdir -p {params.outdir} $(dirname {log})

        # Pre-indexed directory or index prefix (.tar.gz archives are rejected in get_bwa_index)
        if [[ -d "{input.bwa_index}" ]]; then
            idx_prefix=$(ls -1 {input.bwa_index}/*.bwt 2>/dev/null | head -1 | sed 's/[.]bwt//g')
            if [[ -z "$idx_prefix" ]]; then
                idx_prefix="{input.bwa_index}"
            fi
        else
            idx_prefix="{input.bwa_index}"
        fi

        # Execute 4DN standard bwa mem mapping with -SP5M
        bwa mem -t {threads} -SP5M "$idx_prefix" {input.fastq1} {input.fastq2} 2> {log} | samtools view -Shb - > {output.bam} 2>> {log}
        """
