# =============================================================================
# 02_alignment.smk: 4DN BWA-MEM Alignment (-SP5M)
# Uses official 4DN v43 run-bwa-mem.sh inside official container
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
        outdir = lambda wildcards, output: os.path.dirname(output.bam),
        prefix = "{sample}_{run}"
    threads:
        8
    container:
        CONTAINER_4DN
    shell:
        """
        mkdir -p {params.outdir} $(dirname {log})
        
        # Check if bwa_index is a tar.gz / tgz archive (standard 4DN wrapper format)
        if [[ -f "{input.bwa_index}" && ( "{input.bwa_index}" == *.tar.gz || "{input.bwa_index}" == *.tgz ) ]]; then
            /usr/local/bin/run-bwa-mem.sh \
                {input.fastq1} \
                {input.fastq2} \
                {input.bwa_index} \
                {params.outdir} \
                {params.prefix} \
                {threads} > {log} 2>&1
        else
            # Pre-indexed directory or index prefix
            if [[ -d "{input.bwa_index}" ]]; then
                idx_prefix=$(ls -1 {input.bwa_index}/*.bwt 2>/dev/null | head -1 | sed 's/\.bwt//g')
                if [[ -z "$idx_prefix" ]]; then
                    idx_prefix="{input.bwa_index}"
                fi
            else
                idx_prefix="{input.bwa_index}"
            fi
            
            # Execute 4DN standard bwa mem mapping with -SP5M
            bwa mem -t {threads} -SP5M "$idx_prefix" {input.fastq1} {input.fastq2} 2> {log} | samtools view -Shb - > {output.bam} 2>> {log}
        fi
        """
