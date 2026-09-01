# =============================================================================
# 03_pairtools.smk: 4DN Pairtools Contact Processing & Filtering
# Uses official 4DN v43 scripts inside official container
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Pair Parsing & Coordinate Sorting
# -----------------------------------------------------------------------------
rule pairsam_parse_sort:
    input:
        bam = f"{OUTDIR}/bam/{{sample}}_{{run}}.bam",
        chrom_sizes = get_chrom_sizes()
    output:
        pairsam = temp(f"{OUTDIR}/pairsam/{{sample}}_{{run}}.sam.pairs.gz")
    params:
        outdir = f"{OUTDIR}/pairsam",
        prefix = "{sample}_{run}",
        compress_program = "gzip"
    threads:
        16
    container:
        CONTAINER_4DN
    shell:
        """
        mkdir -p {params.outdir}
        /usr/local/bin/run-pairsam-parse-sort.sh \
            {input.bam} \
            {input.chrom_sizes} \
            {params.outdir} \
            {params.prefix} \
            {threads} \
            {params.compress_program}
        """

# -----------------------------------------------------------------------------
# 2. Sequencing Replicate Merging (Merged BEFORE duplicate removal)
# -----------------------------------------------------------------------------
def get_pairsam_runs_for_sample(wildcards):
    runs = SAMPLE_RUNS[wildcards.sample]
    return [f"{OUTDIR}/pairsam/{wildcards.sample}_{r}.sam.pairs.gz" for r in runs]

rule pairsam_merge_seq_reps:
    input:
        pairsams = get_pairsam_runs_for_sample
    output:
        merged_pairsam = temp(f"{OUTDIR}/pairsam_merged/{{sample}}.merged.sam.pairs.gz")
    params:
        outprefix = f"{OUTDIR}/pairsam_merged/{{sample}}"
    threads:
        8
    container:
        CONTAINER_4DN
    shell:
        """
        mkdir -p $(dirname {output.merged_pairsam})
        /usr/local/bin/run-pairsam-merge.sh \
            {params.outprefix} \
            {threads} \
            {input.pairsams}
        """

# -----------------------------------------------------------------------------
# 3. Duplicate Marking
# -----------------------------------------------------------------------------
rule pairsam_markasdup:
    input:
        pairsam = get_pairsam_for_dedup
    output:
        marked_pairsam = temp(f"{OUTDIR}/pairsam_marked/{{sample}}.marked.sam.pairs.gz"),
        marked_px2 = temp(f"{OUTDIR}/pairsam_marked/{{sample}}.marked.sam.pairs.gz.px2")
    params:
        outprefix = f"{OUTDIR}/pairsam_marked/{{sample}}"
    threads:
        4
    container:
        CONTAINER_4DN
    shell:
        """
        mkdir -p $(dirname {output.marked_pairsam})
        /usr/local/bin/run-pairsam-markasdup.sh \
            {input.pairsam} \
            {params.outprefix}
        """

# -----------------------------------------------------------------------------
# 4. Contact Filtering & Valid Pair Extraction
# -----------------------------------------------------------------------------
rule pairsam_filter:
    input:
        marked_pairsam = f"{OUTDIR}/pairsam_marked/{{sample}}.marked.sam.pairs.gz",
        chrom_sizes = get_chrom_sizes()
    output:
        dedup_pairs = f"{OUTDIR}/pairs_filtered/{{sample}}.dedup.pairs.gz",
        dedup_px2 = f"{OUTDIR}/pairs_filtered/{{sample}}.dedup.pairs.gz.px2",
        lossless_bam = f"{OUTDIR}/bam_lossless/{{sample}}.lossless.bam" if SAVE_BAM else temp(f"{OUTDIR}/bam_lossless/{{sample}}.lossless.bam"),
        unmapped_pairsam = temp(f"{OUTDIR}/pairs_filtered/{{sample}}.unmapped.sam.pairs.gz")
    params:
        outprefix = f"{OUTDIR}/pairs_filtered/{{sample}}"
    threads:
        4
    container:
        CONTAINER_4DN
    shell:
        """
        mkdir -p $(dirname {output.dedup_pairs})
        mkdir -p $(dirname {output.lossless_bam})
        
        /usr/local/bin/run-pairsam-filter.sh \
            {input.marked_pairsam} \
            {params.outprefix} \
            {input.chrom_sizes}
            
        # Move lossless BAM to dedicated directory if retained
        if [ -f "{params.outprefix}.lossless.bam" ]; then
            mv {params.outprefix}.lossless.bam {output.lossless_bam}
        fi
        """

# -----------------------------------------------------------------------------
# 5. Fragment Assignment (Micro-C restriction-free handling)
# -----------------------------------------------------------------------------
rule addfrag_or_pass:
    input:
        pairs = f"{OUTDIR}/pairs_filtered/{{sample}}.dedup.pairs.gz",
        px2 = f"{OUTDIR}/pairs_filtered/{{sample}}.dedup.pairs.gz.px2"
    output:
        ff_pairs = f"{OUTDIR}/pairs/{{sample}}.ff.pairs.gz",
        ff_px2 = f"{OUTDIR}/pairs/{{sample}}.ff.pairs.gz.px2"
    params:
        outprefix = f"{OUTDIR}/pairs/{{sample}}",
        re_arg = f"-r {config['restriction_file']}" if (config.get("restriction_file") and EXPERIMENT_TYPE not in ["micro_c", "dnase_hic", "chia_pet"]) else ""
    threads:
        2
    container:
        CONTAINER_4DN
    shell:
        """
        mkdir -p $(dirname {output.ff_pairs})
        /usr/local/bin/run-addfrag2pairs.sh \
            -i {input.pairs} \
            -o {params.outprefix} \
            {params.re_arg}
        """
