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
    log:
        f"{OUTDIR}/logs/pairtools/parse_sort_{{sample}}_{{run}}.log"
    params:
        outdir = lambda wildcards, output: os.path.dirname(output.pairsam),
        prefix = "{sample}_{run}",
        # Compresses only pairtools sort's temporary chunks (never the output); lz4c is the
        # wrapper's own default and much faster than gzip, at the cost of larger temp files
        compress_program = "lz4c"
    threads:
        16
    container:
        CONTAINER_4DN
    shell:
        """
        mkdir -p {params.outdir} $(dirname {log})
        /usr/local/bin/run-pairsam-parse-sort.sh \
            {input.bam} \
            {input.chrom_sizes} \
            {params.outdir} \
            {params.prefix} \
            {threads} \
            {params.compress_program} > {log} 2>&1
        """

# -----------------------------------------------------------------------------
# 2. Sequencing Replicate Merging (Merged BEFORE duplicate removal)
# -----------------------------------------------------------------------------
rule pairsam_merge_seq_reps:
    input:
        pairsams = lambda wildcards: [f"{OUTDIR}/pairsam/{wildcards.sample}_{r}.sam.pairs.gz" for r in SAMPLE_RUNS[wildcards.sample]]
    output:
        merged_pairsam = temp(f"{OUTDIR}/pairsam_merged/{{sample}}.merged.sam.pairs.gz")
    log:
        f"{OUTDIR}/logs/pairtools/merge_{{sample}}.log"
    params:
        outprefix = lambda wildcards, output: output.merged_pairsam.replace('.merged.sam.pairs.gz', '')
    threads:
        8
    container:
        CONTAINER_4DN
    shell:
        """
        mkdir -p $(dirname {output.merged_pairsam}) $(dirname {log})
        /usr/local/bin/run-pairsam-merge.sh \
            {params.outprefix} \
            {threads} \
            {input.pairsams} > {log} 2>&1
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
    log:
        f"{OUTDIR}/logs/pairtools/markasdup_{{sample}}.log"
    params:
        outprefix = lambda wildcards, output: output.marked_pairsam.replace('.marked.sam.pairs.gz', '')
    # The wrapper passes no thread flags, but pairtools starts pbgzip helpers by default
    # (3 decompressing, 8 compressing); 12 cores lets them run side by side
    threads:
        12
    container:
        CONTAINER_4DN
    shell:
        """
        mkdir -p $(dirname {output.marked_pairsam}) $(dirname {log})
        /usr/local/bin/run-pairsam-markasdup.sh \
            {input.pairsam} \
            {params.outprefix} > {log} 2>&1
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
    log:
        f"{OUTDIR}/logs/pairtools/filter_{{sample}}.log"
    params:
        outprefix = lambda wildcards, output: output.dedup_pairs.replace('.dedup.pairs.gz', '')
    # The wrapper passes no thread flags, but each pairtools step starts pbgzip helpers by
    # default (3 decompressing, 8 compressing per output); 12 cores lets them run side by side
    threads:
        12
    # run-pairsam-filter.sh writes fixed-name temp.gz / temp1.gz into the working
    # directory; a per-job shadow directory stops concurrent samples overwriting them
    shadow:
        "minimal"
    container:
        CONTAINER_4DN
    shell:
        """
        mkdir -p $(dirname {output.dedup_pairs}) $(dirname {output.lossless_bam}) $(dirname {log})
        
        /usr/local/bin/run-pairsam-filter.sh \
            {input.marked_pairsam} \
            {params.outprefix} \
            {input.chrom_sizes} > {log} 2>&1
            
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
    log:
        f"{OUTDIR}/logs/pairtools/addfrag_{{sample}}.log"
    params:
        outprefix = lambda wildcards, output: output.ff_pairs.replace('.ff.pairs.gz', ''),
        re_arg = f"-r {config['restriction_file']}" if (config.get("restriction_file") and EXPERIMENT_TYPE not in ["micro_c", "dnase_hic", "chia_pet"]) else ""
    threads:
        2
    container:
        CONTAINER_4DN
    shell:
        """
        mkdir -p $(dirname {output.ff_pairs}) $(dirname {log})
        /usr/local/bin/run-addfrag2pairs.sh \
            -i {input.pairs} \
            -o {params.outprefix} \
            {params.re_arg} > {log} 2>&1
        """

# -----------------------------------------------------------------------------
# 6. Biological Replicate Merging (Official 4DN merge-pairs.cwl)
# -----------------------------------------------------------------------------
rule merge_bio_reps:
    input:
        pairs = lambda wildcards: [f"{OUTDIR}/pairs/{s}.ff.pairs.gz" for s in CONDITION_SAMPLES[wildcards.condition]],
        px2 = lambda wildcards: [f"{OUTDIR}/pairs/{s}.ff.pairs.gz.px2" for s in CONDITION_SAMPLES[wildcards.condition]]
    output:
        merged_ff_pairs = f"{OUTDIR}/pairs_merged/{{condition}}.merged.ff.pairs.gz",
        merged_ff_px2 = f"{OUTDIR}/pairs_merged/{{condition}}.merged.ff.pairs.gz.px2"
    log:
        f"{OUTDIR}/logs/pairtools/merge_bioreps_{{condition}}.log"
    params:
        outprefix = lambda wildcards, output: output.merged_ff_pairs.replace('.pairs.gz', '')
    threads:
        2
    # run-merge-pairs.sh creates fixed-name FIFOs pp.1, pp.2, ... in the working
    # directory; a per-job shadow directory stops concurrent conditions sharing them
    shadow:
        "minimal"
    container:
        CONTAINER_4DN
    shell:
        """
        mkdir -p $(dirname {output.merged_ff_pairs}) $(dirname {log})
        /usr/local/bin/run-merge-pairs.sh \
            {params.outprefix} \
            {input.pairs} > {log} 2>&1
        """
