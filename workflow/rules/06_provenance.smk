# =============================================================================
# 06_provenance.smk: Record Pipeline Provenance & Reproducibility Metadata
# =============================================================================

rule record_provenance:
    input:
        pairs = expand(f"{OUTDIR}/pairs/{{sample}}.ff.pairs.gz", sample=SAMPLES),
        mcool = expand(f"{OUTDIR}/mcool/{{sample}}.mcool", sample=SAMPLES),
        hic = expand(f"{OUTDIR}/hic/{{sample}}.hic", sample=SAMPLES) if config.get("matrix", {}).get("generate_hic", False) else []
    output:
        provenance_json = f"{OUTDIR}/provenance.json"
    log:
        f"{OUTDIR}/logs/provenance/provenance.log"
    params:
        containers = {
            "4dn_hic": CONTAINER_4DN,
            "cooltools": CONTAINER_COOLTOOLS,
            "multiqc": CONTAINER_MULTIQC
        },
        experiment = {
            "experiment_type": EXPERIMENT_TYPE,
            "restriction_file": config.get("restriction_file", None),
            "save_bam": SAVE_BAM
        },
        reference = config.get("reference", {}),
        chrom_sizes = get_chrom_sizes(),
        matrix = {
            "base_resolution": BASE_RES,
            "resolutions": RESOLUTIONS_LIST,
            "balance": config.get("matrix", {}).get("balance", True),
            "add_juicer_norms_to_mcool": config.get("matrix", {}).get("add_juicer_norms_to_mcool", True)
        },
        samples = SAMPLES,
        sample_runs = SAMPLE_RUNS
    script:
        "../../scripts/record_provenance.py"
