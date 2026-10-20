# =============================================================================
# 06_provenance.smk: Record Pipeline Provenance & Reproducibility Metadata
# =============================================================================

localrules: record_provenance

rule record_provenance:
    input:
        pairs = expand(f"{OUTDIR}/pairs/{{sample}}.ff.pairs.gz", sample=SAMPLES),
        mcool = expand(f"{OUTDIR}/mcool/{{sample}}.mcool", sample=SAMPLES),
        hic = expand(f"{OUTDIR}/hic/{{sample}}.hic", sample=SAMPLES) if config.get("matrix", {}).get("generate_hic", False) else []
    output:
        provenance_json = f"{OUTDIR}/provenance.json"
    log:
        f"{OUTDIR}/logs/provenance/provenance.log"
    run:
        import json
        import datetime
        import os

        fasta_path = config.get("reference", {}).get("fasta", "")
        fasta_size = os.path.getsize(fasta_path) if fasta_path and os.path.exists(fasta_path) else "N/A"

        meta = {
            "pipeline": {
                "name": "4DN v43 Hi-C / Micro-C Snakemake Pipeline",
                "specification_url": "https://data.4dnucleome.org/resources/data-analysis/hi_c-processing-pipeline",
                "reference_repo": "https://github.com/4dn-dcic/docker-4dn-hic",
                "target_release": "v43",
                "execution_time_utc": datetime.datetime.utcnow().isoformat()
            },
            "containers": {
                "4dn_hic": CONTAINER_4DN,
                "multiqc": CONTAINER_MULTIQC
            },
            "envmodules": {
                "fastqc": MODULE_FASTQC,
                "cooltools": MODULE_COOLTOOLS
            },
            "experiment": {
                "experiment_type": EXPERIMENT_TYPE,
                "restriction_file": config.get("restriction_file", None),
                "save_bam": SAVE_BAM
            },
            "reference": {
                "assembly": config.get("reference", {}).get("assembly", "unspecified"),
                "fasta": fasta_path,
                "fasta_size_bytes": fasta_size,
                "chrom_sizes": get_chrom_sizes()
            },
            "matrix": {
                "base_resolution": BASE_RES,
                "resolutions": RESOLUTIONS_LIST,
                "balance": config.get("matrix", {}).get("balance", True),
                "add_juicer_norms_to_mcool": config.get("matrix", {}).get("add_juicer_norms_to_mcool", True)
            },
            "conditions": CONDITIONS,
            "condition_samples": CONDITION_SAMPLES,
            "multi_rep_conditions": MULTI_REP_CONDITIONS,
            "samples": SAMPLES,
            "sample_runs": SAMPLE_RUNS
        }

        out_file = str(output[0])
        os.makedirs(os.path.dirname(out_file), exist_ok=True)
        with open(out_file, "w") as f:
            json.dump(meta, f, indent=2)
