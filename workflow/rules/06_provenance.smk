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
        import subprocess

        fasta_path = config.get("reference", {}).get("fasta", "")
        fasta_size = os.path.getsize(fasta_path) if fasta_path and os.path.exists(fasta_path) else "N/A"

        # Exact pipeline code used for this run: commit, nearest tag, and any files in the
        # run folder that differ from that commit (e.g. a dataset-specific config/samples.tsv)
        def run_git(git_arguments):
            try:
                return subprocess.run(["git", "-C", workflow.basedir] + git_arguments, capture_output=True, text=True, check=True).stdout.rstrip("\n")
            except Exception:
                return "unavailable"

        git_status_lines = run_git(["status", "--porcelain", "--untracked-files=no"])

        meta = {
            "pipeline": {
                "name": "4DN v43 Hi-C / Micro-C Snakemake Pipeline",
                "specification_url": "https://data.4dnucleome.org/resources/data-analysis/hi_c-processing-pipeline",
                "reference_repo": "https://github.com/4dn-dcic/docker-4dn-hic",
                "target_release": "v43",
                "execution_time_utc": datetime.datetime.utcnow().isoformat()
            },
            "pipeline_code": {
                "git_commit": run_git(["rev-parse", "HEAD"]),
                "git_describe": run_git(["describe", "--tags", "--always", "--dirty"]),
                "git_remote_origin": run_git(["remote", "get-url", "origin"]),
                "files_differing_from_commit": git_status_lines.splitlines() if git_status_lines not in ("", "unavailable") else ([] if git_status_lines == "" else "unavailable")
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
