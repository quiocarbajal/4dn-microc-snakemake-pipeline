# =============================================================================
# 06_provenance.smk: Record Pipeline Provenance & Reproducibility Metadata
# =============================================================================

import json
import datetime
import hashlib

def calculate_checksum(filepath):
    if not filepath or not os.path.exists(filepath):
        return "N/A"
    try:
        h = hashlib.sha256()
        with open(filepath, "rb") as f:
            while chunk := f.read(8192):
                h.update(chunk)
        return h.hexdigest()
    except Exception:
        return "N/A"

rule record_provenance:
    input:
        pairs = expand(f"{OUTDIR}/pairs/{{sample}}.ff.pairs.gz", sample=SAMPLES),
        mcool = expand(f"{OUTDIR}/mcool/{{sample}}.mcool", sample=SAMPLES),
        hic = expand(f"{OUTDIR}/hic/{{sample}}.hic", sample=SAMPLES) if config.get("matrix", {}).get("generate_hic", False) else []
    output:
        provenance_json = f"{OUTDIR}/provenance.json"
    run:
        meta = {
            "pipeline": {
                "name": "4DN v43 Hi-C / Micro-C Snakemake Pipeline",
                "specification_url": "https://data.4dnucleome.org/resources/data-analysis/hi_c-processing-pipeline",
                "reference_repo": "https://github.com/4dn-dcic/docker-4dn-hic",
                "target_release": "v43",
                "snakemake_version": str(snakemake.__version__) if hasattr(snakemake, '__version__') else "9.25.1",
                "execution_time_utc": datetime.datetime.utcnow().isoformat()
            },
            "containers": {
                "4dn_hic": CONTAINER_4DN,
                "cooltools": CONTAINER_COOLTOOLS,
                "multiqc": CONTAINER_MULTIQC
            },
            "experiment": {
                "experiment_type": EXPERIMENT_TYPE,
                "restriction_file": config.get("restriction_file", None),
                "save_bam": SAVE_BAM
            },
            "reference": {
                "assembly": config.get("reference", {}).get("assembly", "unspecified"),
                "fasta": config.get("reference", {}).get("fasta", ""),
                "fasta_sha256": calculate_checksum(config.get("reference", {}).get("fasta", "")),
                "chrom_sizes": get_chrom_sizes()
            },
            "matrix": {
                "base_resolution": BASE_RES,
                "resolutions": RESOLUTIONS_LIST,
                "balance": config.get("matrix", {}).get("balance", True),
                "add_juicer_norms_to_mcool": config.get("matrix", {}).get("add_juicer_norms_to_mcool", True)
            },
            "samples": SAMPLES,
            "sample_runs": SAMPLE_RUNS
        }
        with open(output.provenance_json, "w") as f:
            json.dump(meta, f, indent=2)
