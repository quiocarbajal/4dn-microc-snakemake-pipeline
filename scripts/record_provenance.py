import json
import datetime
import hashlib
import os

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

meta = {
    "pipeline": {
        "name": "4DN v43 Hi-C / Micro-C Snakemake Pipeline",
        "specification_url": "https://data.4dnucleome.org/resources/data-analysis/hi_c-processing-pipeline",
        "reference_repo": "https://github.com/4dn-dcic/docker-4dn-hic",
        "target_release": "v43",
        "execution_time_utc": datetime.datetime.utcnow().isoformat()
    },
    "containers": snakemake.params.containers,
    "experiment": snakemake.params.experiment,
    "reference": {
        "assembly": snakemake.params.reference.get("assembly", "unspecified"),
        "fasta": snakemake.params.reference.get("fasta", ""),
        "fasta_sha256": calculate_checksum(snakemake.params.reference.get("fasta", "")),
        "chrom_sizes": snakemake.params.chrom_sizes
    },
    "matrix": snakemake.params.matrix,
    "conditions": snakemake.params.conditions,
    "condition_samples": snakemake.params.condition_samples,
    "multi_rep_conditions": snakemake.params.multi_rep_conditions,
    "samples": snakemake.params.samples,
    "sample_runs": snakemake.params.sample_runs
}

os.makedirs(os.path.dirname(snakemake.output.provenance_json), exist_ok=True)
with open(snakemake.output.provenance_json, "w") as f:
    json.dump(meta, f, indent=2)
