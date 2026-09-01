import os
import sys
import pandas as pd
from snakemake.utils import validate

# -----------------------------------------------------------------------------
# 1. Configuration & Sample Sheet Loading
# -----------------------------------------------------------------------------
configfile: "config/config.yaml"

# Load sample sheet
samples_file = config.get("samples", "config/samples.tsv")
if not os.path.exists(samples_file):
    sys.exit(f"Error: Sample sheet '{samples_file}' does not exist.")

samples_df = pd.read_csv(samples_file, sep="\t", dtype=str, comment="#").dropna(how="all")

# Validate required columns
required_cols = {"sample", "run", "fastq1", "fastq2"}
missing_cols = required_cols - set(samples_df.columns)
if missing_cols:
    sys.exit(f"Error: Sample sheet '{samples_file}' is missing required columns: {missing_cols}")

# Support condition/group column for biological replicates (defaults to sample if absent)
if "condition" not in samples_df.columns and "group" in samples_df.columns:
    samples_df["condition"] = samples_df["group"]
elif "condition" not in samples_df.columns:
    samples_df["condition"] = samples_df["sample"]

# Extract unique condition, sample, and run identifiers
CONDITIONS = sorted(list(samples_df["condition"].unique()))
SAMPLES = sorted(list(samples_df["sample"].unique()))
RUNS = sorted(list(samples_df["run"].unique()))

# Map condition -> [sample1, sample2, ...] and identify multi-replicate conditions
CONDITION_SAMPLES = {}
for _, row in samples_df.iterrows():
    c = row["condition"]
    s = row["sample"]
    if c not in CONDITION_SAMPLES:
        CONDITION_SAMPLES[c] = []
    if s not in CONDITION_SAMPLES[c]:
        CONDITION_SAMPLES[c].append(s)

MULTI_REP_CONDITIONS = sorted([c for c, s_list in CONDITION_SAMPLES.items() if len(s_list) > 1])

# Create lookup dictionary mapping (sample, run) -> (fastq1, fastq2)
SAMPLE_RUNS = {}
SAMPLE_RUN_MAP = {}
SAMPLE_CONDITION = {}
for _, row in samples_df.iterrows():
    c = row["condition"]
    s = row["sample"]
    r = row["run"]
    f1 = row["fastq1"]
    f2 = row["fastq2"]
    SAMPLE_CONDITION[s] = c
    if s not in SAMPLE_RUNS:
        SAMPLE_RUNS[s] = []
    SAMPLE_RUNS[s].append(r)
    SAMPLE_RUN_MAP[(s, r)] = (f1, f2)

# Global variables and paths
OUTDIR = config.get("outdir", "results")
CONTAINER_4DN = config.get("containers", {}).get("4dn_hic", "docker://duplexa/4dn-hic:v43")
CONTAINER_COOLTOOLS = config.get("containers", {}).get("cooltools", "docker://open2c/cooltools:v0.6.1")
CONTAINER_MULTIQC = config.get("containers", {}).get("multiqc", "docker://multiqc/multiqc:v1.27")

SAVE_BAM = config.get("save_bam", False)
EXPERIMENT_TYPE = config.get("experiment_type", "micro_c").lower()

# Resolutions formatted as comma-separated string for 4DN scripts
RESOLUTIONS_LIST = config.get("matrix", {}).get("resolutions", [
    1000, 2000, 5000, 10000, 25000, 50000, 100000, 250000, 500000, 1000000, 2500000, 5000000, 10000000
])
RESOLUTIONS_STR = ",".join(str(r) for r in RESOLUTIONS_LIST)
BASE_RES = config.get("matrix", {}).get("base_resolution", 1000)

# -----------------------------------------------------------------------------
# 2. Reference & Chromosome Sizes Resolution
# -----------------------------------------------------------------------------
def find_reference_file(extensions):
    fa = config.get("reference", {}).get("fasta", "")
    if fa and os.path.exists(fa):
        return fa
    
    bwa_idx = config.get("reference", {}).get("bwa_index", "") or config.get("reference", {}).get("bwa_index_tgz", "")
    if bwa_idx and os.path.isdir(bwa_idx):
        try:
            for f in os.listdir(bwa_idx):
                for ext in extensions:
                    if f.endswith(ext):
                        return os.path.join(bwa_idx, f)
            parent = os.path.dirname(os.path.abspath(bwa_idx))
            if os.path.isdir(parent):
                for f in os.listdir(parent):
                    for ext in extensions:
                        if f.endswith(ext):
                            return os.path.join(parent, f)
        except Exception:
            pass
    return ""

def get_chrom_sizes():
    cs = config.get("reference", {}).get("chrom_sizes", "")
    if cs and os.path.exists(cs):
        return cs
    
    found_cs = find_reference_file([".chrom.sizes", ".sizes"])
    if found_cs:
        return found_cs
        
    return f"{OUTDIR}/reference/chrom.sizes"

def get_bwa_index():
    idx = config.get("reference", {}).get("bwa_index", "") or config.get("reference", {}).get("bwa_index_tgz", "")
    if idx:
        return idx
    return f"{OUTDIR}/reference/bwa_index.tar.gz"

def get_chrom_sizes_source():
    fa = config.get("reference", {}).get("fasta", "")
    if fa and os.path.exists(fa):
        return fa
    # Check for .ann, .fai, or FASTA in bwa directory / parent
    found = find_reference_file([".ann", ".fai", ".fasta", ".fa", ".fna"])
    if found:
        return found
    return "reference/genome.fasta"

def get_fastq1(wildcards):
    return SAMPLE_RUN_MAP[(wildcards.sample, wildcards.run)][0]

def get_fastq2(wildcards):
    return SAMPLE_RUN_MAP[(wildcards.sample, wildcards.run)][1]

# Helper to determine pairsam input for duplicate marking
def get_pairsam_for_dedup(wildcards):
    runs = SAMPLE_RUNS[wildcards.sample]
    if len(runs) > 1:
        return f"{OUTDIR}/pairsam_merged/{wildcards.sample}.merged.sam.pairs.gz"
    else:
        r = runs[0]
        return f"{OUTDIR}/pairsam/{wildcards.sample}_{r}.sam.pairs.gz"

