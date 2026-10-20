# Micro-C & Hi-C Processing Pipeline (4DN v43 Official Specification)

A production-ready **Snakemake 9** workflow reproducing the official **4D Nucleome (4DN) Hi-C and Micro-C Processing Pipeline** (Release `v43`) for high-performance computing (HPC) clusters managed by **Slurm** using **Apptainer** (Singularity).

- **Official 4DN Pipeline Documentation**: [4DN Data Portal Hi-C Processing Pipeline](https://data.4dnucleome.org/resources/data-analysis/hi_c-processing-pipeline)
- **Official 4DN Docker/CWL Repository**: [4dn-dcic/docker-4dn-hic (v43)](https://github.com/4dn-dcic/docker-4dn-hic)
- **Official Container Image**: `duplexa/4dn-hic:v43` (Docker Hub)

---

## Key Features

1. **Faithful 4DN v43 Reproduction**:
   - Executes the exact, unmodified 4DN wrapper scripts directly from the official container (`run-pairsam-parse-sort.sh`, `run-pairsam-merge.sh`, `run-pairsam-markasdup.sh`, `run-pairsam-filter.sh`, `run-cooler.sh`, `run-cool2multirescool.sh`, `run-juicebox-pre.sh`, `run-add-hicnormvector-to-mcool.sh`, `run-pairsqc-single.sh`). Alignment runs the same `bwa mem -SP5M` command directly instead of `run-bwa-mem.sh`, which writes fixed-name files into a shared directory; BWA index archives (`.tar.gz`) are therefore not supported.
2. **Micro-C & High-Resolution Support**:
   - Configurable protocol modes (`micro_c`, `dnase_hic`, `in_situ_hic`).
   - Restriction-fragment-free contact extraction and ICE matrix balancing for Micro-C.
   - Multiresolution contact matrices from **1 kb down to 10 Mb** (`.mcool` and `.hic`).
3. **Dual Normalization in `.mcool`**:
   - Generates ICE balanced `.mcool` files and injects Juicer normalization vectors (`VC`, `VC_SQRT`, `KR`) from `.hic` via `run-add-hicnormvector-to-mcool.sh` for seamless HiGlass and Juicebox compatibility.
4. **HPC & Slurm Native (No Root/Admin Privileges Required)**:
   - Automated container pulling and caching via **Apptainer 1.3.6**.
   - Per-rule dynamic Slurm resource allocations and automatic memory escalation.
   - Built-in **Snakemake 9** Slurm profile (`workflow/profiles/slurm`).
5. **Storage Efficient**:
   - Configurable BAM retention (`save_bam: false` by default, saving hundreds of GBs of temporary I/O while preserving contact pairs and `.mcool` files).
6. **Separated Downstream Analysis Layer**:
   - Optional downstream feature calling using `cooltools` ($P(s)$ contact decay curves, multi-scale insulation score TAD boundaries, and A/B compartments).

---

## Workflow DAG & 4DN Architecture

```text
               FASTQ R1 / R2
                     |
         +-----------+-----------+
         |                       |
   [run-fastqc.sh]        [bwa mem -SP5M]
         |                       |
      FastQC Report             BAM
                                 |
                     [run-pairsam-parse-sort.sh]
                                 |
                         .sam.pairs.gz
                                 |
                     [run-pairsam-merge.sh] (Sequencing Replicates)
                                 |
                     [run-pairsam-markasdup.sh]
                                 |
                     [run-pairsam-filter.sh] (UU, UC, CU valid pairs)
                                 |
                         .dedup.pairs.gz
                                 |
                     [run-addfrag2pairs.sh] (Micro-C pass-through)
                                 |
                          .ff.pairs.gz + .px2
                                 |
         +-----------------------+-----------------------+
         |                                               |
  [run-pairsqc.sh]                              [run-cooler.sh] (1kb base res)
         |                                               |
   QC Report Zip                                     .cool file
                                                         |
         +-----------------------------------------------+
         |                                               |
 [run-juicebox-pre.sh]                         [run-cool2multirescool.sh]
         |                                               |
     .hic file                                     .multires.cool (ICE)
         |                                               |
         +-----------------------+-----------------------+
                                 |
               [run-add-hicnormvector-to-mcool.sh]
                                 |
                        Final .mcool (ICE + Juicer Norms)
                                 |
         +-----------------------+-----------------------+  (Optional Downstream)
         |                       |                       |
   [cooltools P(s)]     [cooltools Insulation]  [cooltools Compartments]
```

---

## Directory Structure

```text
microC_pipeline/
├── README.md                          # Pipeline documentation
├── COMMANDS.txt                       # Direct copy-pasteable command cheat sheet
├── initial_README.md                  # Initial user objective
├── config/
│   ├── config.yaml                    # Pipeline parameters & container settings
│   └── samples.tsv                    # Sample sheet (sample, run, fastq1, fastq2)
├── workflow/
│   ├── Snakefile                      # Main workflow entrypoint
│   ├── rules/
│   │   ├── 00_common.smk              # Sample parsing & reference handling
│   │   ├── 01_qc.smk                  # FastQC & pairsqc rules
│   │   ├── 02_alignment.smk           # 4DN BWA-MEM alignment rule
│   │   ├── 03_pairtools.smk           # Parse, sort, dedup, merge, filter rules
│   │   ├── 04_matrices.smk            # Cooler cload, zoomify, juicer pre, add-norm rules
│   │   ├── 05_downstream.smk          # cooltools expected-cis, insulation, eigs-cis rules
│   │   └── 06_provenance.smk          # Provenance & metadata tracking rule
│   └── profiles/
│       └── slurm/
│           └── config.yaml            # Snakemake 9 Slurm profile
├── scripts/
│   └── preflight_check.sh             # Validates Apptainer, container, and wrappers
└── test/
    ├── test_samples.tsv               # Official smoke-test sample sheet
    ├── test_config.yaml               # Smoke-test configuration
    └── run_test.sh                    # End-to-end smoke-test runner
```

---

## Quickstart Guide

### 1. Load Modules on Cluster
```bash
module load apptainer/1.3.6
module load python/3.14.7
```

### 2. Run Container Preflight Check
Before submitting jobs, verify that Apptainer can execute the container and all official 4DN wrappers:
```bash
./scripts/preflight_check.sh
```

### 3. Run the Smoke Test
Verify DAG construction and Snakemake syntax using the official 4DN sample dataset:
```bash
./test/run_test.sh
```

---

## Production Execution

All direct copy-pasteable commands are available in **[`COMMANDS.txt`](file:///Users/carbajala/Antigravity_working_space/microC_pipeline/COMMANDS.txt)**.

### 1. Prepare Sample Sheet (`config/samples.tsv`)
Edit `config/samples.tsv` with your sample IDs, sequencing runs/lanes, and paired FASTQ paths:
```tsv
sample	run	fastq1	fastq2
Sample1	run1	/data/Sample1_L001_R1.fastq.gz	/data/Sample1_L001_R2.fastq.gz
Sample1	run2	/data/Sample1_L002_R1.fastq.gz	/data/Sample1_L002_R2.fastq.gz
Sample2	run1	/data/Sample2_R1.fastq.gz	/data/Sample2_R2.fastq.gz
```
> **Note**: Multiple sequencing runs (lanes) for the same sample are automatically merged before duplicate removal, exactly adhering to 4DN standards.

### 2. Configure Reference Genome (`config/config.yaml`)
```yaml
reference:
  assembly: "mm39"
  fasta: ""
  bwa_index: "/archive/cbds/idx/ensembl/Mus_musculus/GRCm39/sequence_only/bwa/0.7.18"
  chrom_sizes: ""
```

### 3. Launch Pipeline on Slurm HPC

```bash
# Dry run (inspect all shell commands and save output to dry_run.txt):
snakemake -s workflow/Snakefile --configfile config/config.yaml --profile workflow/profiles/slurm -n -p four_dn | tee dry_run.txt

# Run official 4DN pipeline (FastQC, alignment, pairs, .mcool, pairsqc)
snakemake -s workflow/Snakefile --configfile config/config.yaml --profile workflow/profiles/slurm four_dn

# Run 4DN pipeline + downstream cooltools analyses (P(s), TAD boundaries, compartments)
snakemake -s workflow/Snakefile --configfile config/config.yaml --profile workflow/profiles/slurm all
```

---

## Output Deliverables

All deliverables are organized in the `results/` directory:

| Directory | Output Files | Description |
| :--- | :--- | :--- |
| `results/pairs/` | `<sample>.ff.pairs.gz`, `.px2` | Valid, deduplicated 4DN contact pairs with Pairix index |
| `results/mcool/` | `<sample>.mcool` | Multiresolution matrix (1kb to 10Mb) with **ICE** and **Juicer** norm vectors |
| `results/hic/` | `<sample>.hic` | Juicebox multiresolution matrix with VC, VC_SQRT, KR normalization |
| `results/qc/fastqc/` | `*_fastqc.html`, `*_fastqc.zip` | FastQC reports on raw reads |
| `results/qc/pairsqc/` | `<sample>_report.zip` | Comprehensive 4DN contact statistics and distance decay report |
| `results/downstream/`| `*_expected_1000bp.tsv`<br>`*_insulation_10000bp.tsv`<br>`*_eigs_50000bp.cis.vecs.tsv` | $P(s)$ scaling curve, TAD/insulation boundaries, A/B compartments |
| `results/` | `provenance.json` | Complete metadata recording versions, checksums, and parameters |

---

## Provenance & Reproducibility
Every pipeline execution generates a machine-readable `results/provenance.json` containing:
- 4DN pipeline release version (`v43`)
- Container image URI and digest (`duplexa/4dn-hic:v43`)
- Reference FASTA SHA256 checksum
- Exact resolutions, parameters, and Snakemake execution timestamp
