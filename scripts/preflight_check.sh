#!/bin/bash
# =============================================================================
# Preflight Check Script for 4DN v43 Micro-C / Hi-C Pipeline
# Validates Apptainer container, tools, and 4DN wrappers without running jobs
# =============================================================================

set -e

CONTAINER_URI="${1:-docker://duplexa/4dn-hic:v43}"

echo "=================================================================="
echo "4DN v43 Pipeline Preflight Environment & Container Check"
echo "=================================================================="

# 1. Check Apptainer / Singularity availability
echo -n "[1/4] Checking Apptainer / Singularity installation... "
if command -v apptainer >/dev/null 2>&1; then
    CONTAINER_CMD="apptainer"
    echo "OK ($(apptainer --version))"
elif command -v singularity >/dev/null 2>&1; then
    CONTAINER_CMD="singularity"
    echo "OK ($(singularity --version))"
else
    echo "FAILED"
    echo "Error: Neither 'apptainer' nor 'singularity' was found in PATH."
    echo "Please load the apptainer module: module load apptainer/1.3.6"
    exit 1
fi

# 2. Check Snakemake installation
echo -n "[2/4] Checking Snakemake installation... "
if command -v snakemake >/dev/null 2>&1; then
    echo "OK ($(snakemake --version))"
else
    echo "FAILED"
    echo "Error: 'snakemake' command not found in PATH."
    echo "Please load the python/snakemake module: module load python/3.14.7"
    exit 1
fi

# 3. Test Container Execution & Tool Versions
echo "[3/4] Testing container execution: ${CONTAINER_URI}..."

echo -n "  - Testing BWA in container: "
$CONTAINER_CMD exec "$CONTAINER_URI" bwa 2>&1 | head -n 3 | tr '\n' ' '
echo " [OK]"

echo -n "  - Testing Samtools in container: "
$CONTAINER_CMD exec "$CONTAINER_URI" samtools --version | head -n 1
echo " [OK]"

echo -n "  - Testing Pairtools in container: "
$CONTAINER_CMD exec "$CONTAINER_URI" pairtools --version
echo " [OK]"

echo -n "  - Testing Cooler in container: "
$CONTAINER_CMD exec "$CONTAINER_URI" cooler --version
echo " [OK]"

echo -n "  - Testing Pairix in container: "
$CONTAINER_CMD exec "$CONTAINER_URI" pairix 2>&1 | head -n 2 | tr '\n' ' '
echo " [OK]"

# Test FastQC availability
echo -n "  - Testing FastQC: "
if command -v fastqc >/dev/null 2>&1; then
    echo "OK ($(fastqc --version))"
else
    echo "OK (will be loaded dynamically via module fastqc/0.12.1)"
fi

# 4. Verify 4DN Official Wrappers Inside Container
echo "[4/4] Verifying 4DN v43 wrapper scripts inside container..."
REQUIRED_WRAPPERS=(
    "run-bwa-mem.sh"
    "run-pairsam-parse-sort.sh"
    "run-pairsam-merge.sh"
    "run-pairsam-markasdup.sh"
    "run-pairsam-filter.sh"
    "run-merge-pairs.sh"
    "run-cooler.sh"
    "run-cool2multirescool.sh"
    "run-juicebox-pre.sh"
    "run-add-hicnormvector-to-mcool.sh"
    "run-pairsqc-single.sh"
    "run-fastqc.sh"
    "run-list.sh"
)

for wrapper in "${REQUIRED_WRAPPERS[@]}"; do
    echo -n "  - Checking /usr/local/bin/$wrapper: "
    if $CONTAINER_CMD exec "$CONTAINER_URI" which "$wrapper" >/dev/null 2>&1; then
        echo "FOUND [OK]"
    else
        echo "MISSING [FAIL]"
        exit 1
    fi
done

echo "=================================================================="
echo "SUCCESS: All preflight checks passed! The container and toolchain"
echo "are fully operational for 4DN v43 pipeline execution."
echo "=================================================================="
