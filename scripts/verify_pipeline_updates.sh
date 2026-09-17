#!/usr/bin/env bash
# =============================================================================
# Pipeline Updates Verification Script
# Checks whether critical bug fixes and QC toggles are present in the codebase.
# =============================================================================

set -u

ERRORS=0

echo "=== Checking Micro-C Pipeline Updates ==="

# 1. Check Snakefile for conditional pairsqc targets
if grep -q 'run_pairsqc' workflow/Snakefile 2>/dev/null; then
    echo "  [PASS] workflow/Snakefile: pairsqc targets are conditional on qc.run_pairsqc"
else
    echo "  [FAIL] workflow/Snakefile: missing qc.run_pairsqc condition (pairsqc targets will always run)"
    ERRORS=$((ERRORS + 1))
fi

# 2. Check 01_qc.smk does not use the buggy upstream run-pairsqc-merge.sh
if grep -q 'run-pairsqc-merge.sh' workflow/rules/01_qc.smk 2>/dev/null; then
    echo "  [FAIL] workflow/rules/01_qc.smk: still calling buggy upstream run-pairsqc-merge.sh"
    ERRORS=$((ERRORS + 1))
else
    echo "  [PASS] workflow/rules/01_qc.smk: upstream run-pairsqc-merge.sh bypassed"
fi

# 3. Check 01_qc.smk has the self-move collision guard
if grep -q 'output.report_zip.*!=' workflow/rules/01_qc.smk 2>/dev/null || grep -q '!=.*output.report_zip' workflow/rules/01_qc.smk 2>/dev/null; then
    echo "  [PASS] workflow/rules/01_qc.smk: mv self-move collision guard is present"
else
    echo "  [FAIL] workflow/rules/01_qc.smk: missing mv self-move guard in pairsqc_merge"
    ERRORS=$((ERRORS + 1))
fi

# 4. Check config.yaml QC status
echo ""
echo "=== Current QC Configuration (config/config.yaml) ==="
if [ -f "config/config.yaml" ]; then
    grep -E '^\s*run_(pairsqc|fastqc):' config/config.yaml || echo "  [WARN] No run_pairsqc or run_fastqc keys found in config/config.yaml"
else
    echo "  [FAIL] config/config.yaml not found"
    ERRORS=$((ERRORS + 1))
fi

echo ""
if [ $ERRORS -eq 0 ]; then
    echo "Result: ALL NECESSARY FIXES ARE PRESENT."
else
    echo "Result: $ERRORS ISSUE(S) DETECTED. Please sync the latest pipeline files."
fi
exit $ERRORS
