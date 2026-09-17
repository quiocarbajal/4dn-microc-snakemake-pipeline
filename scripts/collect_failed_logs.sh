#!/usr/bin/env bash
# =============================================================================
# Script to collect all logs from failed rules (run: 2026-09-11T213125.136235)
# Combines application rule logs and SLURM stdout/stderr into a single file.
# =============================================================================

set -u

OUTPUT_FILE="${1:-failed_rules_logs.txt}"

# Explicit paths extracted from 2026-09-11T213125.136235.snakemake.log
LOG_FILES=(
    # --- Rule: pairsqc_merge (sample: PL_B6xCAST | job: 96 | slurm: 6206868) ---
    "results/logs/pairsqc/merge_PL_B6xCAST.log"
    "/s/pezzar-lab/chen_2026_microC_pratto/fourDN_snake_pipeline/logs/slurm/rule_pairsqc_merge/PL_B6xCAST/6206868.log"

    # --- Rule: pairsqc_merge (sample: ZY_B6xCAST | job: 97 | slurm: 6206869) ---
    "results/logs/pairsqc/merge_ZY_B6xCAST.log"
    "/s/pezzar-lab/chen_2026_microC_pratto/fourDN_snake_pipeline/logs/slurm/rule_pairsqc_merge/ZY_B6xCAST/6206869.log"

    # --- Rule: pairsam_filter (sample: LE_B6xCAST_rep2 | job: 13 | slurm: 6206864) ---
    "results/logs/pairtools/filter_LE_B6xCAST_rep2.log"
    "/s/pezzar-lab/chen_2026_microC_pratto/fourDN_snake_pipeline/logs/slurm/rule_pairsam_filter/LE_B6xCAST_rep2/6206864.log"

    # --- Rule: pairsam_filter (sample: PA_B6xCAST_rep3 | job: 28 | slurm: 6206866) ---
    "results/logs/pairtools/filter_PA_B6xCAST_rep3.log"
    "/s/pezzar-lab/chen_2026_microC_pratto/fourDN_snake_pipeline/logs/slurm/rule_pairsam_filter/PA_B6xCAST_rep3/6206866.log"
)

# Clear or initialize output file
> "${OUTPUT_FILE}"

echo "Collecting ${#LOG_FILES[@]} log files into: ${OUTPUT_FILE}"

for log_path in "${LOG_FILES[@]}"; do
    printf "\n\n" >> "${OUTPUT_FILE}"
    printf "================================================================================\n" >> "${OUTPUT_FILE}"
    printf "=== LOG FILE: %s\n" "${log_path}" >> "${OUTPUT_FILE}"
    printf "================================================================================\n\n" >> "${OUTPUT_FILE}"

    if [ -f "${log_path}" ]; then
        echo "  [FOUND] ${log_path}"
        cat "${log_path}" >> "${OUTPUT_FILE}" 2>&1
    else
        echo "  [MISSING] ${log_path}"
        printf "[ERROR: File not found or unreadable: %s]\n" "${log_path}" >> "${OUTPUT_FILE}"
    fi
done

printf "\nDone! All collected logs written to: %s\n" "${OUTPUT_FILE}"
