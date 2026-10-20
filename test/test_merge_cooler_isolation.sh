#!/bin/bash
# =============================================================================
# Isolation test for rules merge_bio_reps, cooler_cload and cooler_cload_merged
# (shadow: "minimal")
#
# run-merge-pairs.sh creates fixed-name FIFOs pp.1, pp.2, ... and run-cooler.sh
# copies the chrom sizes to a fixed-name ./tmpchrsize, both in the working
# directory, so concurrent jobs could share them. This test:
#   1. reuses the two dedup pairs files made by test_pairsam_filter_isolation.sh
#      (PL rep1 = SRR31683399, LE rep1 = SRR31683386) and gives each condition
#      two replicates (rep2 is a byte copy of rep1), so the merge takes its
#      FIFO path;
#   2. runs merge_bio_reps for both conditions at the same moment, then
#      cooler_cload for the four samples and cooler_cload_merged for the two
#      conditions, with the fixed rules, while a watcher logs where pp.* and
#      tmpchrsize appear;
#   3. checks that each merged file holds only its own condition's read IDs
#      and exactly twice the pairs of one replicate;
#   4. re-runs the original wrapper calls (no shadow, one job at a time) and
#      checks that merged pairs and every cooler's chroms, bins and pixels are
#      identical to the fixed rules' output.
# The production run directory and the first test's directory are only read.
#
# Usage (inside an interactive allocation with at least 8 cores):
#   srun --cpus-per-task=8 --mem=32G --time=02:00:00 --pty bash
#   module load apptainer/1.3.6
#   module load python/3.14.7
#   bash test_merge_cooler_isolation.sh
#
# Expects next to this script: the fixed 03_pairtools.smk and 04_matrices.smk,
# and ./isolation_test_run from the first test (override with
# FIXED_PAIRTOOLS_RULE_FILE, FIXED_MATRICES_RULE_FILE, FIRST_TEST_DIRECTORY).
# Writes everything into ./merge_cooler_isolation_test_run (override with
# TEST_DIRECTORY), which must not exist yet.
# =============================================================================

# No pipefail: several pipelines stop reading early (head/awk exit) on purpose
set -eu

script_directory="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
production_run_directory="${PRODUCTION_RUN_DIRECTORY:-/s/pezzar-lab/Agustin/chen_2026_microC_pratto/fourDN_snake_pipeline}"
first_test_directory="${FIRST_TEST_DIRECTORY:-$script_directory/isolation_test_run}"
fixed_pairtools_rule_file="${FIXED_PAIRTOOLS_RULE_FILE:-$script_directory/03_pairtools.smk}"
fixed_matrices_rule_file="${FIXED_MATRICES_RULE_FILE:-$script_directory/04_matrices.smk}"
test_directory="${TEST_DIRECTORY:-$script_directory/merge_cooler_isolation_test_run}"
n_cores="${SLURM_CPUS_PER_TASK:-8}"

apptainer_cache_directory="/hpc-prj/pezza/containers/microc_4DN_pipeline_apptainer_cache"
container_image="$apptainer_cache_directory/dfae1079fd3684a3668ad68114f19a33.simg"
apptainer_bind_arguments="--bind /s/pezzar-lab,/hpc-prj/pezza,/Volumes/shared-refs:/Volumes/shared-refs:ro,/archive/cbds:/archive/cbds:ro"

test_conditions=(PL_B6xCAST LE_B6xCAST)
test_samples=(PL_B6xCAST_rep1 PL_B6xCAST_rep2 LE_B6xCAST_rep1 LE_B6xCAST_rep2)
declare -A expected_read_id_prefix=([PL_B6xCAST]=SRR31683399 [LE_B6xCAST]=SRR31683386)
# Test-only bin size: at the pipeline's 1 kb (about 2.7M mouse bins) one cload takes about 28 min;
# the shared-file hazard does not depend on bin size
cooler_bin_size=1000000
cooler_threads=4
cooler_max_split=2

failed_check_count=0
watcher_process_id=""
report_check() {
    local check_description="$1" check_passed="$2"
    if [[ "$check_passed" == "yes" ]]; then
        echo "PASS  $check_description"
    else
        echo "FAIL  $check_description"
        failed_check_count=$((failed_check_count + 1))
    fi
}

run_in_container() {
    apptainer exec $apptainer_bind_arguments "$container_image" "$@"
}

pairs_body_md5() {
    zcat "$1" | grep -a -v '^#' | md5sum | cut -d' ' -f1
}

pair_count() {
    zcat "$1" | grep -a -c -v '^#' || true
}

cooler_table_md5() {
    run_in_container cooler dump --table "$2" "$1" | md5sum | cut -d' ' -f1
}

# Logs, every 0.5 s, any pp.* or tmpchrsize in the test directory root and in shadow directories
start_shared_file_watcher() {
    local watch_log="$1"
    (
        while true; do
            root_shared_files="$(ls -d pp.* tmpchrsize 2>/dev/null | tr '\n' ',' || true)"
            shadow_shared_files="$(find .snakemake/shadow \( -name 'pp.*' -o -name tmpchrsize \) 2>/dev/null | tr '\n' ',' || true)"
            echo "$(date +%T.%N) root:$root_shared_files shadow:$shadow_shared_files" >> "$watch_log"
            sleep 0.5
        done
    ) &
    watcher_process_id=$!
}

stop_shared_file_watcher() {
    kill "$watcher_process_id" 2>/dev/null || true
    wait "$watcher_process_id" 2>/dev/null || true
}

watch_lines_with_root_files() {
    awk '$2 != "root:"' "$1" | wc -l
}

max_simultaneous_shadow_directories() {
    awk '{
        sub(/^shadow:/, "", $3); path_count = split($3, shared_paths, ","); split("", shadow_directories); directory_count = 0
        for (i = 1; i <= path_count; i++) if (shared_paths[i] != "") {
            shadow_directory = shared_paths[i]; sub(/\/[^\/]*$/, "", shadow_directory)
            if (!(shadow_directory in shadow_directories)) { shadow_directories[shadow_directory] = 1; directory_count++ }
        }
        if (directory_count > maximum) maximum = directory_count
    } END { print maximum + 0 }' "$1"
}

echo "== 0. Preconditions"
command -v apptainer >/dev/null || { echo "apptainer not found: run 'module load apptainer/1.3.6'"; exit 1; }
command -v snakemake >/dev/null || { echo "snakemake not found: run 'module load python/3.14.7'"; exit 1; }
[[ -f "$container_image" ]] || { echo "Container image not found: $container_image"; exit 1; }
[[ -f "$fixed_pairtools_rule_file" ]] || { echo "Fixed rule file not found: $fixed_pairtools_rule_file"; exit 1; }
[[ -f "$fixed_matrices_rule_file" ]] || { echo "Fixed rule file not found: $fixed_matrices_rule_file"; exit 1; }
[[ -e "$test_directory" ]] && { echo "Test directory already exists, choose another or remove it: $test_directory"; exit 1; }
(( n_cores >= 8 )) || { echo "Need at least 8 cores so two cooler jobs (4 threads each) run at once; have $n_cores"; exit 1; }
production_chrom_sizes="$production_run_directory/results/reference/chrom.sizes"
[[ -f "$production_chrom_sizes" ]] || { echo "Chromosome sizes not found: $production_chrom_sizes"; exit 1; }
for condition in "${test_conditions[@]}"; do
    first_test_pairs="$first_test_directory/results/pairs_filtered/${condition}_rep1.dedup.pairs.gz"
    [[ -f "$first_test_pairs" && -f "$first_test_pairs.px2" ]] || { echo "First-test pairs or index missing: $first_test_pairs"; exit 1; }
done
snakemake --version

echo "== 1. Build the test workspace from the production workflow plus the fixed rules"
mkdir -p "$test_directory/config" "$test_directory/results/pairs"
cp -r "$production_run_directory/workflow" "$test_directory/workflow"
cp "$production_run_directory/config/config.yaml" "$test_directory/config/config.yaml"
echo "-- diff between production 03_pairtools.smk and the fixed one (expect only the two shadow blocks):"
diff "$production_run_directory/workflow/rules/03_pairtools.smk" "$fixed_pairtools_rule_file" || true
echo "-- diff between production 04_matrices.smk and the fixed one (expect the two shadow blocks, plus any other local edits):"
diff "$production_run_directory/workflow/rules/04_matrices.smk" "$fixed_matrices_rule_file" || true
cp "$fixed_pairtools_rule_file" "$test_directory/workflow/rules/03_pairtools.smk"
cp "$fixed_matrices_rule_file" "$test_directory/workflow/rules/04_matrices.smk"

fastq_directory="$(cd "$production_run_directory/../FASTQs/SRA_standard_download" && pwd)"
{
    printf 'condition\tsample\trun\tfastq1\tfastq2\n'
    for sample in "${test_samples[@]}"; do
        condition="${sample%_rep*}"
        # FASTQ columns are required by the sample sheet but unused: no rule that reads them is allowed to run
        printf '%s\t%s\t1\t%s\t%s\n' "$condition" "$sample" \
            "$fastq_directory/${condition}_BioRep1_1.fastq.gz" "$fastq_directory/${condition}_BioRep1_2.fastq.gz"
    done
} > "$test_directory/config/test_samples.tsv"

cat > "$test_directory/config/test_config.yaml" <<EOF
samples: "config/test_samples.tsv"
outdir: "results"
save_bam: true
reference:
  chrom_sizes: "$production_chrom_sizes"
matrix:
  base_resolution: $cooler_bin_size
EOF

echo "== 2. Stage two replicates per condition from the first test's dedup pairs"
for sample in "${test_samples[@]}"; do
    condition="${sample%_rep*}"
    first_test_pairs="$first_test_directory/results/pairs_filtered/${condition}_rep1.dedup.pairs.gz"
    cp "$first_test_pairs" "$test_directory/results/pairs/$sample.ff.pairs.gz"
    cp "$first_test_pairs.px2" "$test_directory/results/pairs/$sample.ff.pairs.gz.px2"
done

cd "$test_directory"
snakemake_common_arguments=(
    -s workflow/Snakefile --configfile config/test_config.yaml
    --cores "$n_cores" --use-apptainer
    --apptainer-prefix "$apptainer_cache_directory" --apptainer-args "$apptainer_bind_arguments"
    --allowed-rules merge_bio_reps cooler_cload cooler_cload_merged
    --printshellcmds
)
trap 'stop_shared_file_watcher' EXIT

echo "== 3a. Snakemake: merge_bio_reps for both conditions at once, with a FIFO watcher"
merge_watch_log="$test_directory/merge_shared_file_watch.log"
start_shared_file_watcher "$merge_watch_log"
snakemake "${snakemake_common_arguments[@]}" \
    $(for condition in "${test_conditions[@]}"; do echo "results/pairs_merged/$condition.merged.ff.pairs.gz"; done)
stop_shared_file_watcher

echo "== 3b. Snakemake: cooler_cload (4 samples) and cooler_cload_merged (2 conditions), with a tmpchrsize watcher"
cooler_watch_log="$test_directory/cooler_shared_file_watch.log"
start_shared_file_watcher "$cooler_watch_log"
snakemake "${snakemake_common_arguments[@]}" \
    $(for sample in "${test_samples[@]}"; do echo "results/cool/$sample.cool"; done) \
    $(for condition in "${test_conditions[@]}"; do echo "results/cool_merged/$condition.cool"; done)
stop_shared_file_watcher

echo "== 4. Regression: original wrapper calls (no shadow), one job at a time"
mkdir -p regression_original_wrappers
for condition in "${test_conditions[@]}"; do
    (cd regression_original_wrappers && run_in_container /usr/local/bin/run-merge-pairs.sh \
        "$condition.merged.ff" "../results/pairs/${condition}_rep1.ff.pairs.gz" "../results/pairs/${condition}_rep2.ff.pairs.gz" \
        > "merge_$condition.log" 2>&1)
done
for sample in "${test_samples[@]}"; do
    (cd regression_original_wrappers && run_in_container /usr/local/bin/run-cooler.sh \
        "../results/pairs/$sample.ff.pairs.gz" "$production_chrom_sizes" "$cooler_bin_size" "$cooler_threads" "$sample" "$cooler_max_split" \
        > "cload_$sample.log" 2>&1)
done
for condition in "${test_conditions[@]}"; do
    (cd regression_original_wrappers && run_in_container /usr/local/bin/run-cooler.sh \
        "../results/pairs_merged/$condition.merged.ff.pairs.gz" "$production_chrom_sizes" "$cooler_bin_size" "$cooler_threads" "$condition.merged" "$cooler_max_split" \
        > "cload_merged_$condition.log" 2>&1)
done

echo
echo "================================ RESULTS ================================"
merge_max_shadow_directories="$(max_simultaneous_shadow_directories "$merge_watch_log")"
cooler_max_shadow_directories="$(max_simultaneous_shadow_directories "$cooler_watch_log")"
echo "Merge watcher samples: $(wc -l < "$merge_watch_log"); most shadow directories holding pp.* at once: $merge_max_shadow_directories"
echo "Cooler watcher samples: $(wc -l < "$cooler_watch_log"); most shadow directories holding tmpchrsize at once: $cooler_max_shadow_directories"

report_check "no pp.* in the test directory root at any merge watcher sample" \
    "$([[ "$(watch_lines_with_root_files "$merge_watch_log")" -eq 0 ]] && echo yes || echo no)"
report_check "no tmpchrsize in the test directory root at any cooler watcher sample" \
    "$([[ "$(watch_lines_with_root_files "$cooler_watch_log")" -eq 0 ]] && echo yes || echo no)"
report_check "no pp.* or tmpchrsize in the test directory root after the runs" \
    "$([[ -z "$(ls -d pp.* tmpchrsize 2>/dev/null || true)" ]] && echo yes || echo no)"
report_check "both merge jobs overlapped (pp.* seen in 2 shadow directories at once)" \
    "$([[ "$merge_max_shadow_directories" -ge 2 ]] && echo yes || echo no)"
report_check "cooler jobs overlapped (tmpchrsize seen in 2 or more shadow directories at once)" \
    "$([[ "$cooler_max_shadow_directories" -ge 2 ]] && echo yes || echo no)"
report_check "positive control: original run-cooler.sh does write tmpchrsize into its working directory" \
    "$([[ -e regression_original_wrappers/tmpchrsize ]] && echo yes || echo no)"

chrom_sizes_line_count="$(grep -c . "$production_chrom_sizes")"
for condition in "${test_conditions[@]}"; do
    fixed_merged_pairs="results/pairs_merged/$condition.merged.ff.pairs.gz"
    original_merged_pairs="regression_original_wrappers/$condition.merged.ff.pairs.gz"
    echo "-- $condition (merged)"

    read_id_prefix_counts="$(zcat "$fixed_merged_pairs" | awk '!/^#/ {split($1, read_id_parts, "."); prefix_counts[read_id_parts[1]]++} END {for (prefix in prefix_counts) print prefix_counts[prefix], prefix}')"
    echo "   read-ID prefix counts (whole file): $read_id_prefix_counts"
    report_check "$condition: every merged pair carries ${expected_read_id_prefix[$condition]}" \
        "$([[ -n "$read_id_prefix_counts" && "$(echo "$read_id_prefix_counts" | awk '{print $2}')" == "${expected_read_id_prefix[$condition]}" ]] && echo yes || echo no)"

    replicate_pair_count="$(pair_count "results/pairs/${condition}_rep1.ff.pairs.gz")"
    merged_pair_count="$(pair_count "$fixed_merged_pairs")"
    echo "   pairs per replicate: $replicate_pair_count; merged pairs: $merged_pair_count"
    report_check "$condition: merged pairs = 2 x pairs per replicate" \
        "$([[ "$merged_pair_count" -eq $((2 * replicate_pair_count)) ]] && echo yes || echo no)"
    report_check "$condition: merged pairs identical to original wrapper output (header excluded)" \
        "$([[ "$(pairs_body_md5 "$fixed_merged_pairs")" == "$(pairs_body_md5 "$original_merged_pairs")" ]] && echo yes || echo no)"
    report_check "$condition: merged .px2 index present" \
        "$([[ -s "$fixed_merged_pairs.px2" ]] && echo yes || echo no)"
done

compare_cooler_with_original() {
    local cooler_label="$1" fixed_cooler="$2" original_cooler="$3"
    echo "-- $cooler_label"
    local cooler_chrom_count
    cooler_chrom_count="$(run_in_container cooler dump --table chroms "$fixed_cooler" | wc -l)"
    report_check "$cooler_label: cooler has all $chrom_sizes_line_count chromosomes of chrom.sizes (has $cooler_chrom_count)" \
        "$([[ "$cooler_chrom_count" -eq "$chrom_sizes_line_count" ]] && echo yes || echo no)"
    local cooler_table
    for cooler_table in chroms bins pixels; do
        report_check "$cooler_label: $cooler_table table identical to original wrapper output" \
            "$([[ "$(cooler_table_md5 "$fixed_cooler" "$cooler_table")" == "$(cooler_table_md5 "$original_cooler" "$cooler_table")" ]] && echo yes || echo no)"
    done
}
for sample in "${test_samples[@]}"; do
    compare_cooler_with_original "$sample (cool)" "results/cool/$sample.cool" "regression_original_wrappers/$sample.cool"
done
for condition in "${test_conditions[@]}"; do
    compare_cooler_with_original "$condition (merged cool)" "results/cool_merged/$condition.cool" "regression_original_wrappers/$condition.merged.cool"
done

echo "========================================================================="
if [[ "$failed_check_count" -eq 0 ]]; then
    echo "ALL CHECKS PASSED"
else
    echo "$failed_check_count CHECK(S) FAILED"
    exit 1
fi
