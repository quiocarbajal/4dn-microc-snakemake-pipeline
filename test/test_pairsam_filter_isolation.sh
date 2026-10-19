#!/bin/bash
# =============================================================================
# Isolation test for rule pairsam_filter (shadow: "minimal")
#
# run-pairsam-filter.sh writes fixed-name temp.gz / temp1.gz into its working
# directory, so concurrent pairsam_filter jobs could swap samples. This test:
#   1. slices the first READ_PAIRS_PER_SLICE read pairs from two production
#      BAMs whose read IDs differ (PL rep1 = SRR31683399, LE rep1 = SRR31683386);
#   2. runs pairsam_parse_sort -> pairsam_markasdup -> pairsam_filter with the
#      fixed rule in a separate test directory, starting both filter jobs at
#      the same moment, while a watcher logs where temp.gz / temp1.gz appear;
#   3. checks header lineage and whole-file read-ID prefixes of each output;
#   4. re-runs the original wrapper call (no shadow, one sample at a time) on
#      the same marked pairsam and checks that the pairs and lossless BAM
#      records are identical to the fixed rule's output.
# The production run directory is only read, never written.
#
# Usage (inside an interactive allocation with at least 4 cores):
#   srun --cpus-per-task=4 --mem=24G --time=02:00:00 --pty bash
#   module load apptainer/1.3.6
#   module load python/3.14.7
#   bash test_pairsam_filter_isolation.sh
#
# Expects the fixed 03_pairtools.smk next to this script (override with
# FIXED_PAIRTOOLS_RULE_FILE). Writes everything into ./isolation_test_run
# (override with TEST_DIRECTORY), which must not exist yet.
# =============================================================================

# No pipefail: several pipelines stop reading early (head/awk exit) on purpose
set -eu

script_directory="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
production_run_directory="${PRODUCTION_RUN_DIRECTORY:-/s/pezzar-lab/Agustin/chen_2026_microC_pratto/fourDN_snake_pipeline}"
fixed_pairtools_rule_file="${FIXED_PAIRTOOLS_RULE_FILE:-$script_directory/03_pairtools.smk}"
test_directory="${TEST_DIRECTORY:-$script_directory/isolation_test_run}"
read_pairs_per_slice="${READ_PAIRS_PER_SLICE:-1000000}"
n_cores="${SLURM_CPUS_PER_TASK:-4}"

apptainer_cache_directory="/hpc-prj/pezza/containers/microc_4DN_pipeline_apptainer_cache"
container_image="$apptainer_cache_directory/dfae1079fd3684a3668ad68114f19a33.simg"
apptainer_bind_arguments="--bind /s/pezzar-lab,/hpc-prj/pezza,/Volumes/shared-refs:/Volumes/shared-refs:ro,/archive/cbds:/archive/cbds:ro"

test_samples=(PL_B6xCAST_rep1 LE_B6xCAST_rep1)
declare -A expected_read_id_prefix=([PL_B6xCAST_rep1]=SRR31683399 [LE_B6xCAST_rep1]=SRR31683386)

failed_check_count=0
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

header_lineage_tokens() {
    zcat "$1" | awk '/^#/ {print; next} {exit}' \
        | grep -a -o -E 'pairsam_marked/[A-Za-z0-9_]+\.marked|[A-Za-z0-9_]+_BioRep[0-9]+_R?[12]\.fastq\.gz' \
        | sort -u || true
}

pairs_body_md5() {
    zcat "$1" | grep -a -v '^#' | md5sum | cut -d' ' -f1
}

echo "== 0. Preconditions"
command -v apptainer >/dev/null || { echo "apptainer not found: run 'module load apptainer/1.3.6'"; exit 1; }
command -v snakemake >/dev/null || { echo "snakemake not found: run 'module load python/3.14.7'"; exit 1; }
[[ -f "$container_image" ]] || { echo "Container image not found: $container_image"; exit 1; }
[[ -f "$fixed_pairtools_rule_file" ]] || { echo "Fixed rule file not found: $fixed_pairtools_rule_file"; exit 1; }
[[ -e "$test_directory" ]] && { echo "Test directory already exists, choose another or remove it: $test_directory"; exit 1; }
(( n_cores >= 4 )) || { echo "Need at least 4 cores so both filter jobs (2 threads each) run at once; have $n_cores"; exit 1; }
production_chrom_sizes="$production_run_directory/results/reference/chrom.sizes"
[[ -f "$production_chrom_sizes" ]] || { echo "Chromosome sizes not found: $production_chrom_sizes"; exit 1; }
for sample in "${test_samples[@]}"; do
    [[ -f "$production_run_directory/results/bam/${sample}_1.bam" ]] || { echo "Production BAM missing for $sample"; exit 1; }
done
snakemake --version

echo "== 1. Build the test workspace from the production workflow plus the fixed rule"
mkdir -p "$test_directory/config" "$test_directory/results/bam"
cp -r "$production_run_directory/workflow" "$test_directory/workflow"
cp "$production_run_directory/config/config.yaml" "$test_directory/config/config.yaml"
echo "-- diff between production 03_pairtools.smk and the fixed one (expect only the shadow lines):"
diff "$production_run_directory/workflow/rules/03_pairtools.smk" "$fixed_pairtools_rule_file" || true
cp "$fixed_pairtools_rule_file" "$test_directory/workflow/rules/03_pairtools.smk"

fastq_directory="$(cd "$production_run_directory/../FASTQs/SRA_standard_download" && pwd)"
{
    printf 'condition\tsample\trun\tfastq1\tfastq2\n'
    for sample in "${test_samples[@]}"; do
        fastq_stem="${sample/_rep/_BioRep}"
        printf '%s\t%s\t1\t%s\t%s\n' "${sample%_rep*}" "$sample" \
            "$fastq_directory/${fastq_stem}_1.fastq.gz" "$fastq_directory/${fastq_stem}_2.fastq.gz"
    done
} > "$test_directory/config/test_samples.tsv"

cat > "$test_directory/config/test_config.yaml" <<EOF
samples: "config/test_samples.tsv"
outdir: "results"
save_bam: true
reference:
  chrom_sizes: "$production_chrom_sizes"
EOF

echo "== 2. Slice the first $read_pairs_per_slice read pairs of each production BAM"
for sample in "${test_samples[@]}"; do
    source_bam="$production_run_directory/results/bam/${sample}_1.bam"
    slice_bam="$test_directory/results/bam/${sample}_1.bam"
    # bwa -SP5M output is grouped by read name, so cutting at a name boundary keeps every pair whole
    run_in_container samtools view -h "$source_bam" \
        | awk -v maximum_read_names="$read_pairs_per_slice" \
            '/^@/ {print; next} $1 != previous_read_name {read_name_count++; previous_read_name = $1} read_name_count > maximum_read_names {exit} {print}' \
        | run_in_container samtools view -b -o "$slice_bam" -
    slice_prefixes="$(run_in_container samtools view "$slice_bam" | awk '{split($1, read_id_parts, "."); print read_id_parts[1]}' | sort -u | tr '\n' ' ')"
    echo "$sample slice read-ID prefixes: $slice_prefixes"
    [[ "$slice_prefixes" == "${expected_read_id_prefix[$sample]} " ]] \
        || { echo "Slice of $sample does not carry only ${expected_read_id_prefix[$sample]}; stopping"; exit 1; }
done

cd "$test_directory"
snakemake_common_arguments=(
    -s workflow/Snakefile --configfile config/test_config.yaml
    --cores "$n_cores" --use-apptainer
    --apptainer-prefix "$apptainer_cache_directory" --apptainer-args "$apptainer_bind_arguments"
    --allowed-rules pairsam_parse_sort pairsam_markasdup pairsam_filter
    --notemp --printshellcmds
)

echo "== 3a. Snakemake: parse/sort and duplicate marking for both samples"
snakemake "${snakemake_common_arguments[@]}" \
    $(for sample in "${test_samples[@]}"; do echo "results/pairsam_marked/$sample.marked.sam.pairs.gz"; done)

echo "== 3b. Snakemake: both pairsam_filter jobs at once, with a temp-file watcher"
temp_file_watch_log="$test_directory/temp_file_watch.log"
(
    while true; do
        root_temp_files="$(ls temp.gz temp1.gz 2>/dev/null | tr '\n' ',' || true)"
        shadow_temp_files="$(find .snakemake/shadow -name 'temp*.gz' 2>/dev/null | tr '\n' ',' || true)"
        echo "$(date +%T) root:$root_temp_files shadow:$shadow_temp_files" >> "$temp_file_watch_log"
        sleep 2
    done
) &
watcher_process_id=$!
trap 'kill "$watcher_process_id" 2>/dev/null || true' EXIT
snakemake "${snakemake_common_arguments[@]}" \
    $(for sample in "${test_samples[@]}"; do echo "results/pairs_filtered/$sample.dedup.pairs.gz"; done)
kill "$watcher_process_id" 2>/dev/null || true

echo "== 4. Regression: original wrapper call (no shadow), one sample at a time"
mkdir -p regression_original_wrapper
for sample in "${test_samples[@]}"; do
    (cd regression_original_wrapper && run_in_container /usr/local/bin/run-pairsam-filter.sh \
        "../results/pairsam_marked/$sample.marked.sam.pairs.gz" "$sample" "$production_chrom_sizes" \
        > "$sample.filter.log" 2>&1)
done

echo
echo "================================ RESULTS ================================"
watch_lines_with_root_temp_files="$(awk '$2 != "root:"' "$temp_file_watch_log" | wc -l)"
max_simultaneous_shadow_directories="$(awk '{
    sub(/^shadow:/, "", $3); path_count = split($3, temp_paths, ","); split("", shadow_directories); directory_count = 0
    for (i = 1; i <= path_count; i++) if (temp_paths[i] != "") {
        shadow_directory = temp_paths[i]; sub(/\/[^\/]*$/, "", shadow_directory)
        if (!(shadow_directory in shadow_directories)) { shadow_directories[shadow_directory] = 1; directory_count++ }
    }
    if (directory_count > maximum) maximum = directory_count
} END { print maximum + 0 }' "$temp_file_watch_log")"
echo "Watcher samples: $(wc -l < "$temp_file_watch_log"); most shadow directories holding temp files at once: $max_simultaneous_shadow_directories"

report_check "no temp.gz/temp1.gz in the test directory root at any watcher sample" \
    "$([[ "$watch_lines_with_root_temp_files" -eq 0 ]] && echo yes || echo no)"
report_check "no temp.gz/temp1.gz in the test directory root after the run" \
    "$([[ ! -e temp.gz && ! -e temp1.gz ]] && echo yes || echo no)"
report_check "both filter jobs overlapped (temp files seen in 2 shadow directories at once)" \
    "$([[ "$max_simultaneous_shadow_directories" -ge 2 ]] && echo yes || echo no)"
report_check "positive control: original wrapper does write temp.gz into its working directory" \
    "$([[ -e regression_original_wrapper/temp.gz ]] && echo yes || echo no)"

for sample in "${test_samples[@]}"; do
    fixed_dedup_pairs="results/pairs_filtered/$sample.dedup.pairs.gz"
    original_dedup_pairs="regression_original_wrapper/$sample.dedup.pairs.gz"
    fastq_stem="${sample/_rep/_BioRep}"
    echo "-- $sample"

    lineage_tokens="$(header_lineage_tokens "$fixed_dedup_pairs" | tr '\n' ' ')"
    echo "   header lineage tokens: $lineage_tokens"
    foreign_lineage_tokens="$(header_lineage_tokens "$fixed_dedup_pairs" | grep -v -e "pairsam_marked/$sample\.marked" -e "^${fastq_stem}_" || true)"
    report_check "$sample: header lineage names only $sample" \
        "$([[ -n "$lineage_tokens" && -z "$foreign_lineage_tokens" ]] && echo yes || echo no)"

    read_id_prefix_counts="$(zcat "$fixed_dedup_pairs" | awk '!/^#/ {split($1, read_id_parts, "."); prefix_counts[read_id_parts[1]]++} END {for (prefix in prefix_counts) print prefix_counts[prefix], prefix}')"
    echo "   read-ID prefix counts (whole file): $read_id_prefix_counts"
    report_check "$sample: every pair carries ${expected_read_id_prefix[$sample]}" \
        "$([[ -n "$read_id_prefix_counts" && "$(echo "$read_id_prefix_counts" | awk '{print $2}')" == "${expected_read_id_prefix[$sample]}" ]] && echo yes || echo no)"

    report_check "$sample: dedup pairs identical to original wrapper output (header excluded)" \
        "$([[ "$(pairs_body_md5 "$fixed_dedup_pairs")" == "$(pairs_body_md5 "$original_dedup_pairs")" ]] && echo yes || echo no)"
    report_check "$sample: .px2 index present" \
        "$([[ -s "$fixed_dedup_pairs.px2" ]] && echo yes || echo no)"
    report_check "$sample: lossless BAM records identical to original wrapper output" \
        "$([[ "$(run_in_container samtools view "results/bam_lossless/$sample.lossless.bam" | md5sum)" == "$(run_in_container samtools view "regression_original_wrapper/$sample.lossless.bam" | md5sum)" ]] && echo yes || echo no)"
    echo "   pairtools stats, fixed vs original (first 12 lines, expect no difference):"
    diff <(run_in_container pairtools stats "$fixed_dedup_pairs" | head -n 12) \
         <(run_in_container pairtools stats "$original_dedup_pairs" | head -n 12) \
        && echo "   (identical)" || true
done

echo "========================================================================="
if [[ "$failed_check_count" -eq 0 ]]; then
    echo "ALL CHECKS PASSED"
else
    echo "$failed_check_count CHECK(S) FAILED"
    exit 1
fi
