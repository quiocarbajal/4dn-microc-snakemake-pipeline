#!/usr/bin/env python3
"""
4DN Micro-C Pipeline Timing Extractor
-------------------------------------
Extracts and summarizes execution wall-clock time and job statistics for each 
pipeline step leading up to the final merged .mcool file.

Queries:
1. SLURM job logs in logs/slurm/rule_* (extracts JobIDs and start/end timestamps)
2. SLURM accounting database (sacct) for exact cluster wall-clock runtimes
3. Tool logs in results/logs/
4. Snakemake metadata repository (.snakemake/metadata/*)
"""

import os
import sys
import re
import json
import glob
import argparse
import subprocess
from datetime import datetime, timedelta

STAGE_ORDER = [
    "bwa_mem",
    "fastqc_r1",
    "fastqc_r2",
    "pairsam_parse_sort",
    "pairsam_merge_seq_reps",
    "pairsam_markasdup",
    "pairsam_filter",
    "addfrag_or_pass",
    "cooler_cload",
    "cooler_zoomify",
    "pairsqc",
    "merge_bio_reps",
    "cooler_cload_merged",
    "cooler_zoomify_merged",
    "pairsqc_merge"
]

def format_duration(seconds):
    if seconds is None:
        return "N/A"
    td = timedelta(seconds=int(round(seconds)))
    total_sec = int(td.total_seconds())
    hours, remainder = divmod(total_sec, 3600)
    minutes, sec = divmod(remainder, 60)
    if hours > 0:
        return f"{hours:02d}h {minutes:02d}m {sec:02d}s"
    elif minutes > 0:
        return f"{minutes:02d}m {sec:02d}s"
    else:
        return f"{sec:02d}s"

def parse_slurm_logs(workdir, filter_keyword):
    """
    Scans logs/slurm/rule_<rulename>/<target>/<jobid>.log to extract:
    - Rule name
    - Wildcard target (sample / run / condition)
    - Slurm Job ID
    - Start timestamp (from header 'Date: YYYY-MM-DD HH:MM:SS')
    - Finish timestamp (from file mtime)
    - Exit status
    """
    slurm_dir = os.path.join(workdir, "logs", "slurm")
    records = []
    if not os.path.isdir(slurm_dir):
        return records

    date_re = re.compile(r"^\s*Date:\s*(\d{4}-\d{2}-\d{2}\s+\d{2}:\d{2}:\d{2})")
    error_re = re.compile(r"(WorkflowError|failed, SLURM status is: 'FAILED'|exited with non-zero)")

    for root, dirs, files in os.walk(slurm_dir):
        for fname in files:
            if not fname.endswith(".log"):
                continue
            
            jobid_str = os.path.splitext(fname)[0]
            if not jobid_str.isdigit():
                continue

            fpath = os.path.join(root, fname)
            rel_dir = os.path.relpath(root, slurm_dir)
            parts = rel_dir.split(os.sep)
            
            rule = parts[0].replace("rule_", "") if len(parts) >= 1 else "unknown"
            target = parts[1] if len(parts) >= 2 else ""

            if filter_keyword not in target and filter_keyword not in rel_dir:
                continue

            start_dt = None
            has_error = False
            try:
                with open(fpath, "r", encoding="utf-8", errors="ignore") as f:
                    # Scan first 30 lines for Date header
                    for i, line in enumerate(f):
                        if i < 30 and not start_dt:
                            m = date_re.match(line)
                            if m:
                                try:
                                    start_dt = datetime.strptime(m.group(1), "%Y-%m-%d %H:%M:%S")
                                except Exception:
                                    pass
                        if error_re.search(line):
                            has_error = True
            except Exception:
                continue

            # Finish time from file mtime
            mtime = os.path.getmtime(fpath)
            end_dt = datetime.fromtimestamp(mtime)

            dur_sec = None
            if start_dt and end_dt and end_dt >= start_dt:
                dur_sec = (end_dt - start_dt).total_seconds()

            records.append({
                "rule": rule,
                "target": target,
                "jobid": jobid_str,
                "start": start_dt.strftime("%Y-%m-%d %H:%M:%S") if start_dt else "N/A",
                "end": end_dt.strftime("%Y-%m-%d %H:%M:%S") if end_dt else "N/A",
                "duration_sec": dur_sec,
                "duration_fmt": format_duration(dur_sec),
                "state": "FAILED" if has_error else "COMPLETED",
                "source": f"logs/slurm/.../{fname}"
            })
    return records

def query_sacct_batch(jobids):
    """Batch-queries SLURM sacct for exact accounting durations."""
    if not jobids:
        return {}
    
    unique_ids = sorted(list(set(str(j) for j in jobids if str(j).isdigit())))
    if not unique_ids:
        return {}

    id_arg = ",".join(unique_ids)
    try:
        cmd = [
            "sacct", "-j", id_arg,
            "--format=JobID,JobName%35,Start,End,Elapsed,State,MaxRSS",
            "-P", "-n"
        ]
        out = subprocess.check_output(cmd, stderr=subprocess.DEVNULL, universal_newlines=True)
        sacct_map = {}
        for line in out.strip().split("\n"):
            if not line:
                continue
            parts = line.split("|")
            if len(parts) >= 6:
                jid = parts[0]
                if "." in jid: # Skip sub-steps (.batch, .extern)
                    continue
                jobname, start, end, elapsed, state = parts[1], parts[2], parts[3], parts[4], parts[5]
                
                # Parse elapsed HH:MM:SS or D-HH:MM:SS
                sec = parse_elapsed_str(elapsed)
                sacct_map[jid] = {
                    "jobname": jobname,
                    "start": start,
                    "end": end,
                    "elapsed_str": elapsed,
                    "elapsed_sec": sec,
                    "state": state
                }
        return sacct_map
    except Exception:
        return {}

def parse_elapsed_str(elapsed_str):
    try:
        days = 0
        if "-" in elapsed_str:
            d_part, time_part = elapsed_str.split("-")
            days = int(d_part)
        else:
            time_part = elapsed_str
        
        parts = time_part.split(":")
        if len(parts) == 3:
            h, m, s = int(parts[0]), int(parts[1]), int(float(parts[2]))
            return days * 86400 + h * 3600 + m * 60 + s
        elif len(parts) == 2:
            m, s = int(parts[0]), int(float(parts[1]))
            return days * 86400 + m * 60 + s
    except Exception:
        pass
    return None

def parse_rule_logs(workdir, filter_keyword):
    """Checks results/logs/ for files matching filter_keyword."""
    logs_dir = os.path.join(workdir, "results", "logs")
    records = []
    if not os.path.isdir(logs_dir):
        return records

    for root, dirs, files in os.walk(logs_dir):
        for fname in files:
            if not fname.endswith(".log") or filter_keyword not in fname:
                continue
            fpath = os.path.join(root, fname)
            rel_dir = os.path.relpath(root, logs_dir)
            rule = rel_dir.split(os.sep)[0]
            mtime = os.path.getmtime(fpath)
            records.append({
                "rule": rule,
                "target": fname.replace(".log", ""),
                "jobid": "N/A",
                "start": "N/A",
                "end": datetime.fromtimestamp(mtime).strftime("%Y-%m-%d %H:%M:%S"),
                "duration_sec": None,
                "duration_fmt": "N/A",
                "state": "COMPLETED",
                "source": os.path.relpath(fpath, workdir)
            })
    return records

def main():
    parser = argparse.ArgumentParser(description="Extract detailed execution time for each pipeline step")
    parser.add_argument("--condition", "-c", default="PL_B6xCAST", help="Condition or sample identifier (default: PL_B6xCAST)")
    parser.add_argument("--workdir", "-w", default=".", help="Workflow directory (default: .)")
    args = parser.parse_args()

    filter_kw = args.condition
    workdir = os.path.abspath(args.workdir)

    print(f"==========================================================================================")
    print(f" Step-by-Step Execution Time Report: {filter_kw}")
    print(f" Working Directory: {workdir}")
    print(f"==========================================================================================")

    # 1. Gather all Slurm job logs
    slurm_records = parse_slurm_logs(workdir, filter_kw)

    # 2. Query Slurm sacct for accounting data
    jobids = [r["jobid"] for r in slurm_records if r.get("jobid") and r["jobid"] != "N/A"]
    sacct_data = query_sacct_batch(jobids)

    # If sacct has entries, update duration from official Slurm accounting
    for r in slurm_records:
        jid = r["jobid"]
        if jid in sacct_data:
            s_info = sacct_data[jid]
            r["state"] = s_info["state"]
            if s_info["elapsed_sec"] is not None:
                r["duration_sec"] = s_info["elapsed_sec"]
                r["duration_fmt"] = format_duration(s_info["elapsed_sec"])
            r["source"] = f"sacct ({jid})"

    # 3. Fallback to results/logs/ if some steps are missing
    existing_rules = { (r["rule"], r["target"]) for r in slurm_records }
    rule_log_records = parse_rule_logs(workdir, filter_kw)
    for r in rule_log_records:
        if (r["rule"], r["target"]) not in existing_rules:
            slurm_records.append(r)

    # Filter out failed runs if a completed run exists for the same (rule, target)
    best_records = {}
    for r in slurm_records:
        key = (r["rule"], r["target"])
        if key not in best_records:
            best_records[key] = r
        else:
            prev = best_records[key]
            # Prefer COMPLETED over FAILED
            if prev["state"] == "FAILED" and r["state"] == "COMPLETED":
                best_records[key] = r
            # Or prefer whichever has duration
            elif r["duration_sec"] and not prev["duration_sec"]:
                best_records[key] = r

    # Sort records logically along the pipeline flow
    def sort_key(item):
        rule = item["rule"]
        order = STAGE_ORDER.index(rule) if rule in STAGE_ORDER else 99
        return (order, item.get("target", ""))

    sorted_records = sorted(best_records.values(), key=sort_key)

    if not sorted_records:
        print(f"\nNo SLURM or tool logs found for '{filter_kw}'.")
        print("Try running sacct directly:")
        print(f"  sacct --starttime 2026-08-15 --format=JobID,JobName%35,Start,End,Elapsed,State -n -X | grep {filter_kw}")
        return

    # Print Table
    header = f"{'Pipeline Step':<24} | {'Sample / Target':<28} | {'Elapsed Time':<13} | {'State':<10} | {'Source'}"
    print(header)
    print("-" * len(header))

    total_sec = 0
    for r in sorted_records:
        step_str = r['rule'][:24]
        target_str = r['target'][:28]
        dur_str = r['duration_fmt'][:13]
        state_str = r['state'][:10]
        source_str = r.get('source', '')[:25]
        if r.get('duration_sec'):
            total_sec += r['duration_sec']
        print(f"{step_str:<24} | {target_str:<28} | {dur_str:<13} | {state_str:<10} | {source_str}")

    print("-" * len(header))
    if total_sec > 0:
        print(f"Total Cumulative Compute Time: {format_duration(total_sec)}")
    print("==========================================================================================\n")

if __name__ == "__main__":
    main()
