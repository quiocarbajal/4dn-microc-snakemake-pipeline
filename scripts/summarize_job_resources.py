#!/usr/bin/env python3
"""
Summarize Slurm Job Resource Utilization (CPU & Memory) for Micro-C Pipeline
=============================================================================
Extracts JobIDs for all processed samples and rules, runs `seff` (and `sacct`
as verification/fallback), and computes:
  1. Per-sample/per-rule resource metrics (cores, walltime, CPU efficiency,
     memory used vs allocated, memory efficiency).
  2. Aggregated rule-level summary with actionable recommendations for tuning
     `threads` and `mem_mb` in `workflow/profiles/slurm/config.yaml`.

Usage on cluster:
  python3 scripts/summarize_job_resources.py
  python3 scripts/summarize_job_resources.py --tsv resources_summary.tsv
  python3 scripts/summarize_job_resources.py --rule bwa_mem
  python3 scripts/summarize_job_resources.py --jobids 6189336 6189337 6189729
"""

import os
import sys
import re
import csv
import json
import glob
import argparse
import subprocess
from datetime import timedelta
from collections import defaultdict

# -----------------------------------------------------------------------------
# Unit conversion and formatting helpers
# -----------------------------------------------------------------------------

def parse_size_to_mb(val_str):
    """Parses strings like '12.34 GB', '512 MB', '1048576 KB' into megabytes (float)."""
    if not val_str or val_str in ("N/A", "---", "0", ""):
        return 0.0
    val_str = str(val_str).strip()
    m = re.match(r"^([0-9.]+)\s*([A-Za-z]*)", val_str)
    if not m:
        return 0.0
    num = float(m.group(1))
    unit = m.group(2).upper()
    if unit in ("K", "KB", "KIB"):
        return num / 1024.0
    elif unit in ("M", "MB", "MIB", ""):
        return num
    elif unit in ("G", "GB", "GIB"):
        return num * 1024.0
    elif unit in ("T", "TB", "TIB"):
        return num * 1024.0 * 1024.0
    return num

def format_mb(mb):
    """Formats megabytes nicely into MB or GB."""
    if mb is None or mb <= 0:
        return "N/A"
    if mb >= 1024.0:
        return f"{mb / 1024.0:.2f} GB"
    return f"{mb:.1f} MB"

def parse_time_to_sec(time_str):
    """Parses [[D-]HH:]MM:SS or HH:MM:SS into total seconds."""
    if not time_str or time_str in ("N/A", "---", ""):
        return 0
    days = 0
    clean_str = time_str.strip()
    if "-" in clean_str:
        d_part, clean_str = clean_str.split("-", 1)
        try:
            days = int(d_part)
        except ValueError:
            pass
    parts = clean_str.split(":")
    try:
        if len(parts) == 3:
            h, m, s = int(parts[0]), int(parts[1]), int(float(parts[2]))
            return days * 86400 + h * 3600 + m * 60 + s
        elif len(parts) == 2:
            m, s = int(parts[0]), int(float(parts[1]))
            return days * 86400 + m * 60 + s
        elif len(parts) == 1:
            return days * 86400 + int(float(parts[0]))
    except Exception:
        pass
    return 0

def format_duration(seconds):
    """Formats total seconds into HH:MM:SS or MM:SS."""
    if seconds is None or seconds <= 0:
        return "00:00"
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

# -----------------------------------------------------------------------------
# Job Discovery
# -----------------------------------------------------------------------------

def discover_jobs_from_slurm_logs(log_dir):
    """
    Finds job IDs from logs/slurm/rule_<rulename>/<target>/<jobid>.log
    Returns list of dicts: [{'rule': ..., 'target': ..., 'jobid': ...}]
    """
    records = []
    if not os.path.isdir(log_dir):
        return records

    for root, dirs, files in os.walk(log_dir):
        for fname in files:
            if not fname.endswith(".log"):
                continue
            base, _ = os.path.splitext(fname)
            if not base.isdigit():
                continue
            jobid = int(base)
            rel_dir = os.path.relpath(root, log_dir)
            parts = rel_dir.split(os.sep)
            rule = parts[0].replace("rule_", "") if len(parts) >= 1 else "unknown"
            target = parts[1] if len(parts) >= 2 else "global"
            records.append({
                "rule": rule,
                "target": target,
                "jobid": jobid,
                "source": os.path.join(root, fname)
            })
    return records

def discover_jobs_from_text_logs(file_paths):
    """
    Scans Snakemake log files, tmux output, etc. for submitted SLURM job lines:
      Job <N> has been submitted with SLURM jobid <jobid> (log: .../rule_<rule>/<target>/<jobid>.log)
    """
    records = []
    subm_re = re.compile(r"Job \d+ has been submitted with SLURM jobid (\d+) \(log: .*?rule_([^/]+)/([^/]+)/(\d+)\.log\)")
    
    for fpath in file_paths:
        if not os.path.isfile(fpath):
            continue
        try:
            with open(fpath, "r", encoding="utf-8", errors="ignore") as f:
                for line in f:
                    m = subm_re.search(line)
                    if m:
                        jobid = int(m.group(1))
                        rule = m.group(2)
                        target = m.group(3)
                        records.append({
                            "rule": rule,
                            "target": target,
                            "jobid": jobid,
                            "source": fpath
                        })
        except Exception:
            pass
    return records

# -----------------------------------------------------------------------------
# Live Slurm Querying (seff & sacct)
# -----------------------------------------------------------------------------

def query_seff(jobid):
    """
    Runs `seff <jobid>` and parses its standard output.
    Returns dict of parsed metrics or None if seff failed/unavailable.
    """
    try:
        cmd = ["seff", str(jobid)]
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=15)
        if res.returncode != 0 and not res.stdout:
            return None
        output = res.stdout
    except Exception:
        return None

    data = {
        "jobid": str(jobid),
        "state": "UNKNOWN",
        "cores": None,
        "cpu_utilized_str": "N/A",
        "cpu_utilized_sec": 0,
        "cpu_efficiency_pct": None,
        "walltime_str": "N/A",
        "walltime_sec": 0,
        "memory_utilized_str": "N/A",
        "memory_utilized_mb": 0.0,
        "memory_allocated_str": "N/A",
        "memory_allocated_mb": 0.0,
        "memory_efficiency_pct": None
    }

    for line in output.splitlines():
        line = line.strip()
        if line.startswith("State:"):
            data["state"] = line.replace("State:", "").strip()
        elif line.startswith("Cores:") or line.startswith("Cores per node:"):
            m = re.search(r"(\d+)", line)
            if m:
                data["cores"] = int(m.group(1))
        elif line.startswith("CPU Utilized:"):
            val = line.replace("CPU Utilized:", "").strip()
            data["cpu_utilized_str"] = val
            data["cpu_utilized_sec"] = parse_time_to_sec(val)
        elif line.startswith("CPU Efficiency:"):
            # e.g. "CPU Efficiency: 85.20% of 02:39:00 core-walltime"
            m = re.search(r"([0-9.]+)%", line)
            if m:
                data["cpu_efficiency_pct"] = float(m.group(1))
        elif line.startswith("Job Wall-clock time:"):
            val = line.replace("Job Wall-clock time:", "").strip()
            data["walltime_str"] = val
            data["walltime_sec"] = parse_time_to_sec(val)
        elif line.startswith("Memory Utilized:"):
            val = line.replace("Memory Utilized:", "").strip()
            data["memory_utilized_str"] = val
            data["memory_utilized_mb"] = parse_size_to_mb(val)
        elif line.startswith("Memory Efficiency:"):
            # e.g. "Memory Efficiency: 37.66% of 32.77 GB"
            m_pct = re.search(r"([0-9.]+)%", line)
            if m_pct:
                data["memory_efficiency_pct"] = float(m_pct.group(1))
            m_alloc = re.search(r"of\s+([0-9.]+\s*[A-Za-z]+)", line)
            if m_alloc:
                data["memory_allocated_str"] = m_alloc.group(1)
                data["memory_allocated_mb"] = parse_size_to_mb(m_alloc.group(1))

    return data

def query_sacct_batch(jobids):
    """
    Fallback/accelerator: batch queries `sacct` for detailed metrics if seff is slow.
    """
    if not jobids:
        return {}
    id_list = sorted(list(set(str(j) for j in jobids)))
    sacct_dict = {}
    
    # Process in chunks of 100 to avoid command length limits
    chunk_size = 100
    for i in range(0, len(id_list), chunk_size):
        chunk = id_list[i:i+chunk_size]
        cmd = [
            "sacct", "-j", ",".join(chunk),
            "--format=JobID,JobName%35,State,ExitCode,Elapsed,AllocCPUS,TotalCPU,ReqMem,MaxRSS",
            "-P", "-n"
        ]
        try:
            res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=30)
            for line in res.stdout.splitlines():
                if not line.strip():
                    continue
                parts = line.strip().split("|")
                if len(parts) >= 9:
                    raw_jid, jobname, state, exitcode, elapsed, cpus, total_cpu, req_mem, max_rss = parts[:9]
                    parent_jid = raw_jid.split(".")[0]
                    if parent_jid not in sacct_dict:
                        sacct_dict[parent_jid] = {
                            "jobid": parent_jid,
                            "state": state,
                            "elapsed": elapsed,
                            "elapsed_sec": parse_time_to_sec(elapsed),
                            "cpus": int(cpus) if cpus.isdigit() else None,
                            "total_cpu": total_cpu,
                            "total_cpu_sec": parse_time_to_sec(total_cpu),
                            "req_mem": req_mem,
                            "req_mem_mb": parse_size_to_mb(req_mem),
                            "max_rss": max_rss,
                            "max_rss_mb": parse_size_to_mb(max_rss)
                        }
                    else:
                        # Update peak RSS from sub-steps (e.g. .batch or .0)
                        rss_mb = parse_size_to_mb(max_rss)
                        if rss_mb > sacct_dict[parent_jid]["max_rss_mb"]:
                            sacct_dict[parent_jid]["max_rss_mb"] = rss_mb
                            sacct_dict[parent_jid]["max_rss"] = max_rss
        except Exception:
            pass
    return sacct_dict

# -----------------------------------------------------------------------------
# Recommendation Engine
# -----------------------------------------------------------------------------

def generate_rule_recommendation(alloc_cores, mean_cpu_eff, alloc_mem_mb, peak_mem_mb, mean_mem_eff):
    """
    Produces actionable rule-level tuning advice.
    """
    recs = []

    # Memory recommendation
    if peak_mem_mb and alloc_mem_mb and alloc_mem_mb > 0:
        ratio = peak_mem_mb / alloc_mem_mb
        if ratio < 0.25:
            # Over-allocated by 4x+
            safe_alloc_mb = int(max(4096, peak_mem_mb * 1.6))
            safe_alloc_mb = ((safe_alloc_mb + 1023) // 1024) * 1024
            recs.append(f"RAM over-allocated ({ratio*100:.1f}% used). Can reduce from {format_mb(alloc_mem_mb)} to ~{format_mb(safe_alloc_mb)}.")
        elif ratio > 0.85:
            safe_alloc_mb = int(alloc_mem_mb * 1.5)
            safe_alloc_mb = ((safe_alloc_mb + 1023) // 1024) * 1024
            recs.append(f"RAM near limit ({ratio*100:.1f}% peak). Recommend increasing to ~{format_mb(safe_alloc_mb)} to prevent OOM.")
        else:
            recs.append(f"RAM sizing is well-balanced ({ratio*100:.1f}% peak).")
    
    # CPU recommendation
    if alloc_cores and mean_cpu_eff is not None:
        if alloc_cores > 4 and mean_cpu_eff < 35.0:
            suggested_cores = max(2, int(alloc_cores * (mean_cpu_eff / 100.0) * 1.5))
            recs.append(f"CPU under-utilized ({mean_cpu_eff:.1f}% eff across {alloc_cores} cores). Reduce threads to {suggested_cores} to decrease queue wait time.")
        elif alloc_cores > 1 and mean_cpu_eff > 85.0:
            recs.append(f"CPU highly saturated ({mean_cpu_eff:.1f}% eff). Scaling threads up may improve throughput if I/O permits.")
        else:
            recs.append(f"CPU allocation is appropriate ({mean_cpu_eff:.1f}% eff).")

    return " | ".join(recs) if recs else "Sufficient data pending."

# -----------------------------------------------------------------------------
# Main Execution Flow
# -----------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Summarize Slurm job CPU and memory utilization across Micro-C rules and samples."
    )
    parser.add_argument(
        "--slurm-logs", default="logs/slurm",
        help="Path to logs/slurm directory (default: logs/slurm)"
    )
    parser.add_argument(
        "--snakemake-logs", nargs="*", default=["microC_tmux_output.txt", "2026-09-11T213125.136235.snakemake.log"],
        help="Additional Snakemake log files to scan for job IDs"
    )
    parser.add_argument(
        "--jobids", nargs="*", type=int, default=None,
        help="Explicit list of Slurm Job IDs to query"
    )
    parser.add_argument(
        "--rule", default=None,
        help="Filter summary to a specific rule (e.g. bwa_mem, cooler_zoomify)"
    )
    parser.add_argument(
        "--sample", default=None,
        help="Filter summary to a specific sample or biological replicate"
    )
    parser.add_argument(
        "--tsv", default=None,
        help="Save detailed sample metrics to a TSV file"
    )
    parser.add_argument(
        "--csv", default=None,
        help="Save detailed sample metrics to a CSV file"
    )
    parser.add_argument(
        "--summary-tsv", default=None,
        help="Save aggregated rule summary table to a TSV file"
    )
    parser.add_argument(
        "--use-sacct-only", action="store_true",
        help="Skip calling `seff` and rely exclusively on fast batch `sacct`"
    )
    args = parser.parse_args()

    print("================================================================================")
    print(" 4DN Micro-C Pipeline: Slurm Resource Utilization & Optimization Auditor")
    print("================================================================================")

    # 1. Discover jobs
    discovered = []
    if args.jobids:
        for jid in args.jobids:
            discovered.append({"rule": "manual", "target": "manual", "jobid": jid, "source": "CLI"})
    else:
        # Scan logs/slurm
        from_dir = discover_jobs_from_slurm_logs(args.slurm_logs)
        discovered.extend(from_dir)
        
        # Scan text logs
        from_txt = discover_jobs_from_text_logs(args.snakemake_logs)
        discovered.extend(from_txt)

    # Deduplicate by (rule, target, jobid)
    unique_map = {}
    for d in discovered:
        k = (d["rule"], d["target"], d["jobid"])
        if k not in unique_map:
            unique_map[k] = d

    all_jobs = list(unique_map.values())

    if not all_jobs:
        print("\n[!] No Slurm job IDs found in logs/slurm or specified log files.")
        print("    Ensure you run this script in the root directory of the pipeline on the cluster.")
        sys.exit(0)

    # Apply filters
    if args.rule:
        all_jobs = [j for j in all_jobs if args.rule.lower() in j["rule"].lower()]
    if args.sample:
        all_jobs = [j for j in all_jobs if args.sample.lower() in j["target"].lower()]

    print(f"[*] Discovered {len(all_jobs)} unique job execution records across pipeline rules.")

    # 2. Query resource data (seff and/or sacct)
    print("[*] Querying cluster accounting database (`seff` / `sacct`)...")
    
    # Query sacct batch first as a base/fallback
    all_jids = [j["jobid"] for j in all_jobs]
    sacct_data = query_sacct_batch(all_jids)

    job_records = []
    for idx, item in enumerate(all_jobs, 1):
        jid = item["jobid"]
        rule = item["rule"]
        target = item["target"]

        metric = None
        if not args.use_sacct_only:
            metric = query_seff(jid)
        
        # Merge with sacct if seff was unavailable or incomplete
        s_item = sacct_data.get(str(jid))
        if metric is None or metric["state"] == "UNKNOWN":
            if s_item:
                cores = s_item["cpus"]
                wall_sec = s_item["elapsed_sec"]
                cpu_sec = s_item["total_cpu_sec"]
                cpu_eff = (cpu_sec / (wall_sec * cores) * 100.0) if (wall_sec and cores and cores > 0) else None
                req_mb = s_item["req_mem_mb"]
                max_mb = s_item["max_rss_mb"]
                mem_eff = (max_mb / req_mb * 100.0) if (req_mb and req_mb > 0) else None
                
                metric = {
                    "jobid": str(jid),
                    "state": s_item["state"],
                    "cores": cores,
                    "cpu_utilized_str": s_item["total_cpu"],
                    "cpu_utilized_sec": cpu_sec,
                    "cpu_efficiency_pct": cpu_eff,
                    "walltime_str": s_item["elapsed"],
                    "walltime_sec": wall_sec,
                    "memory_utilized_str": format_mb(max_mb),
                    "memory_utilized_mb": max_mb,
                    "memory_allocated_str": format_mb(req_mb),
                    "memory_allocated_mb": req_mb,
                    "memory_efficiency_pct": mem_eff
                }
            else:
                metric = {
                    "jobid": str(jid),
                    "state": "UNRETRIEVED",
                    "cores": None,
                    "cpu_utilized_str": "N/A",
                    "cpu_utilized_sec": 0,
                    "cpu_efficiency_pct": None,
                    "walltime_str": "N/A",
                    "walltime_sec": 0,
                    "memory_utilized_str": "N/A",
                    "memory_utilized_mb": 0.0,
                    "memory_allocated_str": "N/A",
                    "memory_allocated_mb": 0.0,
                    "memory_efficiency_pct": None
                }

        rec = {
            "rule": rule,
            "target": target,
            "jobid": jid,
            "state": metric.get("state", "UNKNOWN"),
            "cores": metric.get("cores"),
            "walltime_sec": metric.get("walltime_sec", 0),
            "walltime_str": metric.get("walltime_str", "N/A"),
            "cpu_utilized_str": metric.get("cpu_utilized_str", "N/A"),
            "cpu_eff_pct": metric.get("cpu_efficiency_pct"),
            "mem_used_mb": metric.get("memory_utilized_mb", 0.0),
            "mem_used_str": metric.get("memory_utilized_str", "N/A"),
            "mem_alloc_mb": metric.get("memory_allocated_mb", 0.0),
            "mem_alloc_str": metric.get("memory_allocated_str", "N/A"),
            "mem_eff_pct": metric.get("memory_efficiency_pct")
        }
        job_records.append(rec)

    # 3. Display Detail Table
    print("\n" + "="*120)
    print(f"{'RULE':<24} {'SAMPLE / TARGET':<25} {'JOBID':<9} {'STATE':<11} {'WALLTIME':<10} {'CORES':<6} {'CPU EFF':<9} {'MEM USED':<11} {'MEM ALLOC':<11} {'MEM EFF'}")
    print("="*120)
    
    # Sort by rule, target
    job_records.sort(key=lambda x: (x["rule"], x["target"], x["jobid"]))

    for r in job_records:
        cpu_eff_str = f"{r['cpu_eff_pct']:.1f}%" if r["cpu_eff_pct"] is not None else "N/A"
        mem_eff_str = f"{r['mem_eff_pct']:.1f}%" if r["mem_eff_pct"] is not None else "N/A"
        cores_str = str(r["cores"]) if r["cores"] is not None else "N/A"
        print(f"{r['rule'][:23]:<24} {r['target'][:24]:<25} {r['jobid']:<9} {r['state'][:10]:<11} {r['walltime_str']:<10} {cores_str:<6} {cpu_eff_str:<9} {r['mem_used_str']:<11} {r['mem_alloc_str']:<11} {mem_eff_str}")

    # 4. Aggregated Summary per Rule
    rule_groups = defaultdict(list)
    for r in job_records:
        rule_groups[r["rule"]].append(r)

    print("\n" + "="*120)
    print(" AGGREGATED RULE RESOURCE SUMMARY & SIZING RECOMMENDATIONS")
    print("="*120)
    print(f"{'RULE':<24} {'JOBS':<6} {'CORES':<6} {'AVG CPU EFF':<13} {'ALLOC RAM':<11} {'AVG RAM USED':<14} {'PEAK RAM':<12} {'AVG RAM EFF':<12} {'MAX WALLTIME'}")
    print("-"*120)

    summary_rows = []
    for rule, items in sorted(rule_groups.items()):
        completed = [x for x in items if "COMPLETED" in x["state"].upper()]
        eval_set = completed if completed else items
        
        n_jobs = len(items)
        cores = eval_set[0]["cores"] if eval_set and eval_set[0]["cores"] else "N/A"
        alloc_ram_mb = max((x["mem_alloc_mb"] for x in eval_set), default=0.0)
        
        cpu_effs = [x["cpu_eff_pct"] for x in eval_set if x["cpu_eff_pct"] is not None]
        avg_cpu_eff = sum(cpu_effs) / len(cpu_effs) if cpu_effs else None
        
        used_mems = [x["mem_used_mb"] for x in eval_set if x["mem_used_mb"] > 0]
        avg_used_mem = sum(used_mems) / len(used_mems) if used_mems else 0.0
        peak_used_mem = max(used_mems) if used_mems else 0.0

        mem_effs = [x["mem_eff_pct"] for x in eval_set if x["mem_eff_pct"] is not None]
        avg_mem_eff = sum(mem_effs) / len(mem_effs) if mem_effs else None

        walltimes = [x["walltime_sec"] for x in eval_set if x["walltime_sec"] > 0]
        max_wall_sec = max(walltimes) if walltimes else 0

        avg_cpu_str = f"{avg_cpu_eff:.1f}%" if avg_cpu_eff is not None else "N/A"
        alloc_ram_str = format_mb(alloc_ram_mb) if alloc_ram_mb > 0 else "N/A"
        avg_used_str = format_mb(avg_used_mem) if avg_used_mem > 0 else "N/A"
        peak_used_str = format_mb(peak_used_mem) if peak_used_mem > 0 else "N/A"
        avg_mem_str = f"{avg_mem_eff:.1f}%" if avg_mem_eff is not None else "N/A"
        max_wall_str = format_duration(max_wall_sec)

        cores_int = cores if isinstance(cores, int) else None
        recommendation = generate_rule_recommendation(cores_int, avg_cpu_eff, alloc_ram_mb, peak_used_mem, avg_mem_eff)

        print(f"{rule[:23]:<24} {n_jobs:<6} {str(cores):<6} {avg_cpu_str:<13} {alloc_ram_str:<11} {avg_used_str:<14} {peak_used_str:<12} {avg_mem_str:<12} {max_wall_str}")
        print(f"  └─► RECOMMENDATION: {recommendation}\n")

        summary_rows.append({
            "rule": rule,
            "total_jobs": n_jobs,
            "completed_jobs": len(completed),
            "allocated_cores": cores,
            "avg_cpu_efficiency_pct": round(avg_cpu_eff, 2) if avg_cpu_eff is not None else "N/A",
            "allocated_ram_mb": round(alloc_ram_mb, 1),
            "avg_used_ram_mb": round(avg_used_mem, 1),
            "peak_used_ram_mb": round(peak_used_mem, 1),
            "avg_mem_efficiency_pct": round(avg_mem_eff, 2) if avg_mem_eff is not None else "N/A",
            "max_walltime_sec": max_wall_sec,
            "max_walltime_str": max_wall_str,
            "recommendation": recommendation
        })

    # 5. Export files if requested
    if args.tsv:
        save_tsv(job_records, args.tsv)
        print(f"[+] Saved detailed sample resource records to: {args.tsv}")
    if args.csv:
        save_csv(job_records, args.csv)
        print(f"[+] Saved detailed sample resource records to: {args.csv}")
    if args.summary_tsv:
        save_tsv(summary_rows, args.summary_tsv)
        print(f"[+] Saved aggregated rule recommendations to: {args.summary_tsv}")

    print("\n[✔] Audit complete. Use the recommendations above to update `workflow/profiles/slurm/config.yaml`.")

def save_tsv(records, path):
    if not records:
        return
    keys = list(records[0].keys())
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=keys, delimiter="\t")
        writer.writeheader()
        writer.writerows(records)

def save_csv(records, path):
    if not records:
        return
    keys = list(records[0].keys())
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(records)

if __name__ == "__main__":
    main()
