# Micro-C & High-Resolution Contact Mapping Pipeline Specification

## 1. Objective
Design an snakemake pipeline to process **Micro-C** (default values, if any, should be those most appropiate for microC experiments) and related high-resolution chromosome conformation capture experiments (e.g., DNase Hi-C, HiChIP) on an HPC cluster managed by Slurm, without root/admin privileges.

The workflow must follow as closely as possible that described at https://data.4dnucleome.org/resources/data-analysis/hi_c-processing-pipeline and websites cited within (use latest versions whenever possible)

## 2. Modules available on the cluster
- apptainer/1.3.6
- python/3.14.7
- snakemake 9.25.1 (inside python module)