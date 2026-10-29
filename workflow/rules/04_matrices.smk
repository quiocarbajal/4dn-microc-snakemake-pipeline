# =============================================================================
# 04_matrices.smk: Multiresolution Matrix Generation (.mcool and .hic)
# Uses official 4DN v43 scripts inside official container
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Base Cooler Generation (1kb Base Resolution)
# -----------------------------------------------------------------------------
rule cooler_cload:
    input:
        pairs = f"{OUTDIR}/pairs/{{sample}}.ff.pairs.gz",
        px2 = f"{OUTDIR}/pairs/{{sample}}.ff.pairs.gz.px2",
        chrom_sizes = get_chrom_sizes()
    output:
        cool = f"{OUTDIR}/cool/{{sample}}.cool"
    log:
        f"{OUTDIR}/logs/cooler/cload_{{sample}}.log"
    params:
        outprefix = lambda wildcards, output: output.cool.replace('.cool', ''),
        bin_size = BASE_RES,
        max_split = 2
    threads:
        4
    # run-cooler.sh copies the chrom sizes to a fixed-name ./tmpchrsize; a per-job
    # shadow directory stops a concurrent job reading it while it is being rewritten
    shadow:
        "minimal"
    container:
        CONTAINER_4DN
    shell:
        """
        mkdir -p $(dirname {output.cool}) $(dirname {log})
        /usr/local/bin/run-cooler.sh \
            {input.pairs} \
            {input.chrom_sizes} \
            {params.bin_size} \
            {threads} \
            {params.outprefix} \
            {params.max_split} > {log} 2>&1
        """

# -----------------------------------------------------------------------------
# 2. Multiresolution Cooler (.mcool) with ICE Balancing
# -----------------------------------------------------------------------------
if config.get("matrix", {}).get("generate_hic", False) and config.get("matrix", {}).get("add_juicer_norms_to_mcool", False):
    rule cooler_zoomify:
        input:
            cool = f"{OUTDIR}/cool/{{sample}}.cool"
        output:
            mcool_raw = f"{OUTDIR}/mcool_raw/{{sample}}.multires.cool"
        log:
            f"{OUTDIR}/logs/cooler/zoomify_{{sample}}.log"
        params:
            outprefix = lambda wildcards, output: output.mcool_raw.replace('.multires.cool', ''),
            resolutions = RESOLUTIONS_STR,
            chunksize = 10000000,
            no_balance_flag = "" if config.get("matrix", {}).get("balance", True) else "-B"
        threads:
            8
        container:
            CONTAINER_4DN
        shell:
            """
            mkdir -p $(dirname {output.mcool_raw}) $(dirname {log})
            /usr/local/bin/run-cool2multirescool.sh \
                -i {input.cool} \
                -o {params.outprefix} \
                -p {threads} \
                -c {params.chunksize} \
                -u {params.resolutions} \
                {params.no_balance_flag} > {log} 2>&1
            """

    rule add_hic_norm_to_mcool:
        input:
            hic = f"{OUTDIR}/hic/{{sample}}.hic",
            mcool_raw = f"{OUTDIR}/mcool_raw/{{sample}}.multires.cool"
        output:
            mcool = f"{OUTDIR}/mcool/{{sample}}.mcool"
        log:
            f"{OUTDIR}/logs/cooler/add_hic_norm_{{sample}}.log"
        params:
            outdir = lambda wildcards, output: f"{os.path.dirname(output.mcool)}/tmp_{wildcards.sample}"
        threads:
            2
        container:
            CONTAINER_4DN
        shell:
            """
            mkdir -p {params.outdir} $(dirname {log})
            /usr/local/bin/run-add-hicnormvector-to-mcool.sh \
                {input.hic} \
                {input.mcool_raw} \
                {params.outdir} > {log} 2>&1
            mv {params.outdir}/output.mcool {output.mcool}
            rm -rf {params.outdir}
            """
else:
    rule cooler_zoomify:
        input:
            cool = f"{OUTDIR}/cool/{{sample}}.cool"
        output:
            mcool = f"{OUTDIR}/mcool/{{sample}}.mcool"
        log:
            f"{OUTDIR}/logs/cooler/zoomify_{{sample}}.log"
        params:
            outprefix = lambda wildcards, output: f"{os.path.dirname(output.mcool)}/{wildcards.sample}",
            resolutions = RESOLUTIONS_STR,
            chunksize = 10000000,
            no_balance_flag = "" if config.get("matrix", {}).get("balance", True) else "-B"
        threads:
            8
        container:
            CONTAINER_4DN
        shell:
            """
            mkdir -p $(dirname {output.mcool}) $(dirname {log})
            /usr/local/bin/run-cool2multirescool.sh \
                -i {input.cool} \
                -o {params.outprefix} \
                -p {threads} \
                -c {params.chunksize} \
                -u {params.resolutions} \
                {params.no_balance_flag} > {log} 2>&1
            mv {params.outprefix}.multires.cool {output.mcool}
            """

# -----------------------------------------------------------------------------
# 3. Juicebox Matrix (.hic) Generation with Juicer Tools
# -----------------------------------------------------------------------------
rule juicebox_pre:
    input:
        pairs = f"{OUTDIR}/pairs/{{sample}}.ff.pairs.gz",
        chrom_sizes = get_chrom_sizes()
    output:
        hic = f"{OUTDIR}/hic/{{sample}}.hic"
    log:
        f"{OUTDIR}/logs/juicebox/pre_{{sample}}.log"
    params:
        outprefix = lambda wildcards, output: output.hic.replace('.hic', ''),
        min_res = BASE_RES,
        resolutions = RESOLUTIONS_STR,
        maxmem = lambda wildcards, resources: f"{max(1, int(resources.mem_mb // 1024))}g",
        no_balance_flag = "" if config.get("matrix", {}).get("balance", True) else "-B"
    threads:
        2
    # Default so params.maxmem works without a profile; the Slurm profile's set-resources overrides it
    resources:
        mem_mb = 65536
    container:
        CONTAINER_4DN
    shell:
        """
        mkdir -p $(dirname {output.hic}) $(dirname {log})
        /usr/local/bin/run-juicebox-pre.sh \
            -i {input.pairs} \
            -c {input.chrom_sizes} \
            -o {params.outprefix} \
            -r {params.min_res} \
            -u {params.resolutions} \
            -m {params.maxmem} \
            {params.no_balance_flag} > {log} 2>&1
        """

# -----------------------------------------------------------------------------
# 4. Merged Biological Replicate Matrices (Official 4DN hi-c-processing-pairs.cwl)
# -----------------------------------------------------------------------------
rule cooler_cload_merged:
    input:
        pairs = f"{OUTDIR}/pairs_merged/{{condition}}.merged.ff.pairs.gz",
        px2 = f"{OUTDIR}/pairs_merged/{{condition}}.merged.ff.pairs.gz.px2",
        chrom_sizes = get_chrom_sizes()
    output:
        cool = f"{OUTDIR}/cool_merged/{{condition}}.cool"
    log:
        f"{OUTDIR}/logs/cooler/cload_merged_{{condition}}.log"
    params:
        outprefix = lambda wildcards, output: output.cool.replace('.cool', ''),
        bin_size = BASE_RES,
        max_split = 2
    threads:
        4
    # run-cooler.sh copies the chrom sizes to a fixed-name ./tmpchrsize; a per-job
    # shadow directory stops a concurrent job reading it while it is being rewritten
    shadow:
        "minimal"
    container:
        CONTAINER_4DN
    shell:
        """
        mkdir -p $(dirname {output.cool}) $(dirname {log})
        /usr/local/bin/run-cooler.sh \
            {input.pairs} \
            {input.chrom_sizes} \
            {params.bin_size} \
            {threads} \
            {params.outprefix} \
            {params.max_split} > {log} 2>&1
        """

rule cooler_zoomify_merged:
    input:
        cool = f"{OUTDIR}/cool_merged/{{condition}}.cool"
    output:
        mcool = f"{OUTDIR}/mcool_merged/{{condition}}.mcool"
    log:
        f"{OUTDIR}/logs/cooler/zoomify_merged_{{condition}}.log"
    params:
        outprefix = lambda wildcards, output: f"{os.path.dirname(output.mcool)}/{wildcards.condition}",
        resolutions = RESOLUTIONS_STR,
        chunksize = 10000000,
        no_balance_flag = "" if config.get("matrix", {}).get("balance", True) else "-B"
    threads:
        8
    container:
        CONTAINER_4DN
    shell:
        """
        mkdir -p $(dirname {output.mcool}) $(dirname {log})
        /usr/local/bin/run-cool2multirescool.sh \
            -i {input.cool} \
            -o {params.outprefix} \
            -p {threads} \
            -c {params.chunksize} \
            -u {params.resolutions} \
            {params.no_balance_flag} > {log} 2>&1
        mv {params.outprefix}.multires.cool {output.mcool}
        """

rule juicebox_pre_merged:
    input:
        pairs = f"{OUTDIR}/pairs_merged/{{condition}}.merged.ff.pairs.gz",
        chrom_sizes = get_chrom_sizes()
    output:
        hic = f"{OUTDIR}/hic_merged/{{condition}}.hic"
    log:
        f"{OUTDIR}/logs/juicebox/pre_merged_{{condition}}.log"
    params:
        outprefix = lambda wildcards, output: output.hic.replace('.hic', ''),
        min_res = BASE_RES,
        resolutions = RESOLUTIONS_STR,
        maxmem = lambda wildcards, resources: f"{max(1, int(resources.mem_mb // 1024))}g",
        no_balance_flag = "" if config.get("matrix", {}).get("balance", True) else "-B"
    threads:
        2
    # Default so params.maxmem works without a profile; the Slurm profile's set-resources overrides it
    resources:
        mem_mb = 65536
    container:
        CONTAINER_4DN
    shell:
        """
        mkdir -p $(dirname {output.hic}) $(dirname {log})
        /usr/local/bin/run-juicebox-pre.sh \
            -i {input.pairs} \
            -c {input.chrom_sizes} \
            -o {params.outprefix} \
            -r {params.min_res} \
            -u {params.resolutions} \
            -m {params.maxmem} \
            {params.no_balance_flag} > {log} 2>&1
        """

# -----------------------------------------------------------------------------
# 5. Raichu Normalization (In-Place)
# -----------------------------------------------------------------------------
rule raichu_normalization:
    input:
        mcool = f"{OUTDIR}/mcool/{{sample}}.mcool"
    output:
        done = touch(f"{OUTDIR}/mcool/{{sample}}.raichu_done")
    log:
        f"{OUTDIR}/logs/raichu/raichu_{{sample}}.log"
    params:
        resolutions = config.get("matrix", {}).get("raichu_resolutions", [5000, 10000]),
        venv = config.get("matrix", {}).get("raichu_venv", "/hpc-prj/pezza/Agustin/environments/raichu")
    # Matches the lab's standalone raichu jobs: raichu -p 2
    # (memory: 5 GB via set-resources in the Slurm profile)
    threads:
        2
    shell:
        """
        module purge
        module load slurm python/3.10.2
        unset PYTHONPATH PYTHONHOME
        export PYTHONNOUSERSITE=1
        export HDF5_USE_FILE_LOCKING=FALSE
        source {params.venv}/bin/activate
        
        for res in {params.resolutions}; do
            echo "Running raichu for resolution $res on {input.mcool}" >> {log}
            raichu --cool-uri "{input.mcool}::/resolutions/$res" --window-size 200 -p {threads} -n raichu_weight -f >> {log} 2>&1
        done
        """

rule raichu_normalization_merged:
    input:
        mcool = f"{OUTDIR}/mcool_merged/{{condition}}.mcool"
    output:
        done = touch(f"{OUTDIR}/mcool_merged/{{condition}}.raichu_done")
    log:
        f"{OUTDIR}/logs/raichu/raichu_merged_{{condition}}.log"
    params:
        resolutions = config.get("matrix", {}).get("raichu_resolutions", [5000, 10000]),
        venv = config.get("matrix", {}).get("raichu_venv", "/hpc-prj/pezza/Agustin/environments/raichu")
    # Matches the lab's standalone raichu jobs: raichu -p 2
    # (memory: 32 GB via set-resources in the Slurm profile; 5 GB ran out on a merged condition)
    threads:
        2
    shell:
        """
        module purge
        module load slurm python/3.10.2
        unset PYTHONPATH PYTHONHOME
        export PYTHONNOUSERSITE=1
        export HDF5_USE_FILE_LOCKING=FALSE
        source {params.venv}/bin/activate
        
        for res in {params.resolutions}; do
            echo "Running raichu for resolution $res on {input.mcool}" >> {log}
            raichu --cool-uri "{input.mcool}::/resolutions/$res" --window-size 200 -p {threads} -n raichu_weight -f >> {log} 2>&1
        done
        """
