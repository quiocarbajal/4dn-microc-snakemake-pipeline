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
    params:
        outprefix = f"{OUTDIR}/cool/{{sample}}",
        bin_size = BASE_RES,
        max_split = 2
    threads:
        8
    container:
        CONTAINER_4DN
    shell:
        """
        mkdir -p $(dirname {output.cool})
        /usr/local/bin/run-cooler.sh \
            {input.pairs} \
            {input.chrom_sizes} \
            {params.bin_size} \
            {threads} \
            {params.outprefix} \
            {params.max_split}
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
        params:
            outprefix = f"{OUTDIR}/mcool_raw/{{sample}}",
            resolutions = RESOLUTIONS_STR,
            chunksize = 10000000,
            no_balance_flag = "" if config.get("matrix", {}).get("balance", True) else "-B"
        threads:
            8
        container:
            CONTAINER_4DN
        shell:
            """
            mkdir -p $(dirname {output.mcool_raw})
            /usr/local/bin/run-cool2multirescool.sh \
                -i {input.cool} \
                -o {params.outprefix} \
                -p {threads} \
                -c {params.chunksize} \
                -u {params.resolutions} \
                {params.no_balance_flag}
            """

    rule add_hic_norm_to_mcool:
        input:
            hic = f"{OUTDIR}/hic/{{sample}}.hic",
            mcool_raw = f"{OUTDIR}/mcool_raw/{{sample}}.multires.cool"
        output:
            mcool = f"{OUTDIR}/mcool/{{sample}}.mcool"
        params:
            outdir = f"{OUTDIR}/mcool/tmp_{{sample}}"
        threads:
            4
        container:
            CONTAINER_4DN
        shell:
            """
            mkdir -p {params.outdir}
            /usr/local/bin/run-add-hicnormvector-to-mcool.sh \
                {input.hic} \
                {input.mcool_raw} \
                {params.outdir}
            mv {params.outdir}/output.mcool {output.mcool}
            rm -rf {params.outdir}
            """
else:
    rule cooler_zoomify:
        input:
            cool = f"{OUTDIR}/cool/{{sample}}.cool"
        output:
            mcool = f"{OUTDIR}/mcool/{{sample}}.mcool"
        params:
            outprefix = f"{OUTDIR}/mcool/{{sample}}",
            resolutions = RESOLUTIONS_STR,
            chunksize = 10000000,
            no_balance_flag = "" if config.get("matrix", {}).get("balance", True) else "-B"
        threads:
            8
        container:
            CONTAINER_4DN
        shell:
            """
            mkdir -p $(dirname {output.mcool})
            /usr/local/bin/run-cool2multirescool.sh \
                -i {input.cool} \
                -o {params.outprefix} \
                -p {threads} \
                -c {params.chunksize} \
                -u {params.resolutions} \
                {params.no_balance_flag}
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
    params:
        outprefix = f"{OUTDIR}/hic/{{sample}}",
        min_res = BASE_RES,
        resolutions = RESOLUTIONS_STR,
        maxmem = lambda wildcards, resources: f"{max(1, int(resources.mem_mb // 1024))}g",
        no_balance_flag = "" if config.get("matrix", {}).get("balance", True) else "-B"
    threads:
        8
    container:
        CONTAINER_4DN
    shell:
        """
        mkdir -p $(dirname {output.hic})
        /usr/local/bin/run-juicebox-pre.sh \
            -i {input.pairs} \
            -c {input.chrom_sizes} \
            -o {params.outprefix} \
            -r {params.min_res} \
            -u {params.resolutions} \
            -m {params.maxmem} \
            {params.no_balance_flag}
        """
