# =============================================================================
# 05_downstream.smk: Optional High-Resolution Downstream Analyses (cooltools)
# Operates independently on .mcool output
# =============================================================================

EXP_RES = config.get("downstream", {}).get("expected_resolution", 1000)
INS_RES = config.get("downstream", {}).get("insulation_resolution", 10000)
INS_WINDOWS = " ".join(str(w) for w in config.get("downstream", {}).get("insulation_window_bp", [10000, 25000, 50000, 100000]))
COMP_RES = config.get("downstream", {}).get("compartment_resolution", 50000)

# -----------------------------------------------------------------------------
# 1. Contact Probability Decay Curves P(s) vs Distance
# -----------------------------------------------------------------------------
rule cooltools_expected:
    input:
        mcool = f"{OUTDIR}/mcool/{{sample}}.mcool"
    output:
        expected_tsv = f"{OUTDIR}/downstream/expected/{{sample}}_expected_{EXP_RES}bp.tsv"
    params:
        cooler_uri = f"{OUTDIR}/mcool/{{sample}}.mcool::/resolutions/{EXP_RES}"
    threads:
        4
    container:
        CONTAINER_COOLTOOLS
    shell:
        """
        mkdir -p $(dirname {output.expected_tsv})
        cooltools expected-cis \
            -p {threads} \
            -o {output.expected_tsv} \
            {params.cooler_uri}
        """

# -----------------------------------------------------------------------------
# 2. Multi-Scale Insulation Score & Boundary Calling
# -----------------------------------------------------------------------------
rule cooltools_insulation:
    input:
        mcool = f"{OUTDIR}/mcool/{{sample}}.mcool"
    output:
        insulation_tsv = f"{OUTDIR}/downstream/insulation/{{sample}}_insulation_{INS_RES}bp.tsv"
    params:
        cooler_uri = f"{OUTDIR}/mcool/{{sample}}.mcool::/resolutions/{INS_RES}",
        windows = INS_WINDOWS
    threads:
        4
    container:
        CONTAINER_COOLTOOLS
    shell:
        """
        mkdir -p $(dirname {output.insulation_tsv})
        cooltools insulation \
            {params.cooler_uri} \
            {params.windows} \
            --output {output.insulation_tsv}
        """

# -----------------------------------------------------------------------------
# 3. A/B Compartment Eigenvector Decomposition
# -----------------------------------------------------------------------------
rule cooltools_compartments:
    input:
        mcool = f"{OUTDIR}/mcool/{{sample}}.mcool",
        fasta = config.get("reference", {}).get("fasta", "")
    output:
        cis_eigs_tsv = f"{OUTDIR}/downstream/compartments/{{sample}}_eigs_{COMP_RES}bp.cis.vecs.tsv"
    params:
        cooler_uri = f"{OUTDIR}/mcool/{{sample}}.mcool::/resolutions/{COMP_RES}",
        outprefix = f"{OUTDIR}/downstream/compartments/{{sample}}_eigs_{COMP_RES}bp",
        fasta_arg = f"--phasing-track {config.get('reference', {}).get('fasta')}" if (config.get("reference", {}).get("fasta") and os.path.exists(config.get("reference", {}).get("fasta", ""))) else ""
    threads:
        4
    container:
        CONTAINER_COOLTOOLS
    shell:
        """
        mkdir -p $(dirname {output.cis_eigs_tsv})
        cooltools eigs-cis \
            -o {params.outprefix} \
            {params.fasta_arg} \
            {params.cooler_uri}
        """
