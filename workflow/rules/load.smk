rule load_seurat:
    input:
        cellranger=lambda wc: get_meta(
            wc.sample,
            "CellRanger_Dir"
        )

    output:
        rds=(
            f"{OUTDIR}/01_loaded/"
            "{sample}_raw.rds"
        )

    params:
        sample="{sample}",
        species=lambda wc: get_meta(
            wc.sample,
            "Species",
            config.get("species", "human")
        ),
        use_filt=as_bool(
            config.get("use_filtered", True),
            default=True
        ),
        use_h5=as_bool(
            config.get("use_h5", True),
            default=True
        ),
        outdir=f"{OUTDIR}/01_loaded"

    log:
        f"{LOGDIR}/01_load/{{sample}}.log"

    benchmark:
        f"{LOGDIR}/benchmarks/01_load_{{sample}}.txt"

    threads:
        4

    resources:
        mem_mb=16000,
        runtime=120

    script:
        "../scripts/01_load_seurat.R"