rule decontaminate:
    input:
        rds=(
            f"{OUTDIR}/03_qc/"
            "{sample}_qc.rds"
        ),
        raw_counts=lambda wc: get_meta(
            wc.sample,
            "CellRanger_Dir"
        )

    output:
        rds=(
            f"{OUTDIR}/04_clean/"
            "{sample}_clean.rds"
        ),
        summary=(
            f"{OUTDIR}/04_clean/"
            "{sample}_decontamination_summary.tsv"
        ),
        per_cell=(
            f"{OUTDIR}/04_clean/"
            "{sample}_contamination_per_cell.tsv"
        )

    params:
        sample="{sample}",

        raw_dir=lambda wc: get_meta(
            wc.sample,
            "CellRanger_Dir"
        ),

        assays=config.get(
            "contam_assays",
            ["RNA", "DecontX"]
        ),

        input_assay=config.get(
            "contam_input_assay",
            "RNA"
        ),

        preferred_assay=config.get(
            "preferred_contam_assay",
            "DecontX"
        ),

        create_clean_alias=as_bool(
            config.get(
                "create_clean_assay",
                True
            )
        ),

        normalize_corrected=as_bool(
            config.get(
                "normalize_corrected_assay",
                False
            )
        ),

        batch_col=config.get(
            "decontx_batch_col",
            ""
        ),

        cluster_col=config.get(
            "decontx_cluster_col",
            "seurat_clusters"
        ),

        outdir=f"{OUTDIR}/04_clean"

    log:
        f"{LOGDIR}/04_contam/{{sample}}.log"

    benchmark:
        f"{LOGDIR}/benchmarks/04_contam_{{sample}}.txt"

    threads:
        4

    resources:
        mem_mb=32000,
        runtime=360

    script:
        "../scripts/04_decontaminate.R"