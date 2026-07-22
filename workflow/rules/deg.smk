rule run_deg:
    input:
        rds=(
            f"{OUTDIR}/05_integrated/"
            "integrated.rds"
        ),
        contrasts=lambda wc: config[
            "contrasts_file"
        ]

    output:
        report=(
            f"{OUTDIR}/06_deg/"
            "deg_report.html"
        ),
        deg_dir=directory(
            f"{OUTDIR}/06_deg/tables"
        ),
        plot_dir=directory(
            f"{OUTDIR}/06_deg/plots"
        )

    params:
        method=config.get(
            "deg_method",
            "pseudo_bulk"
        ),

        donor_col=config.get(
            "donor_col",
            "donor_id"
        ),

        fdr_thr=float(
            config.get(
                "fdr_threshold",
                0.05
            )
        ),

        lfc_thr=float(
            config.get(
                "lfc_threshold",
                0.5
            )
        ),

        top_n=int(
            config.get(
                "top_n_genes",
                20
            )
        ),

        run_pathway=as_bool(
            config.get(
                "run_pathway",
                True
            )
        ),

        species=config.get(
            "species",
            "human"
        ),

        outdir=f"{OUTDIR}/06_deg"

    log:
        f"{LOGDIR}/06_deg/deg.log"

    benchmark:
        f"{LOGDIR}/benchmarks/06_deg.txt"

    threads:
        8

    resources:
        mem_mb=48000,
        runtime=720

    script:
        "../scripts/06_run_deg.R"