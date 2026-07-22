rule integrate:
    input:
        rds_list=lambda wc: integration_inputs()

    output:
        rds=(
            f"{OUTDIR}/05_integrated/"
            "integrated.rds"
        ),
        umap=(
            f"{OUTDIR}/05_integrated/"
            "umap_overview.pdf"
        )

    params:
        method=config.get(
            "integration_method",
            "Harmony"
        ),

        batch_col=config.get(
            "integration_batch_col",
            "sample"
        ),

        dims=int(
            config.get(
                "pca_dims",
                30
            )
        ),

        npcs=int(
            config.get(
                "npcs",
                50
            )
        ),

        nfeatures=int(
            config.get(
                "nfeatures",
                3000
            )
        ),

        resolution=float(
            config.get(
                "cluster_resolution",
                0.5
            )
        ),

        seed=int(
            config.get(
                "seed",
                1234
            )
        ),

        outdir=f"{OUTDIR}/05_integrated"

    log:
        f"{LOGDIR}/05_integrate/integrate.log"

    benchmark:
        f"{LOGDIR}/benchmarks/05_integrate.txt"

    threads:
        8

    resources:
        mem_mb=64000,
        runtime=720

    script:
        "../scripts/05_integrate.R"