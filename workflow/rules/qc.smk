rule run_qc:
    input:
        rds=lambda wc: qc_input(wc.sample)

    output:
        rds=(
            f"{OUTDIR}/03_qc/"
            "{sample}_qc.rds"
        ),
        summary=(
            f"{OUTDIR}/03_qc/"
            "QC_summary_{sample}.tsv"
        )

    params:
        sample="{sample}",

        min_feat=lambda wc: int(
            get_meta(
                wc.sample,
                "Min_Feat",
                config.get(
                    "default_min_feat",
                    200
                )
            )
        ),

        min_umi=lambda wc: int(
            get_meta(
                wc.sample,
                "Min_UMI",
                config.get(
                    "default_min_umi",
                    500
                )
            )
        ),

        max_mito=lambda wc: float(
            get_meta(
                wc.sample,
                "Max_Mito",
                config.get(
                    "default_max_mito",
                    5
                )
            )
        ),

        mad_n=lambda wc: int(
            get_meta(
                wc.sample,
                "MAD_N",
                config.get(
                    "default_mad_n",
                    5
                )
            )
        ),

        rm_dbl=lambda wc: as_bool(
            get_meta(
                wc.sample,
                "Rm_Dbl",
                config.get(
                    "default_rm_dbl",
                    False
                )
            )
        ),

        dbl_score=lambda wc: float(
            get_meta(
                wc.sample,
                "Dbl_Score",
                config.get(
                    "default_dbl_score",
                    0.5
                )
            )
        ),

        species=lambda wc: get_meta(
            wc.sample,
            "Species",
            config.get("species", "human")
        ),

        calc_ribo=as_bool(
            config.get("calc_ribo", False)
        ),

        max_ribo=float(
            config.get(
                "default_max_ribo",
                30
            )
        ),

        calc_drop=as_bool(
            config.get("calc_drop", False)
        ),

        max_drop=float(
            config.get(
                "default_max_drop",
                0.95
            )
        ),

        outdir=f"{OUTDIR}/03_qc"

    log:
        f"{LOGDIR}/03_qc/{{sample}}.log"

    benchmark:
        f"{LOGDIR}/benchmarks/03_qc_{{sample}}.txt"

    threads:
        4

    resources:
        mem_mb=32000,
        runtime=240

    script:
        "../scripts/03_run_qc.R"