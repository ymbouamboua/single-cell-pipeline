rule add_demux:
    input:
        rds=(
            f"{OUTDIR}/01_loaded/"
            "{sample}_raw.rds"
        )

    output:
        rds=(
            f"{OUTDIR}/02_demuxed/"
            "{sample}_demuxed.rds"
        ),
        summary=(
            f"{OUTDIR}/02_demuxed/"
            "{sample}_demux_summary.tsv"
        ),
        assignments=(
            f"{OUTDIR}/02_demuxed/"
            "{sample}_assignment_counts.tsv"
        )

    params:
        sample="{sample}",
        demuxafy_dir=lambda wc: get_meta(
            wc.sample,
            "Demuxafy_Dir",
            ""
        ),
        bulk_tsv=lambda wc: get_meta(
            wc.sample,
            "BulkMapping_TSV",
            ""
        ),
        remove_doublets=as_bool(
            config.get(
                "remove_demux_doublets",
                False
            )
        ),
        remove_unassigned=as_bool(
            config.get(
                "remove_demux_unassigned",
                True
            )
        ),
        assignment_preference=config.get(
            "demux_assignment_preference",
            "bulk"
        ),
        outdir=f"{OUTDIR}/02_demuxed"

    wildcard_constraints:
        sample="|".join(
            map(str, DEMUX_SAMPLES)
        ) if DEMUX_SAMPLES else r"(?!x)x"

    log:
        f"{LOGDIR}/02_demux/{{sample}}.log"

    benchmark:
        f"{LOGDIR}/benchmarks/02_demux_{{sample}}.txt"

    threads:
        2

    resources:
        mem_mb=8000,
        runtime=120

    script:
        "../scripts/02_add_demux.R"