OUTDIR = "out"
SMKDIR = f"{OUTDIR}/smk"
AFDB50 = f"{SMKDIR}/afdb50"
#AFDB50 = "/cbscratch/michel/prefilter_benchmark2/martin_db/tDB"

rule all:
    input:
        f"{SMKDIR}/plot_gscore_firsthit_sub2000.png",

rule download_query:
    output:
        tgz   = f"{SMKDIR}/FESNov_families.pdb.tar.gz",
        plddt = f"{SMKDIR}/pLDDT_values.tab",
    shell:
        "curl -L -o {output.tgz} "
        "  'https://zenodo.org/records/10242439/files/FESNov_families.pdb.tar.gz?download=1' && "
        "curl -L -o {output.plddt} "
        "  'https://zenodo.org/records/10242439/files/pLDDT_values.tab?download=1'"

rule subset_query:
    input:
        tgz   = f"{SMKDIR}/FESNov_families.pdb.tar.gz",
        plddt = f"{SMKDIR}/pLDDT_values.tab",
    output: f"{SMKDIR}/FESNov_families.pdb.sub2000.tar.gz"
    shell:  "python scripts/make_subset_tar.py {input.tgz} {input.plddt} {output}"

# --mask-bfactor-threshold 50 lower-cases positions with pLDDT<50 so gscore
# is computed only over confident positions. Prefilter must pass
# --mask-lower-case 0 so masking does not affect prefiltering.
rule createdb_query:
    input:  f"{SMKDIR}/FESNov_families.pdb.tar.gz"
    output: f"{SMKDIR}/qDB_all"
    shell:  "foldseek createdb --ss-12st 1 --mask-bfactor-threshold 50 {input} {output}"

rule createdb_query_subset:
    input:  f"{SMKDIR}/FESNov_families.pdb.sub2000.tar.gz"
    output: f"{SMKDIR}/qDB_sub2000"
    shell:  "foldseek createdb --ss-12st 1 --mask-bfactor-threshold 50 {input} {output}"

rule download_target:
    output: AFDB50
    threads: 32
    shell:
        f"foldseek databases Alphafold/UniProt50 {{output}} {SMKDIR}/tmp && "
        "foldseek add12st --threads {threads} {output}"

rule search_foldseek1:
    input:
        qdb = f"{SMKDIR}/qDB_{{tag}}",
        tdb = AFDB50,
    output: f"{SMKDIR}/search_foldseek1_{{tag}}.tsv"
    threads: 32
    shell:
        """# (This could be also done with one easy-search command.)
        mkdir -p {output}_tmp
        foldseek prefilter {input.qdb}_ss {input.tdb}_ss {output}_tmp/prefDB \
            -s 9.5 -k 6 --max-seqs 2000 --comp-bias-corr 1 --comp-bias-corr-scale 0.15 \
            --mask-lower-case 0 --aux-score 0 --threads {threads}
        foldseek structurealign {input.qdb} {input.tdb} {output}_tmp/prefDB {output}_tmp/alnDB \
            -e 10000 --sort-by-structure-bits 0 --ss-12st 0 -a --threads {threads}
        # TODO: reindex a query db with the full structure based on {input.qdb} -> querydb_ca
        # and use this in convertalis to calculate gscores
        foldseek convertalis {input.qdb} {input.tdb} {output}_tmp/alnDB {output} \
            --format-output 'query,target,evalue,bits,qstart,tstart,cigar,gscore' \
            --threads {threads}
        rm -rf {output}_tmp
        """

rule search_foldseek2_pred:
    input:
        pred_qdb = f"{OUTDIR}/pred_db/{{name}}/db",
        gt_qdb   = f"{SMKDIR}/qDB_{{tag}}",
        tdb      = AFDB50,
    output: f"{SMKDIR}/search_foldseek2_{{tag}}_pred_{{name}}.tsv"
    threads: 32
    shell:
        """
        mkdir -p {output}_tmp
        python scripts/compare_dbs.py reindex \
            {input.pred_qdb} {input.gt_qdb} {output}_tmp/gtReidx
        foldseek prefilter {input.pred_qdb}_ss {input.tdb}_ss {output}_tmp/prefDB \
            -s 9.5 -k 6 --max-seqs 2000 --comp-bias-corr 1 --comp-bias-corr-scale 0.15 \
            --mask-lower-case 0 --aux-score 1 --threads {threads}
        # pred-side alignment + tsv (gives the ranking / top hit per query)
        foldseek structurealign {input.pred_qdb} {input.tdb} {output}_tmp/prefDB {output}_tmp/predAlnDB \
            -e 10000 --sort-by-structure-bits 0 --ss-12st 1 \
            --use-reverse-score 0 \
            --gap-open aa:14,nucl:14 --gap-extend aa:2,nucl:2 \
            -a --threads {threads}
        # pred has no _ca, so we cannot compute gscore here; only need query/target for the merge.
        foldseek convertalis {input.pred_qdb} {input.tdb} {output}_tmp/predAlnDB {output}_tmp/pred.tsv \
            --format-output 'query,target,bits' \
            --threads {threads}
        # GT-side alignment + tsv (gives the per-pair score on ground-truth coords).
        # -e 10000 so pred's top target is rarely filtered out before the merge.
        foldseek structurealign {output}_tmp/gtReidx {input.tdb} {output}_tmp/prefDB {output}_tmp/gtAlnDB \
            -e 10000 --sort-by-structure-bits 0 --ss-12st 1 \
            --use-reverse-score 0 \
            --gap-open aa:14,nucl:14 --gap-extend aa:2,nucl:2 \
            -a --threads {threads}
        foldseek convertalis {output}_tmp/gtReidx {input.tdb} {output}_tmp/gtAlnDB {output}_tmp/gt.tsv \
            --format-output 'query,target,evalue,bits,qstart,tstart,cigar,gscore' \
            --threads {threads}
        # For each query take the GT row whose target matches pred's top target.
        awk -F'\\t' 'NR==FNR{{if(!seen[$1]++) top[$1]=$2; next}} ($1 in top) && top[$1]==$2' \
            {output}_tmp/pred.tsv {output}_tmp/gt.tsv > {output}
        rm -rf {output}_tmp
        """

rule plot_gscore_firsthit:
    input:
        fs1  = f"{SMKDIR}/search_foldseek1_{{tag}}.tsv",
        fs2  = f"{SMKDIR}/search_foldseek2_{{tag}}.tsv",
        pred = f"{SMKDIR}/search_foldseek2_{{tag}}_pred_mprost12B.tsv",
        qdb  = f"{SMKDIR}/qDB_{{tag}}",
    output: f"{SMKDIR}/plot_gscore_firsthit_{{tag}}.png"
    shell:
        "python scripts/plot_gscore_firsthit.py "
        "--queries {input.qdb}.lookup {output} "
        "{input.fs1} {input.fs2} {input.pred}"

rule search_foldseek2:
    input:
        qdb = f"{SMKDIR}/qDB_{{tag}}",
        tdb = AFDB50,
    output: f"{SMKDIR}/search_foldseek2_{{tag}}.tsv"
    threads: 32
    shell:
        """# (This could be also done with one easy-search command.)
        mkdir -p {output}_tmp
        foldseek prefilter {input.qdb}_ss {input.tdb}_ss {output}_tmp/prefDB \
            -s 9.5 -k 6 --max-seqs 2000 --comp-bias-corr 1 --comp-bias-corr-scale 0.15 \
            --mask-lower-case 0 --aux-score 1 --threads {threads}
        foldseek structurealign {input.qdb} {input.tdb} {output}_tmp/prefDB {output}_tmp/alnDB \
            -e 10000 --sort-by-structure-bits 0 --ss-12st 1 \
            --use-reverse-score 0 \
            --gap-open aa:14,nucl:14 --gap-extend aa:2,nucl:2 \
            -a --threads {threads}
        # TODO: reindex a query db with the full structure based on {input.qdb} -> querydb_ca
        # and use this in convertalis to calculate gscores
        foldseek convertalis {input.qdb} {input.tdb} {output}_tmp/alnDB {output} \
            --format-output 'query,target,evalue,bits,qstart,tstart,cigar,gscore' \
            --threads {threads}
        rm -rf {output}_tmp
        """
