OUTDIR = "out"
SMKDIR = f"{OUTDIR}/smk"

# Foldseek binary. Override with `snakemake --config foldseek=/path/to/foldseek`.
FOLDSEEK = config.get("foldseek", "foldseek")

# Query subset size. Override with `snakemake --config sub=2` for quick testing.
SUB = int(config.get("sub", 2000))
TAG = f"sub{SUB}"

# Target databases. The key names the local db under {SMKDIR}/ and is the
# `_{db}` suffix on all search/gscore/evalue outputs; the value is the
# Foldseek database name passed to `foldseek databases`.
DB_NAMES = {
    "afdb50":    "Alphafold/UniProt50",
    "swissprot": "Alphafold/Swiss-Prot",
}

wildcard_constraints:
    db   = "|".join(DB_NAMES),
    name = "[A-Za-z0-9]+",

rule all:
    input:
        f"{SMKDIR}/plot_gscore_firsthit_{TAG}_swissprot.png",
        f"{SMKDIR}/plot_tmscore_firsthit_{TAG}_swissprot.png",
        f"{SMKDIR}/evalue_firsthit_fractions_{TAG}_swissprot.txt",

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
    output: f"{SMKDIR}/FESNov_families.pdb.{TAG}.tar.gz"
    shell:  f"python scripts/make_subset_tar.py {{input.tgz}} {{input.plddt}} {{output}} {SUB}"

# --mask-bfactor-threshold 50 lower-cases positions with pLDDT<50 so gscore
# is computed only over confident positions. Prefilter must pass
# --mask-lower-case 0 so masking does not affect prefiltering.
rule createdb_query:
    input:  f"{SMKDIR}/FESNov_families.pdb.tar.gz"
    output: f"{SMKDIR}/qDB_all"
    shell:  f"{FOLDSEEK} createdb --ss-12st 1 --mask-bfactor-threshold 50 {{input}} {{output}}"

rule createdb_query_subset:
    input:  f"{SMKDIR}/FESNov_families.pdb.{TAG}.tar.gz"
    output: f"{SMKDIR}/qDB_{TAG}"
    shell:  f"{FOLDSEEK} createdb --ss-12st 1 --mask-bfactor-threshold 50 {{input}} {{output}}"

rule download_target:
    output: f"{SMKDIR}/{{db}}"
    threads: 32
    params:
        fs_name = lambda w: DB_NAMES[w.db],
    shell:
        f"{FOLDSEEK} databases {{params.fs_name}} {{output}} {SMKDIR}/tmp_{{wildcards.db}} && "
        f"{FOLDSEEK} add12st --threads {{threads}} {{output}}"

rule search_foldseek1:
    input:
        qdb = f"{SMKDIR}/qDB_{TAG}",
        tdb = f"{SMKDIR}/{{db}}",
    output: f"{SMKDIR}/search_foldseek1_{TAG}_{{db}}.tsv"
    threads: 32
    shell:
        f"""
        mkdir -p {{output}}_tmp
        {FOLDSEEK} prefilter {{input.qdb}}_ss {{input.tdb}}_ss {{output}}_tmp/prefDB \
            -s 9.5 -k 6 --max-seqs 2000 --comp-bias-corr 1 --comp-bias-corr-scale 0.15 \
            --mask-lower-case 0 --aux-score 0 --threads {{threads}}
        {FOLDSEEK} structurealign {{input.qdb}} {{input.tdb}} {{output}}_tmp/prefDB {{output}}_tmp/alnDB \
            -e 10000 --sort-by-structure-bits 0 --ss-12st 0 -a --threads {{threads}}
        {FOLDSEEK} convertalis {{input.qdb}} {{input.tdb}} {{output}}_tmp/alnDB {{output}}_tmp/full.tsv \
            --format-output 'query,target' --threads {{threads}}
        awk -F'\\t' '!seen[$1]++' {{output}}_tmp/full.tsv > {{output}}
        rm -rf {{output}}_tmp
        """

rule search_foldseek1_pred:
    input:
        pred_qdb = f"{OUTDIR}/pred_db/{{name}}/db",
        tdb      = f"{SMKDIR}/{{db}}",
    output: f"{SMKDIR}/search_foldseek1_{TAG}_pred_{{name}}_{{db}}.tsv"
    threads: 32
    shell:
        f"""
        mkdir -p {{output}}_tmp
        {FOLDSEEK} prefilter {{input.pred_qdb}}_ss {{input.tdb}}_ss {{output}}_tmp/prefDB \
            -s 9.5 -k 6 --max-seqs 2000 --comp-bias-corr 1 --comp-bias-corr-scale 0.15 \
            --mask-lower-case 0 --aux-score 0 --threads {{threads}}
        {FOLDSEEK} structurealign {{input.pred_qdb}} {{input.tdb}} {{output}}_tmp/prefDB {{output}}_tmp/alnDB \
            -e 10000 --sort-by-structure-bits 0 --ss-12st 0 -a --threads {{threads}}
        {FOLDSEEK} convertalis {{input.pred_qdb}} {{input.tdb}} {{output}}_tmp/alnDB {{output}}_tmp/full.tsv \
            --format-output 'query,target' --threads {{threads}}
        awk -F'\\t' '!seen[$1]++' {{output}}_tmp/full.tsv > {{output}}
        rm -rf {{output}}_tmp
        """

rule search_foldseek2:
    input:
        qdb = f"{SMKDIR}/qDB_{TAG}",
        tdb = f"{SMKDIR}/{{db}}",
    output: f"{SMKDIR}/search_foldseek2_{TAG}_{{db}}.tsv"
    threads: 32
    shell:
        f"""
        mkdir -p {{output}}_tmp
        {FOLDSEEK} prefilter {{input.qdb}}_ss {{input.tdb}}_ss {{output}}_tmp/prefDB \
            -s 9.5 -k 6 --max-seqs 2000 --comp-bias-corr 1 --comp-bias-corr-scale 0.15 \
            --mask-lower-case 0 --aux-score 1 --threads {{threads}}
        {FOLDSEEK} structurealign {{input.qdb}} {{input.tdb}} {{output}}_tmp/prefDB {{output}}_tmp/alnDB \
            -e 10000 --sort-by-structure-bits 0 --ss-12st 1 \
            --use-reverse-score 0 \
            --gap-open aa:14,nucl:14 --gap-extend aa:2,nucl:2 \
            -a --threads {{threads}}
        {FOLDSEEK} convertalis {{input.qdb}} {{input.tdb}} {{output}}_tmp/alnDB {{output}}_tmp/full.tsv \
            --format-output 'query,target' --threads {{threads}}
        awk -F'\\t' '!seen[$1]++' {{output}}_tmp/full.tsv > {{output}}
        rm -rf {{output}}_tmp
        """

rule search_foldseek2_pred:
    input:
        pred_qdb = f"{OUTDIR}/pred_db/{{name}}/db",
        tdb      = f"{SMKDIR}/{{db}}",
    output: f"{SMKDIR}/search_foldseek2_{TAG}_pred_{{name}}_{{db}}.tsv"
    threads: 32
    shell:
        f"""
        mkdir -p {{output}}_tmp
        {FOLDSEEK} prefilter {{input.pred_qdb}}_ss {{input.tdb}}_ss {{output}}_tmp/prefDB \
            -s 9.5 -k 6 --max-seqs 2000 --comp-bias-corr 1 --comp-bias-corr-scale 0.15 \
            --mask-lower-case 0 --aux-score 1 --threads {{threads}}
        {FOLDSEEK} structurealign {{input.pred_qdb}} {{input.tdb}} {{output}}_tmp/prefDB {{output}}_tmp/alnDB \
            -e 10000 --sort-by-structure-bits 0 --ss-12st 1 \
            --use-reverse-score 0 \
            --gap-open aa:14,nucl:14 --gap-extend aa:2,nucl:2 \
            -a --threads {{threads}}
        {FOLDSEEK} convertalis {{input.pred_qdb}} {{input.tdb}} {{output}}_tmp/alnDB {{output}}_tmp/full.tsv \
            --format-output 'query,target' --threads {{threads}}
        awk -F'\\t' '!seen[$1]++' {{output}}_tmp/full.tsv > {{output}}
        rm -rf {{output}}_tmp
        """

rule search_foldseek2_pred_structbit1:
    input:
        pred_qdb = f"{OUTDIR}/pred_db/{{name}}/db",
        tdb      = f"{SMKDIR}/{{db}}",
    output: f"{SMKDIR}/search_foldseek2_{TAG}_pred_{{name}}_structbit1_{{db}}.tsv"
    threads: 32
    shell:
        f"""
        mkdir -p {{output}}_tmp
        {FOLDSEEK} prefilter {{input.pred_qdb}}_ss {{input.tdb}}_ss {{output}}_tmp/prefDB \
            -s 9.5 -k 6 --max-seqs 2000 --comp-bias-corr 1 --comp-bias-corr-scale 0.15 \
            --mask-lower-case 0 --aux-score 1 --threads {{threads}}
        {FOLDSEEK} structurealign {{input.pred_qdb}} {{input.tdb}} {{output}}_tmp/prefDB {{output}}_tmp/alnDB \
            -e 10000 --sort-by-structure-bits 1 --ss-12st 1 \
            --use-reverse-score 0 \
            --gap-open aa:14,nucl:14 --gap-extend aa:2,nucl:2 \
            -a --threads {{threads}}
        {FOLDSEEK} convertalis {{input.pred_qdb}} {{input.tdb}} {{output}}_tmp/alnDB {{output}}_tmp/full.tsv \
            --format-output 'query,target' --threads {{threads}}
        awk -F'\\t' '!seen[$1]++' {{output}}_tmp/full.tsv > {{output}}
        rm -rf {{output}}_tmp
        """

rule search_mmseqs:
    input:
        qdb = f"{SMKDIR}/qDB_{TAG}",
        tdb = f"{SMKDIR}/{{db}}",
    output: f"{SMKDIR}/search_mmseqs_{TAG}_{{db}}.tsv"
    threads: 32
    shell:
        r"""
        mkdir -p {output}_tmp
        mmseqs search {input.qdb} {input.tdb} {output}_tmp/resDB {output}_tmp/wd \
            -a -s 7.5 -e 10000 --max-seqs 2000 --threads {threads}
        mmseqs convertalis {input.qdb} {input.tdb} {output}_tmp/resDB {output}_tmp/full.tsv \
            --format-output 'query,target' --threads {threads}
        awk -F'\t' '!seen[$1]++' {output}_tmp/full.tsv > {output}
        rm -rf {output}_tmp
        """

# Unified gscore step: rebuild a prefDB keyed to qDB_{TAG} numeric IDs from the
# (qname, tname) tsv, then run lolalign+convertalis so gscore is computed the
# same way for every method. lolalign gives higher-quality alignments than
# structurealign, raising the gscore distribution.
rule gscore_tophit:
    input:
        tsv = f"{SMKDIR}/search_{{stem}}_{{db}}.tsv",
        qdb = f"{SMKDIR}/qDB_{TAG}",
        tdb = f"{SMKDIR}/{{db}}",
    output: f"{SMKDIR}/gscore_{{stem}}_{{db}}.tsv"
    threads: 32
    shell:
        f"""
        mkdir -p {{output}}_tmp
        python scripts/compare_dbs.py build_prefdb \
            {{input.tsv}} {{input.qdb}} {{input.tdb}} {{output}}_tmp/prefDB
        {FOLDSEEK} lolalign {{input.qdb}} {{input.tdb}} {{output}}_tmp/prefDB {{output}}_tmp/alnDB \
            -a --threads {{threads}}
        {FOLDSEEK} convertalis {{input.qdb}} {{input.tdb}} {{output}}_tmp/alnDB {{output}} \
            --format-output 'query,target,evalue,bits,qstart,tstart,cigar,gscore,alntmscore' \
            --threads {{threads}}
        rm -rf {{output}}_tmp
        """

# Like gscore_tophit, but computes the Foldseek1-style evalue
# (structurealign with --ss-12st 0, default gap params, no reverse-score).
# Note: evalues are always computed against the ground-truth qDB_{TAG}
# sequences, regardless of which {stem} produced the top hits.
rule evalue_fs1_tophit:
    input:
        tsv = f"{SMKDIR}/search_{{stem}}_{{db}}.tsv",
        qdb = f"{SMKDIR}/qDB_{TAG}",
        tdb = f"{SMKDIR}/{{db}}",
    output: f"{SMKDIR}/evalue_{{stem}}_{{db}}.tsv"
    threads: 32
    shell:
        f"""
        mkdir -p {{output}}_tmp
        python scripts/compare_dbs.py build_prefdb \
            {{input.tsv}} {{input.qdb}} {{input.tdb}} {{output}}_tmp/prefDB
        {FOLDSEEK} structurealign {{input.qdb}} {{input.tdb}} {{output}}_tmp/prefDB {{output}}_tmp/alnDB \
            -e 10000 --sort-by-structure-bits 0 --ss-12st 0 \
            -a --threads {{threads}}
        {FOLDSEEK} convertalis {{input.qdb}} {{input.tdb}} {{output}}_tmp/alnDB {{output}} \
            --format-output 'query,target,evalue,bits,qstart,tstart,cigar' \
            --threads {{threads}}
        rm -rf {{output}}_tmp
        """

rule plot_gscore_firsthit:
    input:
        fs1      = f"{SMKDIR}/gscore_foldseek1_{TAG}_{{db}}.tsv",
        fs1_pred = f"{SMKDIR}/gscore_foldseek1_{TAG}_pred_prostt5_{{db}}.tsv",
        fs2      = f"{SMKDIR}/gscore_foldseek2_{TAG}_{{db}}.tsv",
        fs2_pred = f"{SMKDIR}/gscore_foldseek2_{TAG}_pred_mprost12B_{{db}}.tsv",
        fs2_pred_sb1 = f"{SMKDIR}/gscore_foldseek2_{TAG}_pred_mprost1Bvqvae_structbit1_{{db}}.tsv",
        mmseqs   = f"{SMKDIR}/gscore_mmseqs_{TAG}_{{db}}.tsv",
        qdb      = f"{SMKDIR}/qDB_{TAG}",
    output:
        gscore  = f"{SMKDIR}/plot_gscore_firsthit_{TAG}_{{db}}.png",
        tmscore = f"{SMKDIR}/plot_tmscore_firsthit_{TAG}_{{db}}.png",
    shell:
        "python scripts/plot_gscore_firsthit.py --col 7 --ylabel 'gscore of top hit' "
        "--queries {input.qdb}.lookup {output.gscore} "
        "{input.fs1} {input.fs1_pred} {input.fs2} {input.fs2_pred} {input.fs2_pred_sb1} {input.mmseqs} && "
        "python scripts/plot_gscore_firsthit.py --col 8 --ylabel 'TM-score of top hit' "
        "--queries {input.qdb}.lookup {output.tmscore} "
        "{input.fs1} {input.fs1_pred} {input.fs2} {input.fs2_pred} {input.fs2_pred_sb1} {input.mmseqs}"

rule print_evalue_firsthit:
    input:
        fs1      = f"{SMKDIR}/evalue_foldseek1_{TAG}_{{db}}.tsv",
        fs1_pred = f"{SMKDIR}/evalue_foldseek1_{TAG}_pred_prostt5_{{db}}.tsv",
        fs2      = f"{SMKDIR}/evalue_foldseek2_{TAG}_{{db}}.tsv",
        fs2_pred = f"{SMKDIR}/evalue_foldseek2_{TAG}_pred_mprost12B_{{db}}.tsv",
        fs2_pred_sb1 = f"{SMKDIR}/evalue_foldseek2_{TAG}_pred_mprost1Bvqvae_structbit1_{{db}}.tsv",
        mmseqs   = f"{SMKDIR}/evalue_mmseqs_{TAG}_{{db}}.tsv",
        qdb      = f"{SMKDIR}/qDB_{TAG}",
    output: f"{SMKDIR}/evalue_firsthit_fractions_{TAG}_{{db}}.txt"
    shell:
        r"""
        N=$(wc -l < {input.qdb}.lookup)
        awk -F'\t' -v n=$N '$3<=1e-3{{c++}} END{{printf "%s\t%d/%d\t%.4f\n", FILENAME, c, n, c/n}}' {input.fs1}          >  {output}
        awk -F'\t' -v n=$N '$3<=1e-3{{c++}} END{{printf "%s\t%d/%d\t%.4f\n", FILENAME, c, n, c/n}}' {input.fs1_pred}     >> {output}
        awk -F'\t' -v n=$N '$3<=1e-3{{c++}} END{{printf "%s\t%d/%d\t%.4f\n", FILENAME, c, n, c/n}}' {input.fs2}          >> {output}
        awk -F'\t' -v n=$N '$3<=1e-3{{c++}} END{{printf "%s\t%d/%d\t%.4f\n", FILENAME, c, n, c/n}}' {input.fs2_pred}     >> {output}
        awk -F'\t' -v n=$N '$3<=1e-3{{c++}} END{{printf "%s\t%d/%d\t%.4f\n", FILENAME, c, n, c/n}}' {input.fs2_pred_sb1} >> {output}
        awk -F'\t' -v n=$N '$3<=1e-3{{c++}} END{{printf "%s\t%d/%d\t%.4f\n", FILENAME, c, n, c/n}}' {input.mmseqs}       >> {output}
        cat {output}
        """
