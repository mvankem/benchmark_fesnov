# FESNov protein-structure search benchmark.
#
# Workflow:
#   snakemake out/smk/results/<label>.tsv --cores N    # one per install
#   snakemake comparison_plot.png --cores N            # plot all results
#
# By default looks for `foldseek` on $PATH. Override the path with:
#   --config foldseek=/path/to/foldseek

import glob

OUTDIR   = "out/smk"
DATA_URL = "https://zenodo.org/records/10242439/files/FESNov_families.pdb.tar.gz"
FOLDSEEK = config.get("foldseek", "foldseek")

DATA_TGZ = f"{OUTDIR}/data/FESNov_families.pdb.tar.gz"
DATA_DIR = f"{OUTDIR}/data/FESNov_families.pdb"
RESULTS  = f"{OUTDIR}/results"
PLOT     = f"{OUTDIR}/comparison_plot.png"

wildcard_constraints:
    label = r"[A-Za-z0-9_.\-]+"

rule all:
    input: PLOT

rule download_dataset:
    output: DATA_TGZ
    shell:
        "mkdir -p $(dirname {output}) && "
        "curl -L -o {output} {DATA_URL!r}"

rule extract_dataset:
    input:  DATA_TGZ
    output: directory(DATA_DIR)
    shell:
        "mkdir -p {output} && "
        "tar -xzf {input} -C $(dirname {output})"

rule search:
    input:  DATA_DIR
    output: f"{RESULTS}/{{label}}.tsv"
    threads: 8
    shell:
        "tmp=$(mktemp -d) && "
        "{FOLDSEEK} easy-search "
        "  {input} {input} {output} $tmp "
        "  --threads {threads} "
        "  --format-output query,target,evalue,bits && "
        "rm -rf $tmp"

rule comparison_plot:
    input:  lambda wc: sorted(glob.glob(f"{RESULTS}/*.tsv"))
    output: PLOT
    shell:  "python scripts/plot.py --output {output} {input}"
