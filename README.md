# benchmark_fesnov

Benchmark protein-structure search methods on the FESNov families dataset
(Nature 2023, [s41586-023-06955-z](https://www.nature.com/articles/s41586-023-06955-z)).

## Reproduce

1. Install conda (once).
2. Create + activate the env:
   ```
   conda env create -f envs/env.yaml
   conda activate benchmark-fesnov
   ```
3. By default the pipeline looks for `foldseek` on `$PATH`. Override the
   path with `--config foldseek=/path/to/foldseek` on the snakemake
   commands below. To install one into the active conda env (no
   system-wide changes), use:
   ```
   ./scripts/install_foldseek <commit-hash | tag | branch>   # default: master
   ```
4. Run a search and label it (one tsv per install you want to compare):
   ```
   snakemake out/smk/results/<label>.tsv --cores 8
   ```
5. Plot whatever results have accumulated:
   ```
   snakemake comparison_plot.png --cores 1
   ```

## Layout

```
Snakefile           pipeline (config inline; OUTDIR=./out/smk)
envs/env.yaml       conda env (snakemake + matplotlib)
scripts/
  install_foldseek  builds foldseek at a git ref into $CONDA_PREFIX
  plot.py           placeholder comparison plot
```

## Requirements

`install_foldseek` calls system `cmake`, `git`, and a C++ compiler. On HPC
systems with `module`, it runs `module load rust`. Anything else needed for
the foldseek build must be available in the environment already.
