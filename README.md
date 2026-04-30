# benchmark_fesnov

Benchmark protein-structure search methods on the FESNov families dataset
(Nature 2023, [s41586-023-06955-z](https://www.nature.com/articles/s41586-023-06955-z)).

## Reproduce

1. Create conda env (once).
   ```
   conda env create -f envs/env.yaml
   ```
2. Activate the env:
   ```
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

## Notes

snakemake --cores 32 out/smk/search_foldseek1_sub2000.tsv
awk '!seen[$1]++' out/smk/search_foldseek1_sub2000.tsv | awk '$3 < 0.001' | wc -l
