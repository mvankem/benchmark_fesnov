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
3. (optional) All data is written under `out/` (downloads, databases, search results;
   tens to hundreds of GB). To put it elsewhere, symlink before running:
   ```
   ln -s /scratch/$USER/fesnov out
   ```
4. Foldseek is expected on `$PATH`. To install a specific version into the
   active conda env, or override the binary used:
   ```
   ./scripts/install_foldseek <commit-hash | tag | branch>   # default: master
   # or: snakemake --config foldseek=/path/to/foldseek ...
   ```
5. Run the pipeline (search + gscore + plot):
   ```
   snakemake -c 32 out/smk/plot_gscore_firsthit_sub2000.png
   ```

## Notes

snakemake --cores 32 out/smk/search_foldseek1_sub2000.tsv
awk '!seen[$1]++' out/smk/search_foldseek1_sub2000.tsv | awk '$3 < 0.001' | wc -l

sbatch -p soeding run.sbatch out/smk/plot_gscore_firsthit_sub2000.png
sbatch -p soeding run.sbatch --config sub=2

<!-- TMP (remove before final): instead of the createdb masking hack
(--mask-bfactor-threshold 50) we drive gscore via foldseek, we could
emit the per-hit pdbs via --format-output and compute gscore ourselves
from CA coords + the alignment. -->
