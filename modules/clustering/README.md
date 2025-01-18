# Run only Clustering module

Use an existing conda environment to avoid creating it again.
```bash
conda activate workflow_environment
nextflow run main.nf --outdir /path/to/outdir --inputs_tsv clustering_module_params.tsv
```
