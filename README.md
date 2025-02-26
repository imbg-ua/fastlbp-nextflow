# Nextflow pipeline for fastLBP-based image segmentation

## Key points
The pipeline supports the following modes:
- **GridSearch** which allows you to specify lists of parameters for each step to run the pipeline using all possible combinations.
- **SingleImage** to process a single image.
- **MultiImage** to process a set of images.

## Requirements
**(to be completed)**

## How to use

The `templates/` folder provides config file templates to run the workflow in different modes. Combined with the `--mode` flag, it automatically determines the execution mode based on the template structure.

## Examples

### Grid mode
```bash
nextflow run main.nf -profile conda --mode grid -params-file templates/grid_search_template.yaml
```

### Process a single image

```bash
nextflow run main.nf -profile conda --mode normal -params-file templates/single_image.yaml
```

### Process multiple images

```bash
nextflow run main.nf -profile conda --mode normal -params-file templates/multiple_images.yaml
```

### fastLBP Only

You can run only `fastLBP` instead of the full analysis pipeline in both Single and Multi Image modes:

```bash
nextflow run main.nf -profile conda --mode lbp_only -params-file templates/multiple_images.yaml
```

```bash
nextflow run main.nf -profile conda --mode lbp_only -params-file templates/single_image.yaml
```

Alternatively, if you want a more flexible option than using the same set of LBP parameters provided in a YAML file to process all images of interest, you can opt for the `lbp_tsv` mode, which accepts a `TSV` with each run's info on a new line.

```bash
nextflow run main.nf -profile conda --mode lbp_tsv --lbp_runs_tsv modules/feature_extraction/lbp_only_template.tsv --outdir /path/to/outdir
```


### LSF
Use LSF configuration file to run pipeline on an LSF cluster.
```bash
nextflow run main.nf -profile conda -params-file templates/multiple_images.yaml -c lsf_config.config --mode lbp_only
```
