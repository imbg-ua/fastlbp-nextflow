# Nextflow pipeline for fastLBP-based image segmentation

## Key points
The pipeline supports the following modes:
- **GridSearch** which allows you to specify lists of parameters for each step to run the pipeline using all possible combinations.
- **SingleImage** to process a single image.
- **MultiImage** to process a folder of images.

## Requirements
**(to be completed)**

## How to use

Choose a `grid` mode template from the `templates/` folder and run:
```bash
nextflow run main.nf -profile conda --mode grid -params-file templates/grid_search_template.yaml
```

To run the workflow in `Single Image` or `Multi Image` mode, change the `--mode` flag and provide the appropriate template file:
```bash
nextflow run main.nf -profile conda --mode normal -params-file templates/multiple_images.yaml
```

The workflow will automatically determine the execution mode based on the template structure.

### Run only fastLBP

To process images with fastLBP only, check out the structure required for the `TSV` file with input parameters and run:
```bash
nextflow run main.nf -profile conda --mode lbp_only --lbp_runs_tsv modules/feature_extraction/lbp_only_template.tsv
```
