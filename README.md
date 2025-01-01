# Nextflow pipeline for fastLBP-based image segmentation

## Key points
The pipeline supports the following modes:
- **GridSearch** which allows you to specify lists of parameters for each step to run the pipeline using all possible combinations.
- **SingleImage** to process a single image.
- **MultiImage** to process a folder of images.

## Requirements
**(to be completed)**

## How to use

Choose a template from the `templates/` folder and run
```bash
nextflow run main.nf -profile conda -params-file templates/grid_search_template.yaml
```

The workflow will automatically determine the execution mode based on the template structure.
