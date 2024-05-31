# Nextlow pipeline for LBP-based image segmentation

## Key points
The pipeline supports the following modes:
- **GridSearch**, where you can specify lists of parameters for each step to run the pipeline using all possible combinations.
- **SingleImage** to process a single image.
- **MultiImage** to process a folder of images.

## How to use

### Find best parameters combinations
To launch in the **GridSearch** mode, use the corresponding template `grid_search_template.yaml` and type:
```bash
nextflow run grid_search.nf -profile conda -params-file templates/grid_search_template.yaml -entry Pipeline
```

### Just process my images

To run in either **SingleImage** or **MultiImage** mode, choose the corresponding template (e.g. `single_image.yaml`) and run:
```bash
nextflow run main.nf -profile conda -params-file templates/single_image.yaml -entry Pipeline
```


Create parameters file based on one of the templates and run: 
```bash
nextflow run main.nf -profile conda -params-file templates/single_image.yaml -entry Pipeline
```
