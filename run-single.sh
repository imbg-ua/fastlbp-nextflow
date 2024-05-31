#!/usr/bin/env bash
nextflow run main.nf -profile conda -params-file templates/single_image.yaml -entry Pipeline