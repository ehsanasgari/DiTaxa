#!/usr/bin/env bash

python3 ditaxa.py --indir dataset/periodontal/ --fast2label dataset/periodontal/mapping.txt --ext fastq --outdir results_dental/ --dbname periodontal --cores 20 --phenomap diseased:1,healthy:0 --phenoname DvsH --override 0 --heatmap periodontal:healthy
