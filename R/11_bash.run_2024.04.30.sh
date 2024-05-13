#!/bin/bash

source /home/mrd152/.bashrc
cd /data/mrd/scZeroDay/R

# Create output directory if it doesn't exist
output_dir="/data/mrd/scZeroDay/Output"
mkdir -p "$output_dir"

fig_dir="/data/mrd/scZeroDay/Output/Figures/11_spatial"
mkdir -p "$fig_dir"

data_dir="/data/mrd/scZeroDay/Output/Rdata/11_spatial"
mkdir -p "$data_dir"

output_dir_log="/data/mrd/scZeroDay/Output/Rdata/11_spatial/logs"
mkdir -p "$output_dir_log"

# Loop through each Seurat object
for i in {1..18}; do
  export i="$i"
  nohup Rscript 11_spatial.deconvolution_2024.04.30.R "$i" > "$output_dir_log/output_log_$i.txt" 2>&1 # &
  # wait $! # wait is causing it to not run in bkgd
done
