#!/bin/bash

mkdir -p ~/Projects/mammal_longevity_transcriptomics/data/geo

cd ~/Projects/mammal_longevity_transcriptomics/data/geo

wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190756/suppl/GSE190756_RAW.tar

tar -tf GSE190756_RAW.tar | head -50