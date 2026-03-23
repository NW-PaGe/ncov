#!/bin/bash

echo "Installing NCBI datasets CLI..."
curl -fsSL "https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets" -o datasets
curl -fsSL "https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/dataformat" -o dataformat
chmod +x datasets dataformat


echo "Downloading WA SARS-CoV-2 sequences from NCBI..."
datasets download virus genome taxon sars-cov-2 --usa-state WA --filename ./data/wa_sars2.zip

unzip ./data/wa_sars2.zip -d ./data/wa_sars2

cp ./data/wa_sars2/ncbi_dataset/data/genomic.fna ./data/wa_sars2/ncbi_dataset/data/wa_sequences.fasta

#pull all fields in ncbi dataset
dataformat tsv virus-genome --inputfile ./data/wa_sars2/ncbi_dataset/data/data_report.jsonl > ./data/wa_sars2/ncbi_dataset/data/wa_metadata_raw.tsv

echo "NCBI WA metadata and sequences downloaded"
