#!/bin/bash

datasets download virus genome taxon sars-cov-2 --usa-state WA --filename ./data/wa_sars2.zip

unzip ./data/wa_sars2.zip -d ./data/wa_sars2

cp ./data/wa_sars2/ncbi_dataset/data/genomic.fna ./data/wa_sars2/ncbi_dataset/data/wa_sequences.fasta

#pull all fields in ncbi dataset
dataformat tsv virus-genome --inputfile ./data/wa_sars2/ncbi_dataset/data/data_report.jsonl > ./data/wa_sars2/ncbi_dataset/data/wa_metadata_raw.tsv

echo "NCBI WA metadata and sequences downloaded"
