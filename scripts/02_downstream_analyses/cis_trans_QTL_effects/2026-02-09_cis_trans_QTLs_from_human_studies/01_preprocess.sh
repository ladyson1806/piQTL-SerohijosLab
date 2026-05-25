#!/bin/sh

cd $(dirname ${0})

mkdir -p data
mkdir -p out/formatted_tables
mkdir -p out/figures
mkdir -p out/t_test

# Unzip Vosa_cis_eQTLs.tsv.gz to data/Vosa_cis_eQTLs.tsv
if [ ! -f data/Vosa_cis_eQTLs.tsv ]; then
    gunzip -c data/Vosa_cis_eQTLs.tsv.gz > data/Vosa_cis_eQTLs.tsv
fi
