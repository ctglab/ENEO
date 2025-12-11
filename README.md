# ENEO

![GitHub release](https://img.shields.io/github/release/ctglab/ENEO.svg)
[![Snakemake](https://img.shields.io/badge/snakemake->8.0.0-brightgreen.svg)](https://snakemake.github.io)
[![Linting](https://github.com/ctglab/ENEO/actions/workflows/formatting.yml/badge.svg?branch=main)](https://github.com/ctglab/ENEO/actions/workflows/formatting.yml)
[![Testing](https://github.com/ctglab/ENEO/actions/workflows/testing.yml/badge.svg?branch=main)](https://github.com/ctglab/ENEO/actions/workflows/formatting.yml)
![Docker image building](https://github.com/ctglab/ENEO/actions/workflows/docker.yml/badge.svg?branch=main)
[![Conventional Commits](https://img.shields.io/badge/Conventional%20Commits-1.0.0-%23FE5196?logo=conventionalcommits&logoColor=white)](https://conventionalcommits.org)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14931774.svg)](https://doi.org/10.5281/zenodo.14931774)


ENEO (pronunced `/eˈnɛjo/`) is a Snakemake workflow developed for detecting immunogenic neoantigens arising from somatic mutations using only bulk tumor RNA-seq, without requiring matching samples or additional genomic data (WES/WGS). It uses a probabilistic model that leverages population genetics databases and genotype likelihoods to discriminate germline variants from the call set. Additional details are reported in the [publication](https://doi.org/10.1093/nargab/lqaf026).

## Usage

Pipeline documentation is available at [https://ctglab.github.io/ENEO](https://ctglab.github.io/ENEO)

## Cite

If you used this workflow in your work, don't forget to give credit to the authors by citing the original publication

>Danilo Tatoni, Mattia Dalsass, Giulia Brunelli, Mario Chiariello, Guido Grandi, Romina D’Aurizio, Efficient and effective identification of cancer neoantigens from tumor only RNA-seq, NAR Genomics and Bioinformatics, Volume 7, Issue 3, September 2025, lqaf026, https://doi.org/10.1093/nargab/lqaf026






