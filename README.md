# Marine cold-spells after iceberg B09B grounding (Adélie Sill–Commonwealth Bay)

This repository contains a single R script that reproduces the main analysis workflow: preprocessing satellite SST, detecting marine cold-spells using a standardized event framework, aggregating diagnostics, and generating the core time series and figures used in the manuscript.

## Contents

- `run_analysis.R` — workflow (data loading → diagnostics)

## What the script does

At a high level, the script:
1. Loads daily OSTIA foundation SST (and associated sea-ice fields where needed) for 1982–2024 over the Adélie Sill–Commonwealth Bay domain.
2. Detects marine cold-spell events following the event-based framework of Schlegel et al. (2021) (via `heatwaveR`).
3. Computes annual and seasonal summary metrics (with a focus on cumulative intensity).
4. Applies the MCS "ice" flag to identify near-ice cold-spell days.

## Data availability

The script is designed to run with publicly available datasets (e.g., OSTIA and supporting reanalysis/state-estimate products).  
To keep the repository lightweight, **input data are not stored here**. Instead, the script expects local paths to the required files.
