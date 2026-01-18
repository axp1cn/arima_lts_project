# Data Directory

This directory contains the time series data used in the ARIMA modeling project.

## Data Files

- `arms_truncated_index.csv`: French arms and ammunition production index data (2008-2024)
- `arms_index.csv`: Full arms production index dataset
- `jus_index.csv`: Additional index data

## Data Source

The data contains French industrial production indices, specifically focusing on arms and ammunition production from 2008 to 2024.

## Usage

The main analysis uses `arms_truncated_index.csv`. The data should be placed in this directory before running the analysis scripts.

## Note

For privacy and licensing reasons, the actual CSV data files are not included in this repository. To run the analysis:

1. Obtain the data from the original source
2. Place the CSV files in this `data/` directory
3. Ensure the file format matches the expected structure (semicolon-separated with date and index columns)
