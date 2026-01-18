# ARIMA Time Series Analysis: French Arms Production Index

This repositoty implements a time series analysis and forecasting of French arms and ammunition production index using ARIMA models, demonstrating stationarity testing, model selection via information criteria, and residual diagnostics.

## What This Demonstrates

- **Time series preprocessing**: Stationarity testing using ADF and PP tests, first-order differencing to achieve stationarity
- **ARIMA model selection**: Systematic validation of ARMA(p,q) models with significance testing and autocorrelation diagnostics
- **Information criteria optimization**: Model selection using AIC and BIC to identify optimal ARIMA specifications
- **Residual diagnostics**: Ljung-Box tests, Shapiro-Wilk normality tests, and Q-Q plots for model validation
- **Forecasting**: Multi-step ahead predictions with confidence intervals using the selected ARIMA model
- **R programming**: Custom functions for model validation, automated model selection, and statistical testing

## Key Results

- **Best models identified**: ARIMA(4,1,1) minimizes AIC, ARIMA(0,1,1) minimizes BIC for the differenced series
- **Stationarity achieved**: First-order differencing successfully transformed the non-stationary series
- **Model diagnostics**: Residuals pass autocorrelation tests (Ljung-Box) and show acceptable fit
- **Visualizations**: See `reports/` for ACF/PACF plots, residual plots, Q-Q plots, and forecast visualizations
- **Full analysis report**: See `reports/Time_series__axel_jeje_.pdf` for detailed methodology and results

## Repository Layout

```
arima_lts_project/
├── src/                # Analysis scripts
│   └── project_LTS.R   # Main ARIMA analysis script
├── data/               # Time series datasets (see data/README.md)
│   ├── arms_truncated_index.csv
│   └── README.md
├── reports/            # Generated figures and PDF reports
│   ├── Time_series__axel_jeje_.pdf
│   ├── plot_*.png      # Time series and diagnostic plots
│   └── *.png           # Model comparison figures
├── tests/              # (Reserved for future test suites)
├── requirements.txt    # R package dependencies
├── LICENSE             # MIT License
└── README.md           # This file
```

## Quickstart

### Prerequisites

- R (>= 3.6.0)
- Required R packages (see `requirements.txt`)

### Installation

```r
# Install required packages
install.packages(c("zoo", "tseries", "fUnitRoots", "ggplot2", "xtable", "forecast"))
```

### Running the Analysis

```r
# From the project root directory
setwd("src")
source("project_LTS.R")
```

**Note**: Ensure data files are placed in the `data/` directory before running (see `data/README.md`).

### Demo Command

```bash
# Quick demo: Run the main analysis script
cd src && Rscript project_LTS.R
```

## Method Overview

1. **Data Loading**: French arms production index (2008-2024), monthly frequency
2. **Exploratory Analysis**: Visual inspection, ACF/PACF analysis, trend detection
3. **Stationarity Testing**: ADF and PP tests confirm non-stationarity; first-order differencing applied
4. **Model Selection**: 
   - Systematic ARMA(p,q) grid search (p≤4, q≤5)
   - Validation criteria: significant coefficients, no residual autocorrelation
   - Information criteria (AIC/BIC) used for final model selection
5. **Model Estimation**: ARIMA(4,1,1) and ARIMA(0,1,1) estimated on original series
6. **Diagnostics**: Residual autocorrelation tests, normality tests, Q-Q plots
7. **Forecasting**: 6-step ahead predictions with 95% confidence intervals

## Notes / Limitations

- **Data availability**: Original CSV files not included in repository (see `data/README.md` for instructions)
- **Model scope**: Analysis limited to ARIMA models; no comparison with other time series methods (e.g., GARCH, state space models)
- **Sample size**: 193 monthly observations (2008-2024) may limit model complexity
- **External validation**: No out-of-sample validation or walk-forward analysis performed
- **Regime changes**: Analysis assumes stationarity after differencing; potential structural breaks not explicitly tested

## References

- Box, G. E. P., & Jenkins, G. M. (1976). *Time Series Analysis: Forecasting and Control*. Holden-Day.
- Hamilton, J. D. (1994). *Time Series Analysis*. Princeton University Press.
- Hyndman, R. J., & Athanasopoulos, G. (2021). *Forecasting: principles and practice* (3rd ed.). OTexts.
