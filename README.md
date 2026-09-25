# stablr <img src="man/figures/stablr_logo.png" align="right" height="139" alt="" />


## Overview

**stablr** is an R package for kinetic modeling of stability data. It provides a comprehensive framework for fitting and analyzing various degradation reaction models, incorporating the modified Arrhenius equation for temperature-dependent kinetics.

The package supports multiple kinetic reaction models:

- **Kinetic Models (single species)**: Zero-order, first-order, and higher-order kinetics for any one measured species (a decreasing parent or an increasing product)
- **Single Degradation Reactions (A → B)**: With both species A and B measured
- **Šesták–Berggren Model (single species)**: Flexible kinetic modeling with autocatalytic parameter
- **Parallel Reaction Models (A → B, C)**: For products with multiple degradation pathways

## Key Features

- 🔬 **Multiple kinetic orders**: Zero, first, and higher-order kinetics
- 📊 **Flexible model specification**: Estimate rate and intercept by experimental factors (formulation, packaging, batch, etc.)
- ⚡ **Autocatalytic reactions**: Package supports Šesták–Berggren Model with simple degradation and autocatalytic processes
- 🌡️ **Temperature and humidity modeling**:  Arrhenius equation for temperature- and humidity-dependent rate constants
- 📦 **Real data examples**: Includes published and simulated stability datasets
- 📈 **Time out of storage modeling**: Model temperature and/or humidity deviations during product storage

## Installation

You can install the development version of stablr from GitHub:

```r
# install.packages("devtools")
devtools::install_github("MSDLLCPapers/stablr")
```

## Usage

### Basic Example: Generic Kinetic Model

Fit a kinetic model with varying intercepts by batch and rates by packaging:

```r
library(stablr)
library(dplyr)

# Filter and prepare data
df <- evers_supp %>% 
  filter(cqa == "Purity", !is.na(y)) %>% 
  mutate(
    batch = as.numeric(as.factor(batch)),
    pkg = as.numeric(as.factor(pkg))
  ) 

# Fit zero-order kinetic model
fit <- fit_kinetic(
  df,
  kinetic_order = 0, 
  trend = "decreasing",
  intercept_cols = "batch", 
  rate_cols = "pkg"
)

summary(fit)
```

### Single Degradation Model (A → B)

When both the parent compound (A) and degradation product (B) are measured:

```r
library(tidyr)

# Prepare data in wide format
nrcesds_wide <- nrcesds %>% 
  pivot_wider(names_from = "cqa", values_from = "y")

# Fit single degradation model
fit <- fit_singledeg(
  df = nrcesds_wide,
  kinetic_order = 2,
  trend = "decreasing",
  yA_var = "NR Major Peak",
  yB_var = "NR Minor Peaks"
)

summary(fit)
```

### Šesták–Berggren (SB) Model (single species)

The SB model provides maximum flexibility for complex degradation kinetics by estimating both the reaction order (*n*) and an autocatalytic parameter (*m*). It is useful when standard kinetic orders (0, 1, 2) don't fit well, for autocatalytic reactions, and for exploratory kinetic analysis. See `vignette("stablr-introduction", package = "stablr")` for the model equation and worked examples.

```r
# Fit Šesták–Berggren model
fit_sb <- fit_sb(
  df,
  trend = "decreasing",
  intercept_cols = "batch",
  rate_cols = "pkg"
)

summary(fit_sb)
```

### Parallel Reaction Model (A → B, C)

For products with multiple degradation pathways:

```r
# Prepare data
upsec_wide <- upsec %>% 
  pivot_wider(names_from = "cqa", values_from = "y")

# Fit parallel reaction model
fit <- fit_parallel(
  df = upsec_wide,
  kinetic_order = 1,
  trend = "decreasing",
  yA_var = "Total Main",
  yB_var = "HMWS",
  yC_var = "LMWS"
)

summary(fit)
```

## Prediction and Intervals

Every fitted model works with the standard `predict()` interface, which supports
confidence, prediction, and tolerance intervals:

```r
# New conditions to predict at
pred_grid <- data.frame(time = seq(0, 36, by = 3), tempk = 278.15)

# 95% prediction interval for a future observation
predict(fit, newdata = pred_grid, interval = "prediction", level = 0.95)

# Tolerance interval: cover 90% of the population with 95% confidence
predict(fit, newdata = pred_grid, interval = "tolerance",
        coverage = 0.90, level = 0.95, side = "two")
```

See `vignette("intervals-calculation", package = "stablr")` for details.

## Time Out of Storage (TOS) Predictions

Use `predict_tos()` to project an attribute through a multi-segment temperature
schedule (e.g. a cold-chain excursion), where each segment has a start time, end
time, and temperature:

```r
# Temperature schedule: 12 months at 5C, 1 month excursion at 40C, back to 5C
schedule <- data.frame(
  t_start = c(0, 12, 13),
  t_end   = c(12, 13, 24),
  temp    = c(5, 40, 5)
)

predict_tos(fit, schedule = schedule, temp_unit = "C")
```

See `vignette("predict-tos", package = "stablr")` for details.

## Simulating Stability Data

Generate synthetic stability data for model validation and testing:

```r
# Simulate data with formulation and packaging factors
params  <- list(
  A0_target    = 98, 
  A0_target_sd = 0.4,
  Ea           = 80000,
  kref_median  = 5e-4,
  kref_gsd     = 1.05,
  sigma_err    = 0.3
)

sim_data <- simulate_kinetic_data(
  factors = c("conc" = 3, "pkg" = 2),
  intercept_cols = "conc",
  rate_cols = c("conc", "pkg"),
  trend = "decreasing",
  time_points = c(0:12),
  temperatures_C = c(5,25,40),
  reaction_type = "kinetic",
  kinetic_order = 1,
  params=params
)
```

See `vignette("simulating-stability-data", package = "stablr")` for a full walkthrough.

## Included Datasets

The package includes several pharmaceutical stability datasets:

- **`evers_supp`** - HPLC purity and SEC HMW stability data for therapeutic peptide SAR441255 from Evers et al. (2022), *Pharmaceutics*
- **`nrcesds`** - Simulated non-reduced CE-SDS data generated from a single-degradation second-order reaction model
- **`upsec`** - Simulated UPSEC (ultraperformance size-exclusion chromatography) data generated from a parallel first-order reactions model
- **`flavier_supp`** - Colour-change (ΔE\*) data across temperature and relative humidity, used to demonstrate the humidity-modified Arrhenius model

## Model Building Functions

For advanced users who need more control, **stablr** provides lower-level model building functions:

- `build_kinetic()` - Build generic kinetic model components
- `build_singledeg()` - Build single degradation model components
- `build_parallel()` - Build parallel reaction model components
- `expr_intercept()` - Generate intercept expressions with grouping
- `expr_rate()` - Generate rate constant expressions with grouping
- `expr_arrhenius()` - Generate Arrhenius temperature relationship expressions
- `expr_arrhenius_humidity()` - Generate humidity-modified Arrhenius expressions

## Documentation

For detailed examples and theoretical background, see the package vignettes:

```r
# View available vignettes
browseVignettes("stablr")

vignette("stablr-introduction", package = "stablr")         # Start here
vignette("kinetic-model-comparison", package = "stablr")    # Compare kinetic orders
vignette("simulating-stability-data", package = "stablr")   # Simulate datasets
vignette("intervals-calculation", package = "stablr")       # Prediction & tolerance intervals
vignette("predict-tos", package = "stablr")                 # Time-out-of-storage prediction
vignette("humidity-modified-arrhenius", package = "stablr") # Humidity-modified Arrhenius
```


