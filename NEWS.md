# stablr 1.0.1

First public release of **stablr**, a framework for kinetic modeling of
pharmaceutical stability data with Arrhenius temperature dependence and
statistical interval estimation.

## Models

* `fit_kinetic()` — generic single-species kinetics (a decreasing parent or an
  increasing product) for zero-, first-, and higher-order reactions.
* `fit_singledeg()` — single degradation reaction (A → B) where both the parent
  and the degradant are measured.
* `fit_parallel()` — parallel reactions (A → B, C) for a source that degrades
  down two pathways simultaneously.
* `fit_sb()` — Šesták–Berggren α-space ODE model, estimating both the reaction
  order (*n*) and an autocatalytic parameter (*m*) for complex or autocatalytic
  kinetics.

## Temperature and covariate modeling

* Modified Arrhenius temperature dependence of rate constants via `expr_arrhenius()`,
  with a humidity-extended form in `expr_arrhenius_humidity()`. All fitters accept
  a pluggable `temp_expr` argument.
* Rates and intercepts can vary by experimental factors (formulation, packaging,
  batch, etc.) through `intercept_cols` and `rate_cols`.

## Prediction and intervals

* `predict()` methods for every model class, returning fitted values with
  confidence, prediction, and tolerance intervals. Tolerance intervals follow
  Krishnamoorthy & Mathew (2009), with two-sided (Wallis 1951) and one-sided
  variants.
* `predict_tos()` projects an attribute through a multi-segment time-out-of-storage
  temperature schedule (e.g. cold-chain excursions).

## Simulation and data

* `simulate_kinetic_data()` generates synthetic stability data for all model
  types, for validation and study design.
* Lower-level model builders — `build_kinetic()`, `build_singledeg()`,
  `build_parallel()`, `build_alpha_ode()` — and expression builders
  (`expr_intercept()`, `expr_rate()`) are exported for advanced use.
* Bundled datasets: `evers_supp`, `nrcesds`, `upsec`, and `flavier_supp`.

## Documentation

* Six vignettes covering an introduction, kinetic-order comparison, data
  simulation, interval computation, time-out-of-storage prediction, and the
  humidity-modified Arrhenius model.
