# Tests for simulate_kinetic_data()'s redesigned parameter-sampling layer.
#
# The generating parameters are now drawn from distributions matched to their
# physical support (Dirichlet/Beta for compositions, lognormal for rates), so
# closure, positivity, and bounds hold *by construction* — no clamp/rescale/floor
# on the parameters. These tests assert those invariants, target recovery, factor
# effects, the SD->kappa feasibility guard, and that every noise x clamp x
# reaction-type combination runs and returns the documented long format.

# ---- helpers ----------------------------------------------------------------

trajectory_wide <- function(dat) {
  tidyr::pivot_wider(dat, names_from = "cqa", values_from = "y")
}

# ---- closure at the parameter level -----------------------------------------

test_that("parallel intercept parameters sum to A_max within each group", {
  res <- simulate_kinetic_data(
    factors = c(batch = 8), intercept_cols = "batch",
    time_points = c(0, 6, 12), temperatures_C = c(5, 25, 40),
    trefk = 278.15, reaction_type = "parallel", kinetic_order = 1,
    trend = "decreasing",
    params = list(A0_target = 98, B0_target = 0.8, C0_target = 1.2,
                  A0_target_sd = 0.5, EaB = 85000, EaC = 80000,
                  krefB_median = 4e-4, krefB_gsd = 1.05,
                  krefC_median = 2e-4, krefC_gsd = 1.05, sigma_err = 0.3),
    sdpt = 1, return_params = TRUE
  )
  ip <- res$parameters$intercept_params
  expect_equal(ip$A0 + ip$B0 + ip$C0, rep(100, nrow(ip)), tolerance = 1e-10)
})

test_that("single_deg intercept parameters sum to A_max within each group", {
  res <- simulate_kinetic_data(
    factors = c(batch = 8), intercept_cols = "batch",
    time_points = c(0, 6, 12), temperatures_C = c(5, 25),
    trefk = 278.15, reaction_type = "single_deg", kinetic_order = 2,
    trend = "decreasing",
    params = list(A0_target = 98, B0_target = 2, A0_target_sd = 0.4,
                  Ea = 80000, kref_median = 5e-4, kref_gsd = 1.05,
                  sigma_err = 0.3),
    sdpt = 2, return_params = TRUE
  )
  ip <- res$parameters$intercept_params
  expect_equal(ip$A0 + ip$B0, rep(100, nrow(ip)), tolerance = 1e-10)
})

# ---- closure across the whole trajectory ------------------------------------

test_that("parallel species sum to A_max at every time x temp (mass_balance)", {
  dat <- simulate_kinetic_data(
    factors = c(batch = 4), intercept_cols = "batch",
    time_points = c(0, 3, 6, 12, 24), temperatures_C = c(5, 25, 40),
    trefk = 278.15, reaction_type = "parallel", kinetic_order = 1,
    trend = "decreasing",
    params = list(A0_target = 98, B0_target = 0.8, C0_target = 1.2,
                  A0_target_sd = 0.5, EaB = 85000, EaC = 80000,
                  krefB_median = 4e-4, krefB_gsd = 1.05,
                  krefC_median = 2e-4, krefC_gsd = 1.05, sigma_err = 0.4),
    noise = "mass_balance", clamp = FALSE, sdpt = 3
  )
  w <- trajectory_wide(dat)
  s <- w[["Total Main"]] + w[["HMWS"]] + w[["LMWS"]]
  expect_equal(s, rep(100, length(s)), tolerance = 1e-8)
})

test_that("single_deg species sum to A_max at every time x temp (mass_balance)", {
  dat <- simulate_kinetic_data(
    factors = c(batch = 4), intercept_cols = "batch",
    time_points = c(0, 3, 6, 12, 24), temperatures_C = c(5, 25, 40),
    trefk = 278.15, reaction_type = "single_deg", kinetic_order = 2,
    trend = "decreasing",
    params = list(A0_target = 98, B0_target = 2, A0_target_sd = 0.4,
                  Ea = 80000, kref_median = 5e-4, kref_gsd = 1.05,
                  sigma_err = 0.3),
    noise = "mass_balance", clamp = FALSE, sdpt = 4
  )
  w <- trajectory_wide(dat)
  s <- w[["NR Major Peak"]] + w[["NR Minor Peaks"]]
  expect_equal(s, rep(100, length(s)), tolerance = 1e-8)
})

# ---- positivity (no clamp / floor needed) -----------------------------------

test_that("drawn intercepts and rates are positive by construction", {
  res <- simulate_kinetic_data(
    factors = c(batch = 200), intercept_cols = "batch", rate_cols = "batch",
    time_points = c(0, 12), temperatures_C = 25,
    trefk = 278.15, reaction_type = "parallel", kinetic_order = 1,
    trend = "decreasing",
    # Small minor targets + large SD: the old normal+clamp/floor would have
    # produced boundary zeros here; the Dirichlet/lognormal cannot.
    params = list(A0_target = 99, B0_target = 0.4, C0_target = 0.6,
                  A0_target_sd = 1.0, EaB = 85000, EaC = 80000,
                  krefB_median = 5e-4, krefB_gsd = 2.0,
                  krefC_median = 2e-4, krefC_gsd = 2.0, sigma_err = 0.3),
    sdpt = 5, return_params = TRUE
  )
  ip <- res$parameters$intercept_params
  rp <- res$parameters$rate_params
  expect_true(all(ip$A0 > 0) && all(ip$B0 > 0) && all(ip$C0 > 0))
  expect_true(all(rp$kref_B > 0) && all(rp$kref_C > 0))
})

test_that("trajectory components stay positive without clamping", {
  dat <- simulate_kinetic_data(
    factors = c(batch = 6), intercept_cols = "batch",
    time_points = c(0, 3, 6, 12, 24, 36), temperatures_C = c(5, 25, 40),
    trefk = 278.15, reaction_type = "parallel", kinetic_order = 1,
    trend = "decreasing",
    params = list(A0_target = 98, B0_target = 0.8, C0_target = 1.2,
                  A0_target_sd = 0.5, EaB = 85000, EaC = 80000,
                  krefB_median = 4e-4, krefB_gsd = 1.05,
                  krefC_median = 2e-4, krefC_gsd = 1.05, sigma_err = 0.05),
    noise = "mass_balance", clamp = FALSE, sdpt = 6
  )
  expect_true(all(dat$y >= 0))
})

# ---- target-mean and target-SD recovery -------------------------------------

test_that("realised A0 mean and SD recover the requested targets at large N", {
  res <- simulate_kinetic_data(
    factors = c(batch = 4000), intercept_cols = "batch",
    time_points = 0, temperatures_C = 25,
    reaction_type = "kinetic", kinetic_order = 1, trend = "decreasing",
    params = list(A0_target = 90, A0_target_sd = 1.0, Ea = 80000,
                  kref_median = 2.5e-4, kref_gsd = 1.05, sigma_err = 0.3),
    sdpt = 7, return_params = TRUE
  )
  a0 <- res$parameters$intercept_params$A0
  expect_equal(mean(a0), 90, tolerance = 0.1)
  expect_equal(sd(a0), 1.0, tolerance = 0.1)
})

test_that("realised rate median recovers the requested lognormal median", {
  res <- simulate_kinetic_data(
    factors = c(batch = 4000), intercept_cols = "batch", rate_cols = "batch",
    time_points = 0, temperatures_C = 25,
    reaction_type = "kinetic", kinetic_order = 1, trend = "decreasing",
    params = list(A0_target = 90, A0_target_sd = 1.0, Ea = 80000,
                  kref_median = 2.5e-4, kref_gsd = 1.3, sigma_err = 0.3),
    sdpt = 8, return_params = TRUE
  )
  kref <- res$parameters$rate_params$kref
  expect_equal(median(kref), 2.5e-4, tolerance = 0.05)
})

# ---- factor effects on targets ----------------------------------------------

test_that("a factor drives the intercept target", {
  res <- simulate_kinetic_data(
    factors = c(site = 2, batch = 1500), intercept_cols = c("site", "batch"),
    time_points = 0, temperatures_C = 25,
    reaction_type = "kinetic", kinetic_order = 1, trend = "decreasing",
    params = list(A0_target = list(site = c(99, 95)), A0_target_sd = 0.3,
                  Ea = 80000, kref_median = 2.5e-4, kref_gsd = 1.05,
                  sigma_err = 0.3),
    sdpt = 9, return_params = TRUE
  )
  ip <- res$parameters$intercept_params
  m1 <- mean(ip$A0[ip$site == 1])
  m2 <- mean(ip$A0[ip$site == 2])
  expect_equal(m1, 99, tolerance = 0.1)
  expect_equal(m2, 95, tolerance = 0.1)
})

# ---- factor-dependent precision ---------------------------------------------

test_that("a factor can drive batch-to-batch precision", {
  res <- simulate_kinetic_data(
    factors = c(site = 2, batch = 2000), intercept_cols = c("site", "batch"),
    time_points = 0, temperatures_C = 25,
    reaction_type = "kinetic", kinetic_order = 1, trend = "decreasing",
    params = list(A0_target = 96, A0_target_sd = list(site = c(0.3, 1.2)),
                  Ea = 80000, kref_median = 2.5e-4, kref_gsd = 1.05,
                  sigma_err = 0.3),
    sdpt = 10, return_params = TRUE
  )
  ip <- res$parameters$intercept_params
  sd_tight <- sd(ip$A0[ip$site == 1])
  sd_loose <- sd(ip$A0[ip$site == 2])
  expect_equal(sd_tight, 0.3, tolerance = 0.1)
  expect_equal(sd_loose, 1.2, tolerance = 0.15)
  expect_gt(sd_loose, sd_tight)
})

# ---- SD -> kappa feasibility guard ------------------------------------------

test_that("an infeasibly large A0_target_sd is a hard error", {
  expect_error(
    simulate_kinetic_data(
      time_points = 0, temperatures_C = 25,
      reaction_type = "kinetic", kinetic_order = 1, trend = "decreasing",
      params = list(A0_target = 98, A0_target_sd = 50, Ea = 80000,
                    kref_median = 2.5e-4, kref_gsd = 1.05, sigma_err = 0.3)
    ),
    regexp = "too large"
  )
})

test_that("A0_target_sd = 0 pins A0 to the exact target (generic kinetic)", {
  res <- simulate_kinetic_data(
    factors = c(batch = 3), intercept_cols = "batch",
    time_points = 0, temperatures_C = 25,
    reaction_type = "kinetic", kinetic_order = 1, trend = "decreasing",
    params = list(A0_target = 98, A0_target_sd = 0, Ea = 80000,
                  kref_median = 2.5e-4, kref_gsd = 1.05, sigma_err = 0),
    sdpt = 11, return_params = TRUE
  )
  a0 <- res$parameters$intercept_params$A0
  expect_equal(a0, rep(98, length(a0)))
})

test_that("A0_target_sd = 0 pins the whole composition to its exact targets", {
  for (rt in c("single_deg", "parallel")) {
    params <- list(A0_target = 98, A0_target_sd = 0, sigma_err = 0)
    if (rt == "single_deg") {
      params <- c(params, list(B0_target = 2, Ea = 80000,
                               kref_median = 5e-4, kref_gsd = 1.05))
    } else {
      params <- c(params, list(B0_target = 0.8, C0_target = 1.2,
                               EaB = 85000, EaC = 80000,
                               krefB_median = 4e-4, krefB_gsd = 1.05,
                               krefC_median = 2e-4, krefC_gsd = 1.05))
    }
    res <- simulate_kinetic_data(
      factors = c(batch = 3), intercept_cols = "batch",
      time_points = c(0, 6, 12), temperatures_C = c(25, 40),
      reaction_type = rt, kinetic_order = 1,
      params = params, sdpt = 12, return_params = TRUE
    )
    ip <- res$parameters$intercept_params
    expect_equal(ip$A0, rep(98, nrow(ip)))
    expect_equal(ip$B0, rep(if (rt == "single_deg") 2 else 0.8, nrow(ip)))
    if (rt == "parallel") expect_equal(ip$C0, rep(1.2, nrow(ip)))
    expect_false(any(is.nan(res$data$y)))
  }
})

test_that("A0_target_sd mixed per level pins the zero-SD level and draws the rest", {
  res <- simulate_kinetic_data(
    factors = c(site = 2, batch = 2000), intercept_cols = c("site", "batch"),
    time_points = 0, temperatures_C = 25,
    reaction_type = "kinetic", kinetic_order = 1, trend = "decreasing",
    params = list(A0_target = 96, A0_target_sd = list(site = c(0, 0.8)),
                  Ea = 80000, kref_median = 2.5e-4, kref_gsd = 1.05,
                  sigma_err = 0.3),
    sdpt = 13, return_params = TRUE
  )
  ip <- res$parameters$intercept_params
  a0_pinned <- ip$A0[ip$site == 1]
  a0_drawn  <- ip$A0[ip$site == 2]
  expect_equal(a0_pinned, rep(96, length(a0_pinned)))
  expect_gt(sd(a0_drawn), 0)
})

test_that("a negative A0_target_sd is a hard error", {
  expect_error(
    simulate_kinetic_data(
      time_points = 0, temperatures_C = 25,
      reaction_type = "kinetic", kinetic_order = 1, trend = "decreasing",
      params = list(A0_target = 98, A0_target_sd = -1, Ea = 80000,
                    kref_median = 2.5e-4, kref_gsd = 1.05, sigma_err = 0.3)
    ),
    regexp = "must be >= 0"
  )
})

test_that("composition targets that do not sum to A_max are a hard error", {
  expect_error(
    simulate_kinetic_data(
      time_points = 0, temperatures_C = 25,
      reaction_type = "single_deg", kinetic_order = 1, trend = "decreasing",
      params = list(A0_target = 98, B0_target = 5, A0_target_sd = 0.4,
                    Ea = 80000, kref_median = 5e-4, kref_gsd = 1.05,
                    sigma_err = 0.3)
    ),
    regexp = "must sum to A_max"
  )
})

# ---- factor level counts must be positive integers --------------------------

test_that("a non-positive or non-integer factor level count is a hard error", {
  base_params <- list(A0_target = 98, A0_target_sd = 0.5, Ea = 80000,
                      kref_median = 2.5e-4, kref_gsd = 1.05, sigma_err = 0.1)
  for (bad in c(0, -2, 2.7, NA_real_)) {
    expect_error(
      simulate_kinetic_data(
        factors = c(batch = bad), intercept_cols = "batch",
        time_points = c(0, 6), temperatures_C = 25,
        reaction_type = "kinetic", kinetic_order = 1, trend = "decreasing",
        params = base_params
      ),
      regexp = "positive integer"
    )
  }
})

test_that("valid positive-integer level counts (including 1) are accepted", {
  base_params <- list(A0_target = 98, A0_target_sd = 0.5, Ea = 80000,
                      kref_median = 2.5e-4, kref_gsd = 1.05, sigma_err = 0.1)
  for (n in c(1, 3)) {
    res <- simulate_kinetic_data(
      factors = c(batch = n), intercept_cols = "batch",
      time_points = c(0, 6), temperatures_C = 25,
      reaction_type = "kinetic", kinetic_order = 1, trend = "decreasing",
      params = base_params, sdpt = 1
    )
    expect_setequal(levels(res$batch), as.character(seq_len(n)))
  }
})

# ---- sigma_err must be a scalar ---------------------------------------------

test_that("a vector sigma_err is a hard error (no silent recycling)", {
  expect_error(
    simulate_kinetic_data(
      time_points = c(0, 1, 2), temperatures_C = c(25, 40),
      reaction_type = "kinetic", kinetic_order = 1, trend = "decreasing",
      params = list(A0_target = 50, A0_target_sd = 0.1, Ea = 80000,
                    kref_median = 2.5e-4, kref_gsd = 1.05,
                    sigma_err = c(1e-4, 10))
    ),
    regexp = "single non-negative numeric"
  )
})

test_that("a named-list sigma_err is a hard error with a clear message", {
  expect_error(
    simulate_kinetic_data(
      factors = c(batch = 2), intercept_cols = "batch",
      time_points = c(0, 1, 2), temperatures_C = c(25, 40),
      reaction_type = "kinetic", kinetic_order = 1, trend = "decreasing",
      params = list(A0_target = 50, A0_target_sd = 0.1, Ea = 80000,
                    kref_median = 2.5e-4, kref_gsd = 1.05,
                    sigma_err = list(batch = c(0.01, 5)))
    ),
    regexp = "single non-negative numeric"
  )
})

test_that("an NA or negative sigma_err is a hard error", {
  base_params <- list(A0_target = 50, A0_target_sd = 0.1, Ea = 80000,
                      kref_median = 2.5e-4, kref_gsd = 1.05)
  for (bad in list(NA_real_, -1)) {
    expect_error(
      simulate_kinetic_data(
        time_points = 0, temperatures_C = 25,
        reaction_type = "kinetic", kinetic_order = 1, trend = "decreasing",
        params = c(base_params, list(sigma_err = bad))
      ),
      regexp = "single non-negative numeric"
    )
  }
})

test_that("a scalar sigma_err (including 0) is accepted for every reaction type", {
  specs <- list(
    kinetic    = list(A0_target = 98, A0_target_sd = 0.5, Ea = 80000,
                      kref_median = 2.5e-4, kref_gsd = 1.05),
    single_deg = list(A0_target = 98, B0_target = 2, A0_target_sd = 0.4,
                      Ea = 80000, kref_median = 5e-4, kref_gsd = 1.05),
    parallel   = list(A0_target = 98, B0_target = 0.8, C0_target = 1.2,
                      A0_target_sd = 0.5, EaB = 85000, EaC = 80000,
                      krefB_median = 4e-4, krefB_gsd = 1.05,
                      krefC_median = 2e-4, krefC_gsd = 1.05)
  )
  for (rt in names(specs)) {
    for (se in c(0, 0.3)) {
      res <- simulate_kinetic_data(
        time_points = c(0, 6, 12), temperatures_C = c(25, 40),
        reaction_type = rt, kinetic_order = 1,
        params = c(specs[[rt]], list(sigma_err = se)), sdpt = 5
      )
      expect_s3_class(res, "data.frame")
      expect_false(any(is.na(res$y)))
    }
  }
})

# ---- noise modes x clamp x reaction type ------------------------------------

test_that("every noise x clamp x reaction-type combination runs and is well-formed", {
  base <- list(
    kinetic = list(
      reaction_type = "kinetic",
      params = list(A0_target = 98, A0_target_sd = 0.5, Ea = 80000,
                    kref_median = 2.5e-4, kref_gsd = 1.05, sigma_err = 0.4),
      cqas = "Purity"
    ),
    single_deg = list(
      reaction_type = "single_deg",
      params = list(A0_target = 98, B0_target = 2, A0_target_sd = 0.4,
                    Ea = 80000, kref_median = 5e-4, kref_gsd = 1.05,
                    sigma_err = 0.4),
      cqas = c("NR Major Peak", "NR Minor Peaks")
    ),
    parallel = list(
      reaction_type = "parallel",
      params = list(A0_target = 98, B0_target = 0.8, C0_target = 1.2,
                    A0_target_sd = 0.5, EaB = 85000, EaC = 80000,
                    krefB_median = 4e-4, krefB_gsd = 1.05,
                    krefC_median = 2e-4, krefC_gsd = 1.05, sigma_err = 0.4),
      cqas = c("Total Main", "HMWS", "LMWS")
    )
  )

  seed <- 100
  for (rt in names(base)) {
    for (nz in c("mass_balance", "independent")) {
      for (cl in c(FALSE, TRUE)) {
        seed <- seed + 1
        spec <- base[[rt]]
        dat <- simulate_kinetic_data(
          factors = c(batch = 4), intercept_cols = "batch",
          time_points = c(0, 6, 12), temperatures_C = c(5, 25, 40),
          trefk = 278.15, reaction_type = spec$reaction_type,
          kinetic_order = 1, trend = "decreasing", params = spec$params,
          noise = nz, clamp = cl, sdpt = seed
        )
        label <- sprintf("%s / %s / clamp=%s", rt, nz, cl)
        expect_setequal(names(dat), c("batch", "time", "temp", "tempk", "cqa", "y"))
        expect_setequal(unique(dat$cqa), spec$cqas)
        expect_true(all(is.finite(dat$y)), info = label)
        if (cl) {
          expect_true(all(dat$y >= 0 & dat$y <= 100), info = label)
        }
      }
    }
  }
})

# ---- parallel has no valid increasing trend --------------------------------

test_that("explicit trend = \"increasing\" errors for reaction_type = \"parallel\"", {
  params <- list(A0_target = 98, B0_target = 0.8, C0_target = 1.2,
                 A0_target_sd = 0.5, EaB = 85000, EaC = 80000,
                 krefB_median = 4e-4, krefB_gsd = 1.05,
                 krefC_median = 2e-4, krefC_gsd = 1.05, sigma_err = 0.4)
  expect_error(
    simulate_kinetic_data(
      time_points = c(0, 6, 12), temperatures_C = c(25, 40),
      reaction_type = "parallel", kinetic_order = 1,
      trend = "increasing", params = params, sdpt = 1
    ),
    "not supported for reaction_type = \"parallel\""
  )
})

test_that("a bare parallel call defaults to the safe decreasing path", {
  params <- list(A0_target = 98, B0_target = 0.8, C0_target = 1.2,
                 A0_target_sd = 0, EaB = 85000, EaC = 80000,
                 krefB_median = 4e-4, krefB_gsd = 1.05,
                 krefC_median = 2e-4, krefC_gsd = 1.05, sigma_err = 0)
  # No `trend` supplied: the parallel default is forced to decreasing, so the
  # source degrades and every species stays within [0, A_max].
  dat <- simulate_kinetic_data(
    factors = c(batch = 3), intercept_cols = "batch",
    time_points = c(0, 6, 12, 24), temperatures_C = c(25, 40),
    reaction_type = "parallel", kinetic_order = 1,
    params = params, sdpt = 2
  )
  w <- trajectory_wide(dat)
  expect_true(all(w[["Total Main"]] <= 98 + 1e-8))  # source only decays
  expect_true(all(dat$y >= -1e-8))                  # no negative concentrations
})
