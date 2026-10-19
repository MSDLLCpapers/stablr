skip_on_cran()

test_that("predict.fit_kinetic correctly sanitizes intercept_cols and rate_cols", {
  # data("evers_supp")
  df <- evers_supp
  df <- df |> dplyr::filter(cqa=="Purity" & !is.na(y))
  df$batch <- paste("Batch",df$batch)
  df$pkg <- paste("Container-", df$pkg,sep="")
  fit <- fit_kinetic(df = df,
                     kinetic_order = 1,
                     trend = "decreasing",
                     Yinf = 0,
                     intercept_cols = "batch",
                     rate_cols = "pkg")
  
  testthat::expect_no_error(predict(fit))
  
  df_new <- data.frame(time=36,tempk=298.15,pkg="Container-PM2",batch="Batch Combi_F1-PM1")
  testthat::expect_no_error(predict(fit,df_new))
  
  df_new2 <- data.frame(time=36,tempk=298.15,pkg = "ContainerPM1",batch="BatchCombi_F2PM1")
  testthat::expect_error(predict(fit,df_new2))
})


test_that("predict.fit_singledeg correctly sanitizes intercept_cols and rate_cols", {
  # data("nrcesds")
  df <- nrcesds
  df_wide <- df |> tidyr::pivot_wider(id_cols=c("batch","time","temp","tempk"),names_from = "cqa",values_from = "y")
  df_wide$batch <- paste("Batch",df_wide$batch)
  fit <- fit_singledeg(df = df_wide,
                       trend = "decreasing",
                       Yinf = 0,
                       intercept_cols = "batch",
                       yA_var="NR Major Peak",
                       yB_var="NR Minor Peaks")
  
  testthat::expect_no_error(predict(fit))
  
  df_new <- data.frame(time=36,tempk=298.15,batch="Batch 2")
  testthat::expect_no_error(predict(fit,df_new))
  
  df_new2 <- data.frame(time=36,tempk=298.15,pkg = "Batch2")
  testthat::expect_error(predict(fit,df_new2))
})

test_that("predict.fit_parallel correctly sanitizes intercept_cols and rate_cols", {
  set.seed(100)
  # data("upsec")
  df <- upsec
  df_wide <- df |> tidyr::pivot_wider(id_cols=c("batch","time","temp","tempk"),names_from = "cqa",values_from = "y")
  df_wide$batch <- paste("Batch",df_wide$batch)
  df_wide$container <- sample(c("Container-1","Container-2"),size=nrow(df_wide),replace=TRUE)
  
  fit <- fit_parallel(df = df_wide,
                      intercept_cols = "batch",
                      rate_cols = "container",
                      yA_var = "Total Main",
                      yB_var = "HMWS",
                      yC_var = "LMWS",
                      trefk=313.15)
  
  testthat::expect_no_error(predict(fit))
  
  df_new <- data.frame(time=36,tempk=298.15,batch="Batch 2",container="Container-1")
  testthat::expect_no_error(predict(fit,df_new))
  
  df_new2 <- data.frame(time=36,tempk=298.15,batch = "Batch2", container="Container1")
  testthat::expect_error(predict(fit,df_new2))
})


test_that("predict.fit_kinetic correctly sanitizes numeric intercept_cols.", {
  # data("evers_supp")
  df <- evers_supp
  df <- df |> dplyr::filter(cqa=="Purity" & !is.na(y))
  df$batch <- ifelse(df$batch=="Combi_F1-PM1"| df$batch=="Combi_F1-PM2" | df$batch=="Combi_F2-PM1","000000-8302","00000000-8403")
  fit <- fit_kinetic(df = df,
                     kinetic_order = 1,
                     trend = "decreasing",
                     Yinf = 0,
                     intercept_cols = "batch")
  
  testthat::expect_no_error(predict(fit))
  
  df_new <- data.frame(time=36,tempk=298.15,batch="00000000-8403")
  testthat::expect_no_error(predict(fit,df_new))
})

test_that("predict.fit_sb correctly sanitizes intercept_cols and rate_cols", {
  # data("evers_supp")
  df <- evers_supp
  df <- df |> dplyr::filter(cqa=="Purity" & !is.na(y))
  df$batch <- paste("Batch",df$batch)
  df$pkg <- paste("Container-", df$pkg,sep="")
  fit <- fit_sb(df = df,
                     trend = "decreasing",
                     Yinf = 0,
                     intercept_cols = "batch",
                     rate_cols = "pkg")
  
  testthat::expect_no_error(predict(fit))
  
  df_new <- data.frame(time=36,tempk=298.15,pkg="Container-PM2",batch="Batch Combi_F1-PM1")
  testthat::expect_no_error(predict(fit,df_new))
})

test_that("levels differing only by punctuation are rejected, not pooled", {
  df <- evers_supp |> dplyr::filter(cqa=="Purity" & !is.na(y))

  # "Combi_F1-PM1" and "Combi_F1 PM1" both sanitize to "Combi_F1PM1", so they
  # would share one intercept parameter and be silently pooled.
  df_collide <- df
  df_collide$batch[df_collide$batch == "Combi_F2-PM1"] <- "Combi_F1 PM1"

  testthat::expect_error(
    fit_kinetic(df = df_collide, kinetic_order = 1, trend = "decreasing",
                Yinf = 0, intercept_cols = "batch"),
    "collide after removing punctuation"
  )

  # Distinct levels that merely contain punctuation must still fit.
  testthat::expect_no_error(
    fit_kinetic(df = df, kinetic_order = 1, trend = "decreasing",
                Yinf = 0, intercept_cols = "batch")
  )
})
