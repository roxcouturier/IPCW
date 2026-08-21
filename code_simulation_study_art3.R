# ============================================================
# Simulation study: Cox models under informative censoring
# Comparing unweighted, IPCW (full / partial), and standardized
# covariate-adjusted Cox estimators of a treatment log hazard ratio.
#
# Required packages: survival, ipw, parallel, dplyr, smd, cobalt
# ============================================================

library(survival)
library(ipw)
library(parallel)
library(dplyr)
library(smd)
library(cobalt)

# ============================================================
# 1. Data generation
# ------------------------------------------------------------
# Two independent exponential event times are simulated per subject:
#   Tx : time to the event of interest
#   Ta : time to a competing censoring process
# Both hazards depend on treatment (trt) and a binary covariate (Z1),
# which makes censoring potentially informative.
# ============================================================
generate_data_simple <- function(n = 1000, hx0 = 0.115, beta_x = log(0.75),
                                 ha0 = 0.06192308, beta_a = log(1),
                                 betax_z1 = log(0.5), betaa_z1 = log(1),
                                 seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  trt <- rbinom(n, 1, 0.5)
  Z1  <- rbinom(n, 1, 0.5)
  
  # Hazard for the event of interest and for the competing/censoring process
  rate_x <- hx0 * exp(beta_x * trt + betax_z1 * Z1)
  rate_a <- ha0 * exp(beta_a * trt + betaa_z1 * Z1)
  
  Tx <- rexp(n, rate = rate_x)
  Ta <- rexp(n, rate = rate_a)
  
  time   <- pmin(Tx, Ta)
  status <- ifelse(Tx <= Ta, 1, 0)  # 1 = event of interest, 0 = censored
  
  data.frame(trt = trt, Z1 = Z1, time = time, status = status)
}

# ============================================================
# 2. One simulation replication
# ------------------------------------------------------------
# Generates one dataset and fits five Cox models:
#   - unweighted   : no adjustment for censoring
#   - ipw_both     : IPCW, censoring model on trt + Z1
#   - ipw_trtonly  : IPCW, censoring model on trt only
#   - ipw_Z1only   : IPCW, censoring model on Z1 only
#   - cox_adjusted : covariate-adjusted Cox + standardization (marginal HR)
# ============================================================
simulate_cox <- function(n = 500, beta_x, beta_a, betax_z1, betaa_z1) {
  
  data <- generate_data_simple(
    n = n, hx0 = 0.115, beta_x = beta_x, ha0 = 0.06192308,
    beta_a = beta_a, betax_z1 = betax_z1, betaa_z1 = betaa_z1
  )
  
  data$censored <- 1 - data$status
  data$Tstart   <- 0
  
  # ---- Reshape to long (counting-process) format ----
  # One long dataset carries the event indicator, the other the
  # censoring indicator, on the same time grid (required by ipwtm()).
  times <- sort(unique(data$time))
  
  data.long <- survSplit(data, cut = times, end = "time", start = "Tstart",
                         event = "status", id = "id")
  
  data.long.c <- survSplit(data, cut = times, end = "time", start = "Tstart",
                           event = "censored", id = "id")
  data.long$censored <- data.long.c$censored
  
  # ---- IPCW weights: three censoring models ----
  w <- ipwtm(
    exposure = censored, family = "survival",
    numerator = ~ 1, denominator = ~ trt + Z1,
    id = id, tstart = Tstart, timevar = time,
    type = "cens", data = data.long
  )$ipw.weights
  
  w_onlytrt <- ipwtm(
    exposure = censored, family = "survival",
    numerator = ~ 1, denominator = ~ trt,
    id = id, tstart = Tstart, timevar = time,
    type = "cens", data = data.long
  )$ipw.weights
  
  w_onlyZ1 <- ipwtm(
    exposure = censored, family = "survival",
    numerator = ~ 1, denominator = ~ Z1,
    id = id, tstart = Tstart, timevar = time,
    type = "cens", data = data.long
  )$ipw.weights
  
  data.long$ipw        <- w
  data.long$ipw_onlytrt <- w_onlytrt
  data.long$ipw_onlyZ1  <- w_onlyZ1
  
  # ---- Cox models ----
  cox_unweighted <- coxph(Surv(Tstart, time, status) ~ trt + cluster(id),
                          data = data.long)
  
  cox_ipw <- coxph(Surv(Tstart, time, status) ~ trt + cluster(id),
                   data = data.long, weights = ipw)
  
  cox_ipw_trtonly <- coxph(Surv(Tstart, time, status) ~ trt + cluster(id),
                           data = data.long, weights = ipw_onlytrt)
  
  cox_ipw_Z1only <- coxph(Surv(Tstart, time, status) ~ trt + cluster(id),
                          data = data.long, weights = ipw_onlyZ1)
  
  cox_adjusted <- coxph(Surv(time, status) ~ trt + Z1, data = data)
  
  # ---- Standardization (g-computation) from the adjusted Cox model ----
  # Predict survival under two counterfactual scenarios (everyone
  # treated vs. nobody treated), average the survival curves, and
  # convert to a marginal (standardized) hazard ratio.
  data_trt1 <- data; data_trt1$trt <- 1
  data_trt0 <- data; data_trt0$trt <- 0
  t_eval <- median(data$time)
  
  surv1 <- summary(survfit(cox_adjusted, newdata = data_trt1),
                   times = t_eval, extend = TRUE)$surv
  surv0 <- summary(survfit(cox_adjusted, newdata = data_trt0),
                   times = t_eval, extend = TRUE)$surv
  
  log_HR_std <- log(-log(mean(surv1)) / -log(mean(surv0)))
  
  # ---- Return the treatment effect estimates from all five models ----
  c(
    unweighted   = coef(cox_unweighted),
    ipw_both     = coef(cox_ipw),
    ipw_trtonly  = coef(cox_ipw_trtonly),
    ipw_Z1only   = coef(cox_ipw_Z1only),
    cox_adjusted = log_HR_std
  )
}

# ============================================================
# 3. Replication wrapper with retry on failure
# ------------------------------------------------------------
# Model fitting can occasionally fail (e.g. convergence issues in the
# censoring model). This wrapper retries with a fresh draw until it
# succeeds or the retry limit is reached.
# ============================================================
single_simulation_with_retry <- function(sim_id, n, beta_x, beta_a,
                                         betax_z1, betaa_z1,
                                         max_retries = 100) {
  attempt <- 0
  while (attempt < max_retries) {
    out <- tryCatch({
      simulate_cox(n = n, beta_x = beta_x, beta_a = beta_a,
                   betax_z1 = betax_z1, betaa_z1 = betaa_z1)
    }, error = function(e) NULL)
    
    if (!is.null(out)) return(out)
    attempt <- attempt + 1
  }
  warning(paste("Simulation", sim_id, "failed after", max_retries, "attempts"))
  return(NULL)
}

# ============================================================
# 4. Full simulation across a grid of scenarios (parallelized)
# ------------------------------------------------------------
# Scenario grid dimensions:
#   HR_trt_event : effect of treatment on the event of interest
#   HR_trt_allo  : effect of treatment on the competing/censoring process
#   HR_x_z1      : effect of Z1 on the event of interest
#   HR_a_z1      : effect of Z1 on the competing/censoring process
# ============================================================
run_simulation_full_parallel <- function(
    HR_trt_event_vals = c(0.75),
    HR_trt_allo_vals  = c(1, 2, 3),
    HR_x_z1_vals      = c(0.5, 1, 2),
    HR_a_z1_vals      = c(0.2, 0.6, 1, 1.4, 2.2, 3),
    nsim = 10,
    n = 500,
    n_cores = 6
) {
  message("Using ", n_cores, " cores")
  
  results_all <- list()
  scenario_id <- 1
  
  for (hr_trt_event in HR_trt_event_vals) {
    for (hr_trt_allo in HR_trt_allo_vals) {
      for (HR_x in HR_x_z1_vals) {
        for (HR_a in HR_a_z1_vals) {
          
          message("\n===== Scenario ", scenario_id,
                  " : HR_trt_event=", hr_trt_event,
                  ", HR_trt_allo=", hr_trt_allo,
                  ", HR_x_z1=", HR_x,
                  ", HR_a_z1=", HR_a, " =====")
          
          beta_x   <- log(hr_trt_event)
          beta_a   <- log(hr_trt_allo)
          betax_z1 <- log(HR_x)
          betaa_z1 <- log(HR_a)
          
          if (.Platform$OS.type == "unix") {
            out_list <- parallel::mclapply(
              1:nsim,
              function(i) single_simulation_with_retry(i, n, beta_x, beta_a,
                                                       betax_z1, betaa_z1),
              mc.cores = n_cores, mc.set.seed = TRUE
            )
          } else {
            cl <- parallel::makeCluster(n_cores)
            parallel::clusterEvalQ(cl, {
              library(survival); library(ipw); library(dplyr)
              library(smd); library(cobalt)
            })
            parallel::clusterExport(
              cl,
              varlist = c("single_simulation_with_retry", "simulate_cox",
                          "generate_data_simple", "n", "beta_x",
                          "beta_a", "betax_z1", "betaa_z1"),
              envir = environment()
            )
            out_list <- parallel::parLapply(
              cl, 1:nsim,
              function(i) single_simulation_with_retry(i, n, beta_x, beta_a,
                                                       betax_z1, betaa_z1)
            )
            parallel::stopCluster(cl)
          }
          
          n_failures <- sum(sapply(out_list, is.null))
          if (n_failures > 0) {
            warning(paste("Scenario", scenario_id, ":", n_failures,
                          "simulation(s) failed"))
          }
          
          out_list <- Filter(Negate(is.null), out_list)
          res <- as.data.frame(do.call(rbind, out_list))
          res$HR_trt_event <- hr_trt_event
          res$HR_trt_allo  <- hr_trt_allo
          res$HR_x_z1      <- HR_x
          res$HR_a_z1      <- HR_a
          
          message("Successful simulations: ", nrow(res), " / ", nsim)
          
          results_all[[scenario_id]] <- res
          scenario_id <- scenario_id + 1
        }
      }
    }
  }
  
  do.call(rbind, results_all)
}

# ============================================================
# 5. Execution
# ------------------------------------------------------------
# ============================================================
set.seed(2025)
results <- run_simulation_full_parallel(
  HR_trt_event_vals = c(0.75), #changer 1 or 0.5 
  HR_trt_allo_vals  = c(1, 2, 3),
  HR_x_z1_vals      = c(0.5, 1, 2),
  HR_a_z1_vals      = c(0.2, 0.6, 1, 1.4, 2.2, 3),
  nsim    = 10000,
  n       = 500,
  n_cores = 6
)

save(results, file = "results_ipcw_0.75.RData")



# ============================================================
# plot_bias_figures.R
#
# Produces three bias figures from the simulation results
# (columns: unweighted.trt, ipw_both.trt, ipw_trtonly.trt,
#  ipw_Z1only.trt, cox_adjusted, HR_trt_event, HR_trt_allo,
#  HR_x_z1, HR_a_z1), as returned by simulate_cox().
#
#   Figure 1 : Cox (unweighted)  vs  IPCW (full)
#   Figure 2 : Cox (unweighted)  vs  IPCW (full, trt-only, Z1-only)
#   Figure 3 : Cox (unweighted)  vs  IPCW (full)  vs  Standardized adjusted Cox
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggh4x)

# ---- 1. Load simulation results ----
load("results_ipcw_0.75.RData")  # loads an object named `results`

# ---- 2. Reference ("true") marginal HR values ----

# generate_data_simple <- function(n = 1000, hx0 = 0.115, beta_x = log(0.75),
                                 #betax_z1 = log(0.5), seed = NULL) {
  #if (!is.null(seed)) set.seed(seed)
  #trt <- rbinom(n, 1, 0.5)
  #Z1 <- rbinom(n, 1, 0.5)
  
  #rate_x <- hx0 * exp(beta_x * trt + betax_z1 * Z1)
  #Tx <- rexp(n, rate = rate_x)
  #time <- Tx
  #status <- 1
  #data.frame(trt = trt, Z1 = Z1, time = time, status = status)
#}

#set.seed(2025)
#data <- generate_data_simple(n = 10000000, hx0 = 0.115, beta_x = log(0.75),
                             betax_z1 = log(2), seed = NULL)

#fit <- coxph(Surv(time, status)~trt, data=data)  
#fit$coefficients
                    
# The true log-HR used to simulate the data is log(HR_trt_event), but the
# *marginal* HR (what a correctly specified, unbiased estimator recovers)
# differs slightly across HR_x_z1 due to non-collapsibility. These
# reference values should be obtained from a large calibration run
# (e.g. n very large, no censoring) for the specific HR_trt_event used;
# update them here if you re-run the simulation with different settings.
ref_values <- c("0.5" = 0.768204, "1" = 0.75268, "2" = 0.770138)

# ---- 3. Compute relative bias (%) for each estimator ----
bias_data <- results %>%
  mutate(
    ref_value          = ref_values[as.character(HR_x_z1)],
    bias_unweighted    = (exp(unweighted.trt)    - ref_value) / ref_value * 100,
    bias_ipw_both      = (exp(ipw_both.trt)      - ref_value) / ref_value * 100,
    bias_ipw_trtonly   = (exp(ipw_trtonly.trt)   - ref_value) / ref_value * 100,
    bias_ipw_Z1only    = (exp(ipw_Z1only.trt)    - ref_value) / ref_value * 100,
    bias_standardized  = (exp(cox_adjusted)        - ref_value) / ref_value * 100
  )

# ---- 4. Aggregate mean bias + 95% CI by scenario, for one estimator ----
summarise_bias <- function(data, bias_col, method_name) {
  data %>%
    group_by(HR_trt_allo, HR_x_z1, HR_a_z1) %>%
    summarise(
      mean_bias = mean(.data[[bias_col]], na.rm = TRUE),
      se_bias   = sd(.data[[bias_col]], na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    ) %>%
    mutate(
      method = method_name,
      lower  = mean_bias - 1.96 * se_bias,
      upper  = mean_bias + 1.96 * se_bias
    )
}

bias_unweighted   <- summarise_bias(bias_data, "bias_unweighted",   "Cox (unweighted)")
bias_ipw_both     <- summarise_bias(bias_data, "bias_ipw_both",     "IPCW (treatment + Z1)")
bias_ipw_trtonly  <- summarise_bias(bias_data, "bias_ipw_trtonly",  "IPCW (treatment only)")
bias_ipw_Z1only   <- summarise_bias(bias_data, "bias_ipw_Z1only",   "IPCW (Z1 only)")
bias_standardized <- summarise_bias(bias_data, "bias_standardized", "Adjusted Cox (standardized)")

# ---- 5. Reusable plotting function (faceted bias plot) ----
make_bias_plot <- function(plot_data, colors, title) {
  
  plot_data <- plot_data %>%
    mutate(
      beta_trt_event = factor(paste0("HR=", HR_x_z1),
                              levels = c("HR=0.5", "HR=1", "HR=2")),
      beta_trt_allo  = factor(paste0("HR=", HR_trt_allo),
                              levels = c("HR=1", "HR=2", "HR=3")),
      outcome_label  = "Covariate effect on the outcome",
      allo_label     = "Treatment effect on alloSCT"
    )
  
  ggplot(plot_data, aes(x = HR_a_z1, y = mean_bias, color = method, fill = method)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray30", linewidth = 0.6) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, color = NA) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2) +
    ggh4x::facet_nested(outcome_label + beta_trt_event ~ allo_label + beta_trt_allo) +
    scale_x_continuous(breaks = c(0.2, 0.6, 1, 1.4, 2.2, 3), limits = c(0.2, 3)) +
    scale_y_continuous(breaks = seq(-10, 10, by = 5), limits = c(-10, 10)) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    theme_bw(base_size = 12) +
    labs(
      title    = title,
      subtitle = "Mean relative bias \u00b1 95% CI (10,000 simulations, n = 500)",
      x = "Covariate effect on alloSCT (HR)",
      y = "Mean Relative Bias (%)",
      color = "Method", fill = "Method"
    ) +
    theme(
      strip.text      = element_text(size = 9, face = "bold"),
      plot.title      = element_text(size = 14, hjust = 0.5, face = "bold"),
      plot.subtitle   = element_text(size = 11, hjust = 0.5),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
}

# ============================================================
# Figure 1 : Cox (unweighted) vs IPCW (full)
# ============================================================
plot_data_1 <- bind_rows(bias_unweighted, bias_ipw_both)

colors_1 <- c(
  "Cox (unweighted)"      = "#D55E00",
  "IPCW (treatment + Z1)" = "#0072B2"
)

p1 <- make_bias_plot(plot_data_1, colors_1,
                     "Unweighted Cox vs. Full IPCW")
print(p1)

# ============================================================
# Figure 2 : Cox (unweighted) vs IPCW (full, treatment-only, Z1-only)
# ============================================================
plot_data_2 <- bind_rows(bias_unweighted, bias_ipw_both,
                         bias_ipw_trtonly, bias_ipw_Z1only)

colors_2 <- c(
  "Cox (unweighted)"      = "#D55E00",
  "IPCW (treatment + Z1)" = "#0072B2",
  "IPCW (treatment only)" = "#E69F00",
  "IPCW (Z1 only)"        = "#009E73"
)

p2 <- make_bias_plot(plot_data_2, colors_2,
                     "Unweighted Cox vs. IPCW under different censoring models")
print(p2)

# ============================================================
# Figure 3 : Cox (unweighted) vs IPCW (full) vs Standardized adjusted Cox
# ============================================================
plot_data_3 <- bind_rows(bias_unweighted, bias_ipw_both, bias_standardized)

colors_3 <- c(
  "Cox (unweighted)"            = "#D55E00",
  "IPCW (treatment + Z1)"       = "#0072B2",
  "Adjusted Cox (standardized)" = "#009E73"
)

p3 <- make_bias_plot(plot_data_3, colors_3,
                     "Unweighted Cox vs. IPCW vs. Standardized Adjusted Cox")
print(p3)



