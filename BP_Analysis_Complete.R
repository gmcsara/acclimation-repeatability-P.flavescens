# Goal 0 -- Repeatability of boldness-related behaviour under different
# acclimation times, in Pomatoschistus flavescens.
# Martins-Cardoso, S. & Faria, A.M.
#
# Analyses performed:
#   (A) Original-manuscript tests (acclimation x trial interaction,
#       between-trial Spearman) re-run on the reanalysis variables
#   (B) coxme ICC for right-censored latency data
#   (C) rptR Gaussian repeatability on logit-transformed space use (TB2+TB3)
#   (D) Sensitivity check on the logit adjustment parameter (adj)
#   (E) Formal comparison of R across acclimation groups (bootstrap
#       overlap + Fisher's Z)
#   (F) Retrospective power analysis (simulation-based)
#
# Data file required: RawData21.xlsx (place in the same directory as this
# script). Figures are written to the "Figures_Behavioural_Processes"
# subfolder, created automatically if it does not exist.

# Plain-ASCII folder name: accented characters (e.g. "ã") can be garbled on
# Windows, so the figures end up in a folder with a different name.
fig_out <- "Figures_Behavioural_Processes"
dir.create(fig_out, showWarnings = FALSE)
cat("Figures will be saved in:", normalizePath(fig_out), "\n")

# ==== Setup and packages ====
library(tidyverse)
library(readxl)
library(lme4)
library(lmerTest)
library(rptR)
library(survival)
library(coxme)
library(car)
library(patchwork)

# ==== Figure style (shared by all figures) ====
# One definition used everywhere, so that every figure shows the traits in the
# same order, with the same names, the same colours and the same resolution.
trait_levels <- c("Latency to exit shelter", "Near-wall use", "Open-zone use")
accl_cols <- c("2 min" = "#e74c3c", "5 min" = "#3498db", "15 min" = "#2ecc71")
# Colours for 2/5/15 min are reserved for acclimation groups; other groupings
# (logit adj, sex) use greys so the two are never confused.
adj_cols <- c("adj=0.001" = "black", "adj=0.01" = "grey45", "adj=0.05" = "grey70")
sex_cols <- c("Without sex" = "black", "With sex" = "grey60")
theme_fig <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "white"),
        legend.position = "bottom")
FIG_DPI <- 600
save_fig <- function(file, plot = last_plot(), width, height) {
  ggsave(file.path(fig_out, file), plot,
         width = width, height = height, units = "cm", dpi = FIG_DPI)
}

# kableExtra is not used here; kbl()/kable_styling() are thin wrappers
# around knitr::kable() so the original "... %>% kbl() %>% kable_styling() %>% ..."
# pipelines still run and print plain tables to the console.
kbl           <- function(x, ...) knitr::kable(x, ...)
kable_styling <- function(x, ...) x
collapse_rows <- function(x, ...) x
row_spec      <- function(x, ...) x

# ==== Data import and variable construction ====
d <- as.data.frame(read_excel("RawData21.xlsx"))

d <- d %>%
  mutate(
    ID               = factor(ind_ID),
    Sex              = factor(Sex, levels = c(1, 2), labels = c("Male", "Female")),
    trial_number     = factor(trial_number),
    acclimation_time = factor(acclimation_time, levels = c(2, 5, 15)),
    SL               = as.numeric(SL),
    Latency_shelter  = as.numeric(Latency_shelter),
    Shelter          = as.numeric(Shelter),
    TB1              = as.numeric(TB1),
    TB2              = as.numeric(TB2),
    TB3              = as.numeric(TB3)
  )

T_MAX <- 600  # trial duration cap (seconds)

# Logit helper with adjustable floor/ceiling
lg <- function(p, adj = 0.001) {
  p <- pmin(pmax(p, adj), 1 - adj)
  log(p / (1 - p))
}

d <- d %>%
  mutate(
    SLz             = as.numeric(scale(SL)),

    # Emergence: a trial counts as an emergence only if the fish left the
    # shelter before the ceiling AND occupied at least one zone afterwards.
    # Fifteen trials were recorded with latencies of 599.3-599.9 s (the video
    # ended fractionally before the 600 s ceiling) and zero time in all four
    # zones: those fish never emerged, and are right-censored at T_MAX.
    # Genuine emergences reach at most 524.7 s, so the two groups are well
    # separated in the latency distribution.
    zone_sum        = Shelter + TB1 + TB2 + TB3,
    emerged         = zone_sum > 0 & Latency_shelter < T_MAX,
    surv_time       = ifelse(emerged, Latency_shelter, T_MAX),
    surv_event      = as.integer(emerged),

    space_prop      = (TB2 + TB3) / T_MAX,
    space_logit     = lg(space_prop, adj = 0.001),
    space_logit_01  = lg(space_prop, adj = 0.01),
    space_logit_05  = lg(space_prop, adj = 0.05),
    near_prop       = (Shelter + TB1) / T_MAX,
    near_logit      = lg(near_prop, adj = 0.001),
    near_logit_01   = lg(near_prop, adj = 0.01),
    near_logit_05   = lg(near_prop, adj = 0.05)
  )

d2  <- droplevels(subset(d, acclimation_time == "2"))
d5  <- droplevels(subset(d, acclimation_time == "5"))
d15 <- droplevels(subset(d, acclimation_time == "15"))

# ---- Data summary ----
tibble(
  Statistic = c(
    "N observations", "N individuals",
    "Censored latencies (never emerged)",
    "Trials with zero open-zone use"
  ),
  Value = c(
    nrow(d),
    nlevels(d$ID),
    sprintf("%d (%.0f%%)", sum(d$surv_event == 0),
            100 * mean(d$surv_event == 0)),
    sprintf("%d (%.0f%%)", sum(d$space_prop == 0),
            100 * mean(d$space_prop == 0))
  )
) %>%
  kbl(caption = "Data overview") %>%
  kable_styling(bootstrap_options = c("striped","hover"), full_width = FALSE)

(table(d$acclimation_time) / 3) %>%
  as.data.frame() %>%
  setNames(c("Acclimation (min)", "N individuals")) %>%
  kbl(caption = "Sample size per acclimation group") %>%
  kable_styling(bootstrap_options = c("striped","hover"), full_width = FALSE)

# ==== Descriptive statistics ====
desc_stats <- function(x) {
  c(n    = sum(!is.na(x)),
    mean = round(mean(x, na.rm = TRUE), 2),
    SE   = round(sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))), 2),
    min  = round(min(x, na.rm = TRUE), 2),
    max  = round(max(x, na.rm = TRUE), 2))
}

vars_raw   <- c("Latency_shelter","Shelter","TB1","TB2","TB3")
var_labels <- c("Latency to exit shelter (s)","Time in shelter (s)",
                "Time in zone 1 (s)","Time in zone 2 (s)","Time in zone 3 (s)")

# ---- Overall (Table 1) ----
desc_overall <- do.call(rbind, lapply(seq_along(vars_raw), function(i) {
  s <- desc_stats(d[[vars_raw[i]]])
  data.frame(Variable = var_labels[i], t(s))
}))

desc_overall %>%
  kbl(caption = "Table 1. Overall descriptive statistics (n = 31 individuals, 3 trials)",
      digits = 2) %>%
  kable_styling(bootstrap_options = c("striped","hover"))

write.csv(desc_overall, "Table1_descriptive_overall.csv", row.names = FALSE)

# ---- By acclimation group (Table 2) ----
desc_by_accl <- do.call(rbind, lapply(seq_along(vars_raw), function(i) {
  do.call(rbind, lapply(c("2","5","15"), function(ac) {
    dd <- d[d$acclimation_time == ac, ]
    s  <- desc_stats(dd[[vars_raw[i]]])
    data.frame(Variable = var_labels[i],
               Acclimation = paste0(ac, " min"), t(s))
  }))
}))

desc_by_accl %>%
  kbl(caption = "Table 2. Descriptive statistics by acclimation group", digits = 2) %>%
  kable_styling(bootstrap_options = c("striped","hover")) %>%
  collapse_rows(columns = 1, valign = "top")

write.csv(desc_by_accl, "Table2_descriptive_by_acclimation.csv", row.names = FALSE)

# ---- By trial, new composite variables (Figure descriptive) ----
d_trial_new <- d %>%
  mutate(
    `Latency to exit shelter (s)` = Latency_shelter,
    `Near-wall use (s)`           = Shelter + TB1,
    `Open-zone use (s)`           = TB2 + TB3
  ) %>%
  select(ID, acclimation_time, trial_number,
         `Latency to exit shelter (s)`,
         `Near-wall use (s)`,
         `Open-zone use (s)`) %>%
  pivot_longer(cols = -c(ID, acclimation_time, trial_number),
               names_to = "Trait", values_to = "Value") %>%
  mutate(Trait = factor(Trait, levels = c(
    "Latency to exit shelter (s)",
    "Near-wall use (s)",
    "Open-zone use (s)"
  )))

trial_summary_new <- d_trial_new %>%
  group_by(Trait, acclimation_time, trial_number) %>%
  summarise(
    n    = sum(!is.na(Value)),
    mean = round(mean(Value, na.rm = TRUE), 2),
    SE   = round(sd(Value, na.rm = TRUE) / sqrt(sum(!is.na(Value))), 2),
    .groups = "drop"
  )

write.csv(trial_summary_new, "Table_descriptive_by_trial_new_vars.csv", row.names = FALSE)

trial_summary_new %>%
  kbl(caption = "Descriptive statistics (mean +/- SE) by trial and acclimation group, new composite variables", digits = 2) %>%
  kable_styling(bootstrap_options = c("striped","hover")) %>%
  collapse_rows(columns = 1:2, valign = "top")
ggplot(trial_summary_new,
       aes(x = trial_number, y = mean, colour = acclimation_time, group = acclimation_time)) +
  geom_errorbar(aes(ymin = mean - SE, ymax = mean + SE),
                width = 0.12, linewidth = 0.5,
                position = position_dodge(width = 0.2)) +
  geom_point(size = 2.5, position = position_dodge(width = 0.2)) +
  geom_line(position = position_dodge(width = 0.2), linewidth = 0.5) +
  facet_wrap(~ Trait, scales = "free_y") +
  scale_colour_manual(values = setNames(accl_cols, c("2", "5", "15")),
                       labels = names(accl_cols),
                       name = "Acclimation time") +
  labs(x = "Trial", y = "Mean +/- SE") +
  theme_fig

save_fig("Figure_4.png", width = 26, height = 11)

# ==== Boldness (latency to exit shelter): Cox mixed model + profile-likelihood ICC ====
# McCune et al. (2025) tutorial functions (R/rptRsurv.R; Nakagawa et al.).
# Two helpers for a Cox mixed model (coxme):
#   coxme_icc_ci(): 95% CI for the ICC, by profile likelihood.
#   coxme_pval():   p-value for the random effect, by likelihood-ratio test.

# coxme_icc_ci(): confidence interval for the repeatability (ICC).
#
# IDEA (profile likelihood): the model's log-likelihood is highest at the
# fitted among-individual variance. If we FORCE the variance to other values
# and refit, the likelihood drops. The 95% CI is the set of variance values
# whose likelihood is "not significantly worse" than the best fit -- i.e. where
# the likelihood-ratio statistic stays below chi-square(0.95, df = 1) = 3.84.
# We map those variance bounds onto the ICC scale at the end.
coxme_icc_ci <- function(model, upper.multiplier = 10) {

  # Safety checks: must be a coxme model with exactly ONE random effect.
  if (!inherits(model, "coxme")) stop("need a coxme model")
  if (length(summary(model)$random$variance) > 1) stop("only one random effect")

  cut <- 100                                          # grid points on EACH side
  var_point <- summary(model)$random$variance         # the fitted ID variance

  # Build a grid of candidate variances to test:
  #   - lower half: from ~0 up to the fitted value
  #   - upper half: from the fitted value up to (fitted x upper.multiplier)
  # upper.multiplier controls how far ABOVE the estimate we search for the
  # upper CI bound (bigger = look further up).
  estvar <- c(seq(1e-14, var_point, length = cut),
              seq(var_point, var_point * upper.multiplier, length = cut + 1)[-1])

  # For each candidate variance, refit the model with the variance FIXED at that
  # value (vfixed) and record 2 x the log-likelihood gain over the null.
  loglik <- double(cut * 2)
  for (i in seq_len(cut * 2)) {
    tfit <- update(model, vfixed = estvar[i])         # refit, variance held fixed
    loglik[i] <- 2 * diff(tfit$loglik)[1]
  }

  # temp = how much WORSE each fixed-variance fit is than the best (free) fit.
  # This is the likelihood-ratio statistic as a function of the variance.
  temp <- as.numeric(2 * diff(model$loglik)[1]) - loglik

  # Find the variance where temp crosses 3.84 (the 95% threshold), separately
  # below the estimate (lower bound) and above it (upper bound). approx() reads
  # off the variance at that crossing. NOTE: approx() interpolates only — if the
  # curve never reaches 3.84 within the grid (flat likelihood), it returns NA.
  lower <- approx(temp[1:cut],               sqrt(estvar[1:cut]),               qchisq(.95, 1))$y
  upper <- approx(temp[(cut + 1):(2 * cut)], sqrt(estvar[(cut + 1):(2 * cut)]), qchisq(.95, 1))$y

  # Convert the variance bounds to the ICC scale: ICC = var / (var + pi^2/6).
  # (pi^2/6 is the residual variance of the Cox / extreme-value latent scale.)
  ICC_lower <- lower^2 / (lower^2 + pi^2 / 6); if (is.na(ICC_lower)) ICC_lower <- 0
  ICC_point <- var_point / (var_point + pi^2 / 6)
  ICC_upper <- upper^2 / (upper^2 + pi^2 / 6)         # NA when the profile is flat (var ~ 0)

  c(lower = ICC_lower, ICC = ICC_point, upper = ICC_upper)
}

# coxme_pval(): is the among-individual variance > 0? (p-value)
#
# Compares the mixed model (with the ID random effect) to the SAME model
# without it (a plain coxph), via a likelihood-ratio test. Because a variance
# cannot be negative, the test is on a boundary, so we halve the usual p-value
# (the standard 50:50 mixture correction).
coxme_pval <- function(model, data) {
  fixed_formula <- as.formula(model$formulaList$fixed)   # the model's fixed part only
  fit <- survival::coxph(fixed_formula, data = data)     # refit WITHOUT (1|ID)
  # LR statistic = 2 x (loglik_with_RE - loglik_without_RE)
  lr  <- 2 * (as.numeric(model$loglik["Integrated"]) -
              as.numeric(fit$loglik[length(fit$loglik)]))
  setNames(0.5 * pchisq(lr, df = 1, lower.tail = FALSE), "p_LRT")  # boundary-corrected
}

# icc_coxme(): convenience wrapper used below for each acclimation group.
# Fits the Cox mixed model and returns ICC + 95% profile-likelihood CI + LRT p.
#
# Why the .tte_dat / bquote dance: coxme_icc_ci() calls update() internally,
# which re-runs the original model call (it contains `data = ...`). update()
# must be able to FIND that data object by name when it re-evaluates. So we
# park the data in a fixed global name (.tte_dat) and bake that exact name into
# the coxme call with bquote(), guaranteeing update() can re-find it.
icc_coxme <- function(dd, fixed_str, label = "") {
  assign(".tte_dat", droplevels(dd), envir = .GlobalEnv)        # expose data globally
  frm <- as.formula(paste0("Surv(surv_time, surv_event) ~ ", fixed_str, " + (1|ID)"))
  m   <- eval(bquote(coxme(.(frm), data = .tte_dat)))           # fit; call refers to .tte_dat
  v   <- as.numeric(summary(m)$random$variance)                # fitted ID variance

  # ICC + 95% CI (profile likelihood). upper.multiplier = 50 searches far enough
  # up to bound small-but-non-zero variances; tryCatch guards against failures.
  ci  <- tryCatch(coxme_icc_ci(m, upper.multiplier = 50),
                  error = function(e) c(lower = NA, ICC = v / (v + pi^2 / 6), upper = NA))
  # p-value for the random effect (LRT); NA if the test cannot be computed.
  p   <- tryCatch(as.numeric(coxme_pval(m, .tte_dat)), error = function(e) NA_real_)

  data.frame(Trait = "Boldness (latency to exit shelter)", Group = label,
             ICC = round(ci["ICC"], 3),
             CI_lo = round(ci["lower"], 3), CI_hi = round(ci["upper"], 3),
             p_LRT = round(p, 4), sigma2_ID = round(v, 4), row.names = NULL)
}
bold_overall <- icc_coxme(d,   "trial_number + SLz + acclimation_time", label = "Overall")
bold_2min    <- icc_coxme(d2,  "trial_number + SLz", label = "2 min")
bold_5min    <- icc_coxme(d5,  "trial_number + SLz", label = "5 min")
bold_15min   <- icc_coxme(d15, "trial_number + SLz", label = "15 min")

bold_tab <- rbind(bold_overall, bold_2min, bold_5min, bold_15min)
write.csv(bold_tab, "Table_boldness_ICC_coxme.csv", row.names = FALSE)
bold_tab %>%
  kbl(caption = "Boldness (latency to exit shelter): coxme ICC with 95% profile-likelihood CI and LRT p-value (McCune et al. 2025)",
      digits = 3) %>%
  kable_styling(bootstrap_options = c("striped","hover")) %>%
  row_spec(1, bold = TRUE)

# ==== Space use (open zone): rptR Gaussian ====
# Repeatability (R) of a Gaussian trait via rptR, returned as one tidy row.
# rptR fits the Gaussian LMM  y ~ fixed effects + (1|ID)  INTERNALLY with
# lme4::lmer (no direct lmer() call needed), then R = V_among / (V_among+V_resid).
rptR_row <- function(dd, fixed_str, label, yvar = "space_logit", seed = 42) {
  dd  <- droplevels(dd)
  set.seed(seed)
  frm <- as.formula(paste0(yvar, " ~ ", fixed_str, " + (1|ID)"))
  r   <- rpt(frm, grname = "ID", data = dd,
             datatype = "Gaussian", nboot = 1000, npermut = 1000,
             verbose = FALSE)
  data.frame(
    Trait  = yvar, Group = label,
    R      = round(r$R$ID[1], 3),
    CI_lo  = round(quantile(r$R_boot$ID, 0.025, na.rm = TRUE), 3),
    CI_hi  = round(quantile(r$R_boot$ID, 0.975, na.rm = TRUE), 3),
    P_LRT  = round(r$P$LRT_P[1], 4),
    P_perm = round(r$P$P_permut[1], 4),
    row.names = NULL
  )
}

# ---- Open zones ----
space_overall <- rptR_row(d,   "trial_number + SLz + acclimation_time", "Overall")
space_2min    <- rptR_row(d2,  "trial_number + SLz", "2 min")
space_5min    <- rptR_row(d5,  "trial_number + SLz", "5 min")
space_15min   <- rptR_row(d15, "trial_number + SLz", "15 min")

space_tab <- rbind(space_overall, space_2min, space_5min, space_15min)
write.csv(space_tab, "Table_space_use_rptR.csv", row.names = FALSE)
space_tab %>%
  kbl(caption = "Space use (open zone): rptR Gaussian R",
      digits = 3) %>%
  kable_styling(bootstrap_options = c("striped","hover")) %>%
  row_spec(1, bold = TRUE) %>%
  row_spec(which(space_tab$P_perm < 0.05), background = "#f0f8ff")

# ---- Near-wall zone — secondary check ----
near_overall <- rptR_row(d,   "trial_number + SLz + acclimation_time", "Overall", yvar = "near_logit")
near_2min    <- rptR_row(d2,  "trial_number + SLz", "2 min",  yvar = "near_logit")
near_5min    <- rptR_row(d5,  "trial_number + SLz", "5 min",  yvar = "near_logit")
near_15min   <- rptR_row(d15, "trial_number + SLz", "15 min", yvar = "near_logit")

near_tab <- rbind(near_overall, near_2min, near_5min, near_15min)
near_tab$Trait <- "Space use (near wall)"
write.csv(near_tab, "Table_near_zone_rptR.csv", row.names = FALSE)
near_tab %>%
  kbl(caption = "Near-wall zone: rptR Gaussian R",
      digits = 3) %>%
  kable_styling(bootstrap_options = c("striped","hover")) %>%
  row_spec(1, bold = TRUE)

# ==== Variance components (Table 4) ====
fmt <- function(v) sprintf("%.3f (%.3f)", v, sqrt(v))    # variance (SD)

# Each fitted model is also stored in `vc_fits`, so that the fixed effects
# (trial number, SL) reported in the text come from these same models.
vc_fits <- list()

# V_I / V_R from a Gaussian LMM (open-zone, near-wall)
vc_lmer <- function(yvar, dd, label = "") {
  m  <- lmer(as.formula(paste0(yvar, " ~ trial_number + SLz + (1|ID)")),
             data = droplevels(dd), REML = TRUE)
  vc_fits[[paste(yvar, label)]] <<- m
  vc <- as.data.frame(VarCorr(m))
  c(VI = vc$vcov[vc$grp == "ID"], VR = vc$vcov[vc$grp == "Residual"])
}
# V_I / V_R from the Cox mixed model (boldness): V_R fixed at pi^2/6
vc_cox <- function(dd, label = "") {
  m <- coxme(Surv(surv_time, surv_event) ~ trial_number + SLz + (1|ID),
             data = droplevels(dd))
  vc_fits[[paste("latency", label)]] <<- m
  c(VI = as.numeric(VarCorr(m)$ID), VR = pi^2 / 6)
}

grps <- list("2 min" = d2, "5 min" = d5, "15 min" = d15)
build <- function(label, fn) do.call(rbind, lapply(names(grps), function(g) {
  v <- fn(grps[[g]], g)
  data.frame(Variable = label, Acclimation = g,
             `Among-individual (V_I)`  = fmt(v["VI"]),
             `Within-individual (V_R)` = fmt(v["VR"]),
             check.names = FALSE, row.names = NULL)
}))

table4 <- rbind(
  build("Boldness (latency to exit shelter)", function(dd, g) vc_cox(dd, g)),
  build("Near-wall use",                     function(dd, g) vc_lmer("near_logit",  dd, g)),
  build("Open-zone use",                     function(dd, g) vc_lmer("space_logit", dd, g))
)
write.csv(table4, "Table4_variance_components.csv", row.names = FALSE)

table4 %>%
  kbl(caption = "Table 4. Variance components (V_I among-individual, V_R within-individual; SD in parentheses). Boldness: coxme (V_R = pi^2/6, latent scale). Space use: lmer. Scales not comparable across traits.") %>%
  kable_styling(bootstrap_options = c("striped", "hover")) %>%
  collapse_rows(columns = 1, valign = "top")

# ---- Fixed effects (trial number, SL) from the same models ----
# coxme: coefficients on the log-hazard scale (positive = earlier emergence).
# lmer (lmerTest): coefficients on the logit scale, Satterthwaite p-values.
fixed_row <- function(key) {
  m <- vc_fits[[key]]
  parts <- strsplit(key, " ", fixed = TRUE)[[1]]
  if (inherits(m, "coxme")) {
    b <- fixef(m); se <- sqrt(diag(vcov(m))); p <- 2 * pnorm(-abs(b / se))
  } else {
    cf <- summary(m)$coefficients[-1, , drop = FALSE]
    b <- cf[, "Estimate"]; se <- cf[, "Std. Error"]; p <- cf[, "Pr(>|t|)"]
  }
  data.frame(Trait = parts[1], Group = paste(parts[-1], collapse = " "),
             Term = names(b), beta = round(b, 3), SE = round(se, 3),
             p = round(p, 3), row.names = NULL)
}
fixed_tab <- do.call(rbind, lapply(names(vc_fits), fixed_row))
write.csv(fixed_tab, "Table_fixed_effects.csv", row.names = FALSE)
fixed_tab %>%
  kbl(caption = "Fixed effects (trial number, SL) from the variance-component models, per acclimation group.")

# ==== Sensitivity check: logit adjustment parameter ====
run_sensitivity <- function(dd, fixed_str, label, yvar) {
  dd <- droplevels(dd)
  set.seed(42)
  frm <- as.formula(paste0(yvar, " ~ ", fixed_str, " + (1|ID)"))
  r   <- rpt(frm, grname = "ID", data = dd,
             datatype = "Gaussian", nboot = 1000, npermut = 0,
             verbose = FALSE)
  data.frame(Group = label,
             R     = round(r$R$ID[1], 3),
             CI_lo = round(quantile(r$R_boot$ID, 0.025, na.rm = TRUE), 3),
             CI_hi = round(quantile(r$R_boot$ID, 0.975, na.rm = TRUE), 3),
             row.names = NULL)
}

# Both space-use measures, each at three logit floors (adj = 0.001 / 0.01 / 0.05)
adj_sets <- list(
  "Near-wall" = c("adj=0.001" = "near_logit",
                  "adj=0.01"  = "near_logit_01",
                  "adj=0.05"  = "near_logit_05"),
  "Open-zone" = c("adj=0.001" = "space_logit",
                  "adj=0.01"  = "space_logit_01",
                  "adj=0.05"  = "space_logit_05")
)
sens_tab <- do.call(rbind, lapply(names(adj_sets), function(tr) {
  vars <- adj_sets[[tr]]
  do.call(rbind, lapply(seq_along(vars), function(i) {
    out <- rbind(
      run_sensitivity(d2,  "trial_number + SLz", "2 min",  vars[i]),
      run_sensitivity(d5,  "trial_number + SLz", "5 min",  vars[i]),
      run_sensitivity(d15, "trial_number + SLz", "15 min", vars[i]))
    out$Trait <- tr
    out$adj   <- names(vars)[i]
    out
  }))
}))
sens_tab <- sens_tab[, c("Trait", "adj", "Group", "R", "CI_lo", "CI_hi")]
write.csv(sens_tab, "Table_sensitivity_logit_adj.csv", row.names = FALSE)
sens_tab %>%
  kbl(caption = "Sensitivity of R to the logit adjustment parameter, for both space-use measures (near-wall and open-zone).",
      digits = 3) %>%
  kable_styling(bootstrap_options = c("striped","hover")) %>%
  collapse_rows(columns = 1:2, valign = "top")
sens_plot_data <- sens_tab %>%
  mutate(Group = factor(Group, levels = c("2 min","5 min","15 min")),
         adj   = factor(adj, levels = c("adj=0.001","adj=0.01","adj=0.05")),
         Trait = factor(paste(Trait, "use"), levels = trait_levels))

ggplot(sens_plot_data, aes(x = Group, y = R, colour = adj, group = adj)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_errorbar(aes(ymin = CI_lo, ymax = CI_hi), width = 0.12, linewidth = 0.5,
                position = position_dodge(width = 0.3)) +
  geom_point(size = 2.5, position = position_dodge(width = 0.3)) +
  geom_line(position  = position_dodge(width = 0.3), linewidth = 0.5) +
  facet_wrap(~ Trait) +
  scale_colour_manual(values = adj_cols, name = "Logit adj") +
  scale_y_continuous(limits = c(-0.1, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = "Acclimation time", y = "Repeatability R") +
  theme_fig

save_fig("Figure_S1.png", width = 20, height = 12)

# ==== Formal comparison of R across acclimation groups ====

# ---- Bootstrap distributions ----
get_rpt_boot <- function(dd, fixed_str, label,
                          yvar = "space_logit", seed = 42) {
  dd <- droplevels(dd)
  set.seed(seed)
  frm <- as.formula(paste0(yvar, " ~ ", fixed_str, " + (1|ID)"))
  r   <- rpt(frm, grname = "ID", data = dd,
             datatype = "Gaussian", nboot = 1000, npermut = 0,
             verbose = FALSE)
  list(label  = label,
       R      = r$R$ID[1],
       R_boot = r$R_boot$ID,
       R_se   = sd(r$R_boot$ID, na.rm = TRUE))
}

boot_2min  <- get_rpt_boot(d2,  "trial_number + SLz", "2 min")
boot_5min  <- get_rpt_boot(d5,  "trial_number + SLz", "5 min")
boot_15min <- get_rpt_boot(d15, "trial_number + SLz", "15 min")
boot_2min_n  <- get_rpt_boot(d2,  "trial_number + SLz", "2 min",  yvar = "near_logit")
boot_5min_n  <- get_rpt_boot(d5,  "trial_number + SLz", "5 min",  yvar = "near_logit")
boot_15min_n <- get_rpt_boot(d15, "trial_number + SLz", "15 min", yvar = "near_logit")
boot_df <- rbind(
  data.frame(Trait = "Near-wall", Group = "2 min",  R = boot_2min_n$R_boot),
  data.frame(Trait = "Near-wall", Group = "5 min",  R = boot_5min_n$R_boot),
  data.frame(Trait = "Near-wall", Group = "15 min", R = boot_15min_n$R_boot),
  data.frame(Trait = "Open-zone", Group = "2 min",  R = boot_2min$R_boot),
  data.frame(Trait = "Open-zone", Group = "5 min",  R = boot_5min$R_boot),
  data.frame(Trait = "Open-zone", Group = "15 min", R = boot_15min$R_boot)
) %>% filter(!is.na(R)) %>%
  mutate(Group = factor(Group, levels = c("2 min","5 min","15 min")),
         Trait = factor(paste(Trait, "use"), levels = trait_levels))

ggplot(boot_df, aes(x = R, fill = Group, colour = Group)) +
  geom_density(alpha = 0.35, linewidth = 0.6) +
  facet_wrap(~ Trait) +
  scale_fill_manual(values   = accl_cols, name = "Acclimation time") +
  scale_colour_manual(values = accl_cols, name = "Acclimation time") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  scale_x_continuous(limits = c(-0.1, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = "Repeatability R", y = "Bootstrap density") +
  theme_fig

save_fig("Figure_S3.png", width = 20, height = 12)

# ---- Bootstrap overlap P(R~A~ > R~B~) ----
boot_overlap <- function(boot_a, boot_b) {
  ba <- boot_a$R_boot; bb <- boot_b$R_boot
  n  <- min(length(ba), length(bb))
  ba <- ba[!is.na(ba)][1:n]; bb <- bb[!is.na(bb)][1:n]
  data.frame(Group_A  = boot_a$label, Group_B = boot_b$label,
             R_A      = round(boot_a$R, 3), R_B = round(boot_b$R, 3),
             Delta_R  = round(boot_a$R - boot_b$R, 3),
             P_A_gt_B = round(mean(ba > bb, na.rm = TRUE), 3))
}

overlap_trio <- function(b2, b5, b15, trait)
  cbind(Trait = trait,
        rbind(boot_overlap(b2, b5), boot_overlap(b2, b15), boot_overlap(b5, b15)))

overlap_tab <- rbind(
  overlap_trio(boot_2min_n, boot_5min_n, boot_15min_n, "Near-wall"),
  overlap_trio(boot_2min,   boot_5min,   boot_15min,   "Open-zone")
)
write.csv(overlap_tab, "Table_bootstrap_overlap.csv", row.names = FALSE)

overlap_tab %>%
  kbl(caption = "Bootstrap overlap: P(R_A > R_B), near-wall and open-zone.",
      digits = 3) %>%
  kable_styling(bootstrap_options = c("striped","hover")) %>%
  collapse_rows(columns = 1, valign = "top")

# ---- Fisher's Z pairwise test ----
fisher_z <- function(r) 0.5 * log((1 + r) / (1 - r))

fisher_compare <- function(boot_a, boot_b) {
  Za     <- fisher_z(boot_a$R); Zb <- fisher_z(boot_b$R)
  se_Za  <- boot_a$R_se / (1 - boot_a$R^2)
  se_Zb  <- boot_b$R_se / (1 - boot_b$R^2)
  se_diff <- sqrt(se_Za^2 + se_Zb^2)
  z_stat  <- (Za - Zb) / se_diff
  p_val   <- 2 * pnorm(abs(z_stat), lower.tail = FALSE)
  data.frame(Group_A = boot_a$label, Group_B = boot_b$label,
             R_A = round(boot_a$R, 3), R_B = round(boot_b$R, 3),
             Z_A = round(Za, 3), Z_B = round(Zb, 3),
             Z_stat = round(z_stat, 3), P_two_tailed = round(p_val, 4))
}

fisher_trio <- function(b2, b5, b15, trait)
  cbind(Trait = trait,
        rbind(fisher_compare(b2, b5), fisher_compare(b2, b15), fisher_compare(b5, b15)))

fisher_tab <- rbind(
  fisher_trio(boot_2min_n, boot_5min_n, boot_15min_n, "Near-wall"),
  fisher_trio(boot_2min,   boot_5min,   boot_15min,   "Open-zone")
)
write.csv(fisher_tab, "Table_fisherZ_comparison.csv", row.names = FALSE)

fisher_tab %>%
  kbl(caption = "Fisher's Z pairwise comparison of R, near-wall and open-zone.",
      digits = 3) %>%
  kable_styling(bootstrap_options = c("striped","hover")) %>%
  collapse_rows(columns = 1, valign = "top")

# ==== Acclimation × trial interaction (LRT) ====
# --- Gaussian (logit) space-use measures: lmer interaction LRT ---
get_lrt <- function(yvar, label) {
  frm_add <- as.formula(paste0(yvar,
    " ~ trial_number + acclimation_time + SLz + (1|ID)"))
  frm_int <- as.formula(paste0(yvar,
    " ~ trial_number * acclimation_time + SLz + (1|ID)"))
  m_add <- lmer(frm_add, REML = FALSE, data = d)
  m_int <- lmer(frm_int, REML = FALSE, data = d)
  a     <- anova(m_add, m_int)                       # LRT: additive vs interaction
  data.frame(Variable  = label,
             AIC_add   = round(AIC(m_add), 1),
             AIC_int   = round(AIC(m_int), 1),
             Delta_AIC = round(AIC(m_int) - AIC(m_add), 1),
             Chi_sq    = round(a$Chisq[2], 3),
             df        = a$Df[2],
             P         = round(a$`Pr(>Chisq)`[2], 4))
}

# --- Boldness latency (censored): coxme survival interaction LRT ---
get_lrt_surv <- function(label) {
  m_add <- coxme(Surv(surv_time, surv_event) ~ trial_number + acclimation_time + SLz + (1|ID), data = d)
  m_int <- coxme(Surv(surv_time, surv_event) ~ trial_number * acclimation_time + SLz + (1|ID), data = d)
  a     <- anova(m_add, m_int)                       # anova.coxme: p-column is "P(>|Chi|)"
  pcol  <- grep("^P", colnames(a), value = TRUE)[1]  # grab it by pattern (robust)
  data.frame(Variable  = label,
             AIC_add   = NA, AIC_int = NA, Delta_AIC = NA,   # AIC not compared for coxme
             Chi_sq    = round(a$Chisq[2], 3),
             df        = a$Df[2],
             P         = round(a[[pcol]][2], 4))
}

lrt_tab <- rbind(
  get_lrt_surv("Boldness (latency to exit shelter, survival)"),
  get_lrt("near_logit",  "Near-wall use"),
  get_lrt("space_logit", "Open-zone use")
)
write.csv(lrt_tab, "Table_LRT_interaction.csv", row.names = FALSE)
lrt_tab %>%
  kbl(caption = "Table 6. LRT: additive vs. interaction model (acclimation × trial), reanalysis variables. Boldness via coxme (survival); space use via lmer.") %>%
  kable_styling(bootstrap_options = c("striped","hover"))

# ==== Pairwise Spearman correlations between trials ====
spearman_by_group <- function(dd, var, var_label, group_label) {
  dw <- dd %>%
    select(ID, trial_number, all_of(var)) %>%
    pivot_wider(names_from = trial_number, values_from = all_of(var),
                names_prefix = "T")
  pairs <- list(c("T1","T2"), c("T1","T3"), c("T2","T3"))
  do.call(rbind, lapply(pairs, function(p) {
    dp   <- dw %>% select(all_of(p)) %>% drop_na()
    test <- cor.test(dp[[1]], dp[[2]], method = "spearman", exact = FALSE)
    data.frame(Variable    = var_label,
               Acclimation = group_label,
               Pair        = paste(p, collapse = " vs "),
               n           = nrow(dp),
               rho         = round(test$estimate, 3),
               P           = round(test$p.value, 4))
  }))
}

# Reanalysis variables: boldness (latency, capped at T_MAX), open-zone use and
# near-wall use. Spearman is rank-based, so logit vs proportion gives the same rho
# (and the cap on latency only ties the censored fish, which is appropriate here).
spearman_vars <- list(
  list(col = "surv_time",   label = "Boldness (latency to exit shelter)"),
  list(col = "near_logit",  label = "Near-wall use"),
  list(col = "space_logit", label = "Open-zone use")
)

spearman_tab <- do.call(rbind, lapply(spearman_vars, function(v) {
  rbind(
    spearman_by_group(d2,  v$col, v$label, "2 min"),
    spearman_by_group(d5,  v$col, v$label, "5 min"),
    spearman_by_group(d15, v$col, v$label, "15 min")
  )
}))
write.csv(spearman_tab, "Table_spearman_by_acclimation.csv", row.names = FALSE)
spearman_tab %>%
  mutate(sig = ifelse(P < 0.05, "*", "")) %>%
  kbl(caption = "Table 5. Pairwise Spearman correlations between trials by acclimation group") %>%
  kable_styling(bootstrap_options = c("striped","hover")) %>%
  collapse_rows(columns = 1:2, valign = "top") %>%
  row_spec(which(spearman_tab$P < 0.05), background = "#f0f8ff")

# ==== Convergent validity: boldness × space use ====
agg <- d %>%
  mutate(log_lat = log(surv_time)) %>%
  group_by(ID, acclimation_time) %>%
  summarise(mean_log_lat = mean(log_lat,     na.rm = TRUE),
            mean_space   = mean(space_logit, na.rm = TRUE),
            mean_near    = mean(near_logit,  na.rm = TRUE),
            .groups = "drop")

ct_open <- cor.test(agg$mean_log_lat, agg$mean_space, method = "pearson")
ct_near <- cor.test(agg$mean_log_lat, agg$mean_near,  method = "pearson")

conv <- function(ct, pair)
  data.frame(Pair = pair,
             r     = round(ct$estimate, 3),
             CI_lo = round(ct$conf.int[1], 3),
             CI_hi = round(ct$conf.int[2], 3),
             P     = round(ct$p.value, 4),
             row.names = NULL)

rbind(
  conv(ct_near, "Boldness (latency to exit shelter) × Near-wall use"),
  conv(ct_open, "Boldness (latency to exit shelter) × Open-zone use")
) %>%
  kbl(caption = sprintf("Among-individual correlations between boldness and space-use measures (n = %d individuals).", nrow(agg))) %>%
  kable_styling(bootstrap_options = c("striped","hover"), full_width = FALSE)

# ==== Combined repeatability summary and forest plot ====
rep_summary <- rbind(
  bold_tab %>%
    transmute(Trait = "Boldness (latency to exit shelter) [coxme ICC]",
              Group, R = ICC, CI_lo, CI_hi, P_LRT = p_LRT, P_perm = NA),
  near_tab %>%
    transmute(Trait = "Space use – near wall [rptR]",
              Group, R, CI_lo, CI_hi, P_LRT, P_perm),
  space_tab %>%
    transmute(Trait = "Space use – open zone [rptR]",
              Group, R, CI_lo, CI_hi, P_LRT, P_perm)
)
rep_summary$Group <- factor(rep_summary$Group,
                             levels = c("Overall","2 min","5 min","15 min"))
write.csv(rep_summary, "Table_repeatability_combined_reanalysis.csv",
          row.names = FALSE)

rep_summary %>%
  mutate(across(where(is.numeric), ~round(.x, 3))) %>%
  kbl(caption = "Combined repeatability summary (reanalysis)") %>%
  kable_styling(bootstrap_options = c("striped","hover")) %>%
  collapse_rows(columns = 1, valign = "top") %>%
  row_spec(which(rep_summary$Group == "Overall"), bold = TRUE)
fig_data <- rep_summary %>%
  filter(Group != "Overall") %>%
  mutate(
    Trait = factor(Trait, levels = c(
      "Boldness (latency to exit shelter) [coxme ICC]",
      "Space use – near wall [rptR]",
      "Space use – open zone [rptR]"
    )),
    sig = ifelse(!is.na(P_perm) & P_perm < 0.05, "p < 0.05", "ns")
  )

ggplot(fig_data, aes(x = Group, y = R, colour = Group)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_errorbar(aes(ymin = CI_lo, ymax = CI_hi),
                width = 0.15, linewidth = 0.5) +
  geom_point(aes(shape = sig), size = 3) +
  scale_colour_manual(values = accl_cols, guide = "none") +
  scale_shape_manual(values = c("p < 0.05" = 19, "ns" = 1), name = NULL) +
  scale_y_continuous(limits = c(-0.05, 1), breaks = seq(0, 1, 0.25)) +
  facet_wrap(~ Trait, ncol = 3,
             labeller = labeller(Trait = c(
               "Boldness (latency to exit shelter) [coxme ICC]" = "Latency to exit shelter",
               "Space use – open zone [rptR]" = "Open-zone use",
               "Space use – near wall [rptR]" = "Near-wall use"
             ))) +
  labs(x = "Acclimation time",
       y = "Repeatability (R / ICC) with 95% CI") +
  theme_fig

save_fig("Figure_3.png", width = 22, height = 10)

# ==== Retrospective power analysis ====

# ---- 12A — Single-group power: detect R > 0 ----
simulate_power_single <- function(R_true, n_ind, n_trial = 3,
                                  nsim = 500, seed = 42) {
  set.seed(seed)
  sig_count <- 0
  V_among   <- R_true
  V_within  <- 1 - R_true

  for (s in seq_len(nsim)) {
    id_effects <- rnorm(n_ind, 0, sqrt(V_among))
    sim_d <- do.call(rbind, lapply(seq_len(n_ind), function(i) {
      data.frame(ID           = factor(i),
                 trial_number = factor(seq_len(n_trial)),
                 y            = id_effects[i] +
                                  rnorm(n_trial, 0, sqrt(V_within)))
    }))
    r <- tryCatch(
      rpt(y ~ trial_number + (1|ID), grname = "ID", data = sim_d,
          datatype = "Gaussian", nboot = 0, npermut = 0, verbose = FALSE),
      error = function(e) NULL
    )
    if (!is.null(r) && r$P$LRT_P[1] < 0.05) sig_count <- sig_count + 1
  }
  sig_count / nsim
}

n_grid <- c(10, 15, 20, 30, 40, 50, 60)
R_grid <- c(0.10, 0.20, 0.33, 0.42, 0.50)
power_single <- expand.grid(n_ind = n_grid, R_true = R_grid) %>%
  rowwise() %>%
  mutate(power = simulate_power_single(R_true, n_ind, nsim = 500)) %>%
  ungroup()

write.csv(power_single, "Table_power_single_group.csv", row.names = FALSE)

min_n_80 <- power_single %>%
  group_by(R_true) %>%
  summarise(min_n_80 = {
    idx <- which(power >= 0.80)
    if (length(idx) == 0) NA_integer_ else n_grid[min(idx)]
  }, .groups = "drop")
power_single %>%
  mutate(R_label = paste0("R = ", R_true)) %>%
  select(n_ind, R_label, power) %>%
  pivot_wider(names_from = R_label, values_from = power) %>%
  rename(`n per group` = n_ind) %>%
  mutate(across(where(is.numeric), ~round(.x, 2))) %>%
  kbl(caption = "Power to detect R > 0 (single group, 3 trials, alpha = 0.05, nsim = 500)") %>%
  kable_styling(bootstrap_options = c("striped","hover"))

min_n_80 %>%
  rename(`True R` = R_true, `Min n for 80% power` = min_n_80) %>%
  kbl(caption = "Minimum n per group for 80% power (single-group test)") %>%
  kable_styling(bootstrap_options = c("striped","hover"), full_width = FALSE)

# ---- 12B — Two-group comparison power: detect difference between groups ----
simulate_power_comparison <- function(R_A, R_B, n_per_group,
                                      n_trial = 3, nsim = 500, seed = 99) {
  set.seed(seed)
  sig_count <- 0

  sim_group <- function(R_true, n_ind) {
    V_a <- R_true; V_w <- 1 - R_true
    id_fx <- rnorm(n_ind, 0, sqrt(V_a))
    do.call(rbind, lapply(seq_len(n_ind), function(i) {
      data.frame(ID           = factor(i),
                 trial_number = factor(seq_len(n_trial)),
                 y            = id_fx[i] + rnorm(n_trial, 0, sqrt(V_w)))
    }))
  }

  clamp <- function(r) max(min(r, 0.9999), 0.0001)
  fz    <- function(r) 0.5 * log((1 + clamp(r)) / (1 - clamp(r)))
  # Large-sample SE of the ICC for a balanced one-way random-effects design
  # (a individuals x k measures): Var(R) = 2(1-R)^2 (1+(k-1)R)^2 / [k(k-1)(a-1)]
  # (Searle, Casella & McCulloch 1992). Replaces the inner bootstrap so the
  # simulation runs in seconds rather than hours; the trial fixed effect has
  # negligible effect on this approximation.
  icc_se <- function(R, a, k) {
    R <- clamp(R)
    sqrt((2 * (1 - R)^2 * (1 + (k - 1) * R)^2) / (k * (k - 1) * (a - 1)))
  }

  for (s in seq_len(nsim)) {
    dA <- sim_group(R_A, n_per_group)
    dB <- sim_group(R_B, n_per_group)

    rA <- tryCatch(
      rpt(y ~ trial_number + (1|ID), grname = "ID", data = dA,
          datatype = "Gaussian", nboot = 0, npermut = 0, verbose = FALSE),
      error = function(e) NULL)
    rB <- tryCatch(
      rpt(y ~ trial_number + (1|ID), grname = "ID", data = dB,
          datatype = "Gaussian", nboot = 0, npermut = 0, verbose = FALSE),
      error = function(e) NULL)

    if (is.null(rA) || is.null(rB)) next

    R_est_A <- rA$R$ID[1]; R_est_B <- rB$R$ID[1]
    se_Za   <- icc_se(R_est_A, n_per_group, n_trial) / (1 - clamp(R_est_A)^2)
    se_Zb   <- icc_se(R_est_B, n_per_group, n_trial) / (1 - clamp(R_est_B)^2)
    se_diff <- sqrt(se_Za^2 + se_Zb^2)
    if (!is.finite(se_diff) || se_diff == 0) next
    z_stat  <- (fz(R_est_A) - fz(R_est_B)) / se_diff
    p_val   <- 2 * pnorm(abs(z_stat), lower.tail = FALSE)
    if (p_val < 0.05) sig_count <- sig_count + 1
  }
  sig_count / nsim
}

comparisons <- list(
  # Near-wall: observed per-group R = 0.37 (2) / 0.06 (5) / 0.19 (15 min)
  list(label = "Near-wall: 2 vs 5 min (R=0.37 vs 0.06)",  R_A = 0.37, R_B = 0.06),
  list(label = "Near-wall: 2 vs 15 min (R=0.37 vs 0.19)", R_A = 0.37, R_B = 0.19),
  # Open-zone: observed per-group R = 0.42 (2) / 0.33 (5) / 0.01 (15 min)
  list(label = "Open-zone: 2 vs 5 min (R=0.42 vs 0.33)",  R_A = 0.42, R_B = 0.33),
  list(label = "Open-zone: 2 vs 15 min (R=0.42 vs 0.01)", R_A = 0.42, R_B = 0.01)
)
n_comp_grid <- c(10, 15, 20, 30, 40, 50, 60, 80, 100)
power_comp_list <- list()
for (comp in comparisons) {
  pwr <- sapply(n_comp_grid, function(n) {
    simulate_power_comparison(comp$R_A, comp$R_B, n, nsim = 500)
  })
  power_comp_list[[comp$label]] <- data.frame(
    Comparison  = comp$label,
    n_per_group = n_comp_grid,
    power       = round(pwr, 3)
  )
}
power_comp <- do.call(rbind, power_comp_list)
rownames(power_comp) <- NULL
write.csv(power_comp, "Table_power_group_comparison.csv", row.names = FALSE)

min_n_comp <- power_comp %>%
  group_by(Comparison) %>%
  summarise(min_n_80 = {
    idx <- which(power >= 0.80)
    if (length(idx) == 0) paste0(">", max(n_comp_grid))
    else as.character(n_comp_grid[min(idx)])
  }, .groups = "drop")
write.csv(min_n_comp, "Table_power_min_n.csv", row.names = FALSE)
power_comp %>%
  pivot_wider(names_from = Comparison, values_from = power) %>%
  rename(`n per group` = n_per_group) %>%
  kbl(caption = "Power to detect difference between groups (Fisher's Z, alpha = 0.05, nsim = 500)",
      digits = 3) %>%
  kable_styling(bootstrap_options = c("striped","hover"))

min_n_comp %>%
  kbl(caption = "Minimum n per group for 80% power (two-group comparison)") %>%
  kable_styling(bootstrap_options = c("striped","hover"), full_width = FALSE)

# ---- Power curves ----
p_power_single <- ggplot(power_single,
    aes(x = n_ind, y = power,
        colour = factor(R_true), group = factor(R_true))) +
  geom_hline(yintercept = 0.80, linetype = "dashed", colour = "grey40") +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  scale_colour_viridis_d(name = "True R", option = "C", end = 0.85) +
  scale_x_continuous(breaks = n_grid) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
                     labels = scales::percent_format(accuracy = 1)) +
  geom_vline(xintercept = 10, linetype = "dotted",
             colour = "red", linewidth = 0.6) +
  annotate("text", x = 11, y = 0.05, label = "Current\nn = 10",
           hjust = 0, size = 2.8, colour = "red") +
  annotate("text", x = max(n_grid) * 0.95, y = 0.82,
           label = "80% power", hjust = 1, size = 3, colour = "grey30") +
  labs(x = "Number of individuals per group",
       y = "Power (alpha = 0.05)",
       # The panel letter comes from plot_annotation(tag_levels = "A") below;
       # keeping it in the title as well printed it twice.
       title = "Single-group: power to detect R > 0 (3 trials)") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "right")

p_power_comp <- ggplot(power_comp,
    aes(x = n_per_group, y = power,
        colour = Comparison, group = Comparison)) +
  geom_hline(yintercept = 0.80, linetype = "dashed", colour = "grey40") +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  # breaks = legend order only (near-wall first, as in the other figures);
  # the colour assigned to each comparison is unchanged.
  scale_colour_viridis_d(option = "C", end = 0.85, name = NULL,
                         breaks = sapply(comparisons, `[[`, "label")) +
  scale_x_continuous(breaks = n_comp_grid) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
                     labels = scales::percent_format(accuracy = 1)) +
  geom_vline(xintercept = 10, linetype = "dotted",
             colour = "red", linewidth = 0.6) +
  # Label placed at the top: at the bottom it overlapped the low-power curves.
  annotate("text", x = 11, y = 0.95, label = "Current\nn = 10",
           hjust = 0, size = 2.8, colour = "red") +
  annotate("text", x = max(n_comp_grid) * 0.95, y = 0.82,
           label = "80% power", hjust = 1, size = 3, colour = "grey30") +
  # Two-row legend: on one row the four labels were wider than the figure
  # and the first key was cut off.
  guides(colour = guide_legend(nrow = 2)) +
  labs(x = "Number of individuals per group",
       y = "Power (alpha = 0.05)",
       title = "Two-group: power to detect difference in R (Fisher's Z)") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 8))

p_power_single / p_power_comp

fig_power <- p_power_single / p_power_comp +
  plot_annotation(tag_levels = "A")

save_fig("Figure_S2.png", fig_power, width = 24, height = 20)

# ==== Cox proportional hazards assumption check ====
ph_test <- function(dd, fixed_str, label) {
  dd  <- droplevels(dd)
  frm <- as.formula(paste0(
    "Surv(surv_time, surv_event) ~ ", fixed_str,
    " + frailty(ID, distribution = 'gamma')"
  ))
  m <- coxph(frm, data = dd)
  z <- cox.zph(m)
  global_row <- z$table[nrow(z$table), ]
  list(
    label  = label,
    global = data.frame(
      Group   = label,
      chi2    = round(global_row["chisq"], 3),
      df      = global_row["df"],
      p_value = round(global_row["p"], 4)
    ),
    zph = z
  )
}
ph_2min    <- ph_test(d2,  "trial_number + SLz", "2 min")
ph_5min    <- ph_test(d5,  "trial_number + SLz", "5 min")
ph_15min   <- ph_test(d15, "trial_number + SLz", "15 min")

ph_tab <- do.call(rbind, lapply(
  list(ph_2min, ph_5min, ph_15min),
  function(x) x$global
))
ph_tab %>%
  kbl(caption = "Global Schoenfeld test for proportional hazards assumption",
      digits = 3) %>%
  kable_styling(bootstrap_options = c("striped","hover"), full_width = FALSE)

# ==== Sex as covariate: sensitivity of R ====

# ---- LRT for sex effect ----
# Sex is coded 1/2 (relabelled Male/Female), so the coefficient is "SexFemale";
# grab the Sex row by pattern to stay robust to the label. Run for both measures.
sex_lrt_row <- function(yvar, label) {
  m1 <- lmer(as.formula(paste0(yvar, " ~ trial_number + acclimation_time + SLz + Sex + (1|ID)")),
             REML = FALSE, data = d)
  m0 <- lmer(as.formula(paste0(yvar, " ~ trial_number + acclimation_time + SLz + (1|ID)")),
             REML = FALSE, data = d)
  a  <- anova(m0, m1)
  cf <- summary(m1)$coefficients
  sr <- grep("^Sex", rownames(cf), value = TRUE)[1]
  data.frame(Response = label, Term = sr,
             Chi2 = round(a$Chisq[2], 3), df = a$Df[2],
             P_value  = round(a$`Pr(>Chisq)`[2], 4),
             Sex_coef = round(cf[sr, "Estimate"], 3),
             Sex_SE   = round(cf[sr, "Std. Error"], 3), row.names = NULL)
}

# Boldness is censored -> Cox mixed model (coxme); same LRT idea, coef via fixef().
sex_lrt_cox <- function(label) {
  m1 <- coxme(Surv(surv_time, surv_event) ~ trial_number + acclimation_time + SLz + Sex + (1|ID), data = d)
  m0 <- coxme(Surv(surv_time, surv_event) ~ trial_number + acclimation_time + SLz + (1|ID), data = d)
  a  <- anova(m0, m1); pcol <- grep("^P", colnames(a), value = TRUE)[1]
  cf <- fixef(m1); vc <- sqrt(diag(vcov(m1)))
  sr <- grep("^Sex", names(cf), value = TRUE)[1]
  data.frame(Response = label, Term = sr,
             Chi2 = round(a$Chisq[2], 3), df = a$Df[2],
             P_value  = round(a[[pcol]][2], 4),
             Sex_coef = round(unname(cf[sr]), 3),
             Sex_SE   = round(unname(vc[sr]), 3), row.names = NULL)
}

rbind(
  sex_lrt_cox("Boldness (latency to exit shelter)"),
  sex_lrt_row("near_logit",  "Near-wall use"),
  sex_lrt_row("space_logit", "Open-zone use")
) %>%
  kbl(caption = "LRT for the sex effect (overall models): boldness (coxme) and the two space-use measures (lmer).") %>%
  kable_styling(bootstrap_options = c("striped","hover"), full_width = FALSE)

# ---- R with sex as fixed effect ----
rptR_sex <- function(dd, fixed_str, label, yvar = "space_logit", seed = 42) {
  dd <- droplevels(dd); set.seed(seed)
  frm <- as.formula(paste0(yvar, " ~ ", fixed_str, " + (1|ID)"))
  r   <- rpt(frm, grname = "ID", data = dd,
             datatype = "Gaussian", nboot = 1000, npermut = 1000,
             verbose = FALSE)
  data.frame(
    Model = "With sex", Group = label,
    R     = round(r$R$ID[1], 3),
    CI_lo = round(quantile(r$R_boot$ID, 0.025, na.rm = TRUE), 3),
    CI_hi = round(quantile(r$R_boot$ID, 0.975, na.rm = TRUE), 3),
    P_LRT = round(r$P$LRT_P[1], 4),
    P_perm= round(r$P$P_permut[1], 4),
    row.names = NULL
  )
}

mk_sex <- function(yvar, trait) {
  t <- rbind(
    rptR_sex(d2,  "trial_number + SLz + Sex", "2 min",  yvar = yvar),
    rptR_sex(d5,  "trial_number + SLz + Sex", "5 min",  yvar = yvar),
    rptR_sex(d15, "trial_number + SLz + Sex", "15 min", yvar = yvar))
  t$Trait <- trait
  t %>% select(Trait, Model, Group, R, CI_lo, CI_hi, P_LRT, P_perm)
}
without <- function(tab, trait) tab %>%
  mutate(Model = "Without sex", Trait = trait) %>%
  select(Trait, Model, Group, R, CI_lo, CI_hi, P_LRT, P_perm)

# Boldness (censored): ICC with vs without Sex, per group, via coxme.
bold_with <- rbind(
  icc_coxme(subset(d, acclimation_time == "2"),  "trial_number + SLz + Sex", "2 min"),
  icc_coxme(subset(d, acclimation_time == "5"),  "trial_number + SLz + Sex", "5 min"),
  icc_coxme(subset(d, acclimation_time == "15"), "trial_number + SLz + Sex", "15 min"))
bold_cmp <- function(tab, model) tab %>%
  transmute(Trait = "Latency to exit shelter", Model = model, Group,
            R = ICC, CI_lo, CI_hi, P_LRT = p_LRT, P_perm = NA_real_)

sex_comparison <- rbind(
  bold_cmp(bold_tab,  "Without sex"),
  bold_cmp(bold_with, "With sex"),
  without(near_tab,  "Near-wall"),
  mk_sex("near_logit", "Near-wall"),
  without(space_tab, "Open-zone"),
  mk_sex("space_logit", "Open-zone")
) %>% arrange(Trait, Group, Model)

write.csv(sex_comparison, "Table_sex_covariate_comparison.csv", row.names = FALSE)
sex_comparison %>%
  mutate(across(where(is.numeric), ~round(.x, 3))) %>%
  kbl(caption = "R with and without sex as a fixed effect: boldness (coxme ICC, latent scale) and the two space-use measures (lmer). Scales not comparable across traits.") %>%
  kable_styling(bootstrap_options = c("striped","hover")) %>%
  collapse_rows(columns = 1, valign = "top")
ggplot(
  sex_comparison %>%
    # The figure shows the three acclimation groups, as Table S3 does; the
    # pooled "Overall" rows are dropped here (they previously fell outside the
    # factor levels and were drawn as an unlabelled "NA" category).
    filter(Group %in% c("2 min", "5 min", "15 min")) %>%
    mutate(Group = factor(Group, levels = c("2 min","5 min","15 min")),
           Trait = factor(ifelse(Trait == "Latency to exit shelter", Trait,
                                 paste(Trait, "use")), levels = trait_levels),
           Model = factor(Model, levels = c("Without sex", "With sex"))),
  aes(x = Group, y = R, colour = Model, group = Model)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_errorbar(aes(ymin = CI_lo, ymax = CI_hi), width = 0.12, linewidth = 0.5,
                position = position_dodge(width = 0.3)) +
  geom_point(size = 2.5, position = position_dodge(width = 0.3)) +
  geom_line(position = position_dodge(width = 0.3), linewidth = 0.5) +
  facet_wrap(~ Trait) +
  scale_colour_manual(values = sex_cols, name = NULL) +
  scale_y_continuous(limits = c(-0.1, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = "Acclimation time", y = "Repeatability (R / ICC)") +
  theme_fig

save_fig("Figure_S4.png", width = 22, height = 10)

# ==== Adjusted vs unadjusted repeatability ====
get_both_R <- function(dd, fixed_str, label, yvar = "space_logit", seed = 42) {
  dd <- droplevels(dd); set.seed(seed)
  # adjusted R: with fixed effects (their variance excluded from the denominator)
  frm_adj <- as.formula(paste0(yvar, " ~ ", fixed_str, " + (1|ID)"))
  r       <- rpt(frm_adj, grname = "ID", data = dd,
                 datatype = "Gaussian", nboot = 1000, npermut = 0, verbose = FALSE)
  # unadjusted R: rptR has no $R_org, so refit WITHOUT fixed effects (their
  # variance then stays in the denominator) -> the unadjusted repeatability.
  set.seed(seed)
  r_un    <- rpt(as.formula(paste0(yvar, " ~ (1|ID)")), grname = "ID", data = dd,
                 datatype = "Gaussian", nboot = 0, npermut = 0, verbose = FALSE)
  R_adj   <- round(r$R$ID[1], 3)
  R_unadj <- round(r_un$R$ID[1], 3)
  CI_lo   <- round(quantile(r$R_boot$ID, 0.025, na.rm = TRUE), 3)
  CI_hi   <- round(quantile(r$R_boot$ID, 0.975, na.rm = TRUE), 3)
  data.frame(Group = label,
             R_adjusted   = R_adj,
             R_unadjusted = R_unadj,
             Delta        = round(R_adj - R_unadj, 3),
             CI_lo_adj    = CI_lo,
             CI_hi_adj    = CI_hi,
             row.names = NULL)
}

mk_adj <- function(yvar, trait) {
  t <- rbind(
    get_both_R(d2,  "trial_number + SLz", "2 min",  yvar = yvar),
    get_both_R(d5,  "trial_number + SLz", "5 min",  yvar = yvar),
    get_both_R(d15, "trial_number + SLz", "15 min", yvar = yvar))
  cbind(Trait = trait, t)
}
adj_comp <- rbind(
  mk_adj("near_logit",  "Near-wall"),
  mk_adj("space_logit", "Open-zone")
)
write.csv(adj_comp, "Table_adjusted_vs_unadjusted_R.csv", row.names = FALSE)
adj_comp %>%
  kbl(caption = "Adjusted vs unadjusted R, both space-use measures. Delta = R_adjusted − R_unadjusted.",
      digits = 3) %>%
  kable_styling(bootstrap_options = c("striped","hover"), full_width = FALSE) %>%
  collapse_rows(columns = 1, valign = "top")