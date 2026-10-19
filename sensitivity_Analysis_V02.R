############################
### SENSITIVITY ANALYSES ###
############################


# ==============================================================================
# PACKAGE CONFIGURATION ----
# ==============================================================================

## --- Check if Packages are Installed ----------
packages <- c("tidyverse", "simr", "lme4", "lmerTest")

installed <- packages %in% rownames(installed.packages())

if (any(!installed)) {
  install.packages(packages[!installed])
}

lapply(packages, library, character.only = TRUE)


## --- Load Packages -----------------------------
library(tidyverse)
library(simr)
library(lme4)
library(lmerTest)
library(pwr)

## --- Package Citation (for reporting) --------------
citation("simr")


# ==============================================================================
# DATA CONFIGURATION ----
# ==============================================================================

## --- File path -----------------------------------
DATA_PATH <- "~/Desktop/BiG_WSTR/Data/ciara_project/test.csv" # Path to your CSV file

## --- L2 ID variable ----------------
PAIR_ID <- "pairid"  # Column identifying L2 units (e.g., twin pairs). Must match a column name in the data

## --- Predictor of interest (H1) --------------------
# This is a focal predictor whose power you want to evaluate.
# **Must match a column name in your data AND append _std to end if it will eventually be standardized**
PREDICTOR <- "mv_std"  

## --- Outcomes (H1) ----------------------------------
# Each entry defines one outcome model.
# - outcome_var : column name for the dependent variable WITH _std appended to end if it will eventually be standardized
# - label : label for the variable that will be used in the final figure
OUTCOMES <- list(
  list(outcome_var = "connection_std", label = "Connectedness"),
  list(outcome_var = "DunedinPACE_std", label = "DunedinPACE"),
  list(outcome_var = "GrimAgeAA_std", label = "GrimAge2 AA"))

## --- Covariates and Focal Predictors (shared across all H1 models) --------
# List the raw column names PLUS _std for any variables that will eventually be standardized
COVARIATES <- c("mv_std", "ave_age_std", "sexF", "zygF",
                "w_nhp", "smoker", "relstat")

## --- Variables to standardize before modeling ------------------
# List the raw column names; standardized versions will be created with _std suffix.
STANDARDIZE <- c("mv", "connection", "ave_age", "DunedinPACE", "GrimAgeAA")

## --- H1 power simulation settings --------------------------------------------
EFFECT_SIZES <- c(0.05, 0.10, 0.15, 0.20) # Standardized B values to evaluate
NSIM <- 10  # Simulations per effect size 
  # Use 10 simulations to check code is working. Use at least 1,000 simulations for final reporting.
SEED <- 1995  # For reproducibility purposes


# ==============================================================================
# MAIN FUNCTIONS ----
# ==============================================================================

# Standardize named columns and append _std versions to data
add_std_cols <- function(df, vars) {
  for (v in vars) {
    new_col <- paste0(v, "_std")
    df[[new_col]] <- as.numeric(scale(df[[v]]))} 
  df
}

# Filter data to complete pairs for a given set of variables (required for simr)
filter_complete_pairs <- function(df, vars, pair_col, label = "") {
  df$..complete.. <- complete.cases(df[, vars])
  keep <- df %>%
    group_by(across(all_of(pair_col))) %>%
    summarise(ok = all(..complete..), .groups = "drop") %>%
    filter(ok) %>%
    pull(all_of(pair_col))
  out <- df %>% filter(.data[[pair_col]] %in% keep) %>% dplyr::select(-..complete..)
  message(sprintf("[%s] Complete pairs: %d", label, length(unique(out[[pair_col]]))))
  out
}

# Fit mixed-effects model for one outcome
fit_lmer_model <- function(df, outcome_var, predictor, covariates, pair_col) {
  rhs <- paste(c(predictor, setdiff(covariates, predictor), paste0("(1 | ", pair_col, ")")),
                    collapse = " + ")
  formula <- as.formula(paste(outcome_var, "~", rhs))
  eval(bquote(lmer(.(formula), data = df, REML = FALSE,
                   control = lmerControl(optimizer = "bobyqa"))))
}

# Power simulation across vector of effect sizes
power_curve_simr <- function(model,
                             effect_sizes,
                             nsim = 200,
                             param = PREDICTOR,
                             test_var = PREDICTOR) {
  map_dfr(effect_sizes, function(b) {
    m <- model
    fixef(m)[param] <- b
    
    sim <- powerSim(m, test = fixed(test_var),
                    nsim = nsim, progress = TRUE)
    
    p_hat <- mean(sim$pval < sim$alpha, na.rm = TRUE)
    ci <- binom.test(x = sum(sim$pval < sim$alpha, na.rm = TRUE),
                     n = length(sim$pval),
                     conf.level = 0.95)$conf.int
    
    data.frame(effect_size = b, power = p_hat,
               lower = ci[1], upper = ci[2],
               n = sim$nrow, nsim = nsim)
  })
}


# ==============================================================================
# DATA PREPARATION ----
# ==============================================================================

## --- Read in Data File ------------------------------------------
df <- read.csv(DATA_PATH, stringsAsFactors = TRUE)

## --- Standardize variables specified earlier --------------------
df <- add_std_cols(df, STANDARDIZE) 

# ==============================================================================
# H1: FIT MODEL(S) ----
# ==============================================================================

## --- Filter to complete pairs, fit model(s), store results ------------------
h1_models <- map(OUTCOMES, function(spec) {
  
  all_vars <- c(spec$outcome_var, COVARIATES)
  
  df_complete <- filter_complete_pairs(df = df, vars = all_vars,
                                       pair_col = PAIR_ID, label = spec$label)
  
  model <- fit_lmer_model(df = df_complete, outcome_var = spec$outcome_var,
                          predictor = PREDICTOR, covariates = COVARIATES,
                          pair_col = PAIR_ID)
  
  list(label = spec$label, model = model, data = df_complete)
})


## --- Name Lists by Outcomes -------------------------------------

names(h1_models) <- map_chr(OUTCOMES, "label")

# ==============================================================================
# H1: POWER SIMULATION ----
# ==============================================================================

set.seed(SEED)
power_results <- map_dfr(h1_models, function(spec) {
  power_curve_simr(model = spec$model, effect_sizes = EFFECT_SIZES, 
                   nsim = NSIM) %>% 
    mutate(outcome = spec$label)
})


outcome_levels <- map_chr(OUTCOMES, "label")
power_results <- power_results %>%
  mutate(outcome = factor(outcome, levels = outcome_levels))


# ==============================================================================
# H1: PLOT POWER CURVE ----
# ==============================================================================

fig_h1 <- ggplot(power_results, aes(x = effect_size, y = power)) +
  geom_hline(yintercept = 0.80, linetype = "dashed",
             color = "gray50", linewidth = 0.5) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha = 0.2, fill = "steelblue") +
  geom_line(color  = "steelblue", linewidth = 0.8) +
  geom_point(color = "steelblue", size = 2) +
  facet_wrap(~ outcome, nrow = 1) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.10),
                     minor_breaks = seq(0, 1, by = 0.05), 
                     labels = scales::percent_format(accuracy = 1)) +
  scale_x_continuous(breaks = unique(power_results$effect_size)) +
  labs(x = paste0("Standardized ", PREDICTOR, " Effect Size (\u03b2)"),
       y = "Power") +
  theme_bw() +
  theme(text = element_text(family = "Times New Roman"),
        strip.text = element_text(face = "bold", size = 11), 
        axis.title = element_text(size = 11), axis.text = element_text(size = 9),
        panel.grid.minor = element_line(color = "gray90", linewidth = 0.3))

print(fig_h1)
print(power_results) # corresponding values
