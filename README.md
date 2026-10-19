# Sensitivity Analysis

Shiny app for **mixed-model power / sensitivity curves** (clustered or twin data) using [`simr`](https://cran.r-project.org/package=simr) and `lme4`.

**Live app:** [https://abhiramn.shinyapps.io/sensitivity-analysis/](https://abhiramn.shinyapps.io/sensitivity-analysis/)

Lab members upload their own CSV, choose the cluster ID, focal predictor, outcomes, covariates, and which variables to standardize or treat as factors, then run simulation-based power curves across effect sizes.

## Run locally

```r
# from this folder
shiny::runApp("app.R")
```

Install packages once (do **not** rely on runtime installs — they break shinyapps.io):

```r
install.packages(c(
  "shiny", "bslib", "DT", "tidyverse",
  "simr", "lme4", "lmerTest", "pwr", "scales"
))
```

## Data

Upload a CSV with:

- A **cluster / pair ID** (e.g. `pairid`)
- A **focal predictor** (IV)
- One or more **outcomes** (DVs; each gets its own simulation)
- Optional **covariates**
- Optional **ID columns** to exclude from the model (e.g. person IDs)

Character columns can be read as factors. Continuous variables selected for standardization are modeled as `*_std` columns (same convention as `sensitivity_Analysis_V02.R`).

`test.csv` is optional synthetic example data for trying the app. Use your own study data for real analyses.

## Deploy to shinyapps.io

**RStudio (recommended):** open `app.R` → **Publish** → shinyapps.io → select `app.R` → publish.

**Or in R** (after linking your shinyapps account once):

```r
rsconnect::deployApp(
  appFiles = "app.R",
  appName  = "sensitivity-analysis",
  forceUpdate = TRUE,
  launch.browser = TRUE
)
```

## Files

| File | Role |
|------|------|
| `app.R` | Shiny app |
| `sensitivity_Analysis_V02.R` | Original scripted workflow |
| `test.csv` | Optional synthetic example data |
