# Sensitivity Analysis

Shiny app for **mixed-model power / sensitivity curves** (clustered or twin data) using [`simr`](https://cran.r-project.org/package=simr) and `lme4`.

Lab members upload their own CSV, choose the cluster ID, focal predictor, outcomes, covariates, and which variables to standardize or treat as factors, then run simulation-based power curves across effect sizes.

## Run locally

```r
# from this folder
shiny::runApp("app.R")
```

Required packages (installed automatically on first run if missing):  
`shiny`, `bslib`, `DT`, `tidyverse`, `simr`, `lme4`, `lmerTest`, `pwr`, `scales`

## Data

Upload a CSV with:

- A **cluster / pair ID** (e.g. `pairid`)
- A **focal predictor** (IV)
- One or more **outcomes** (DVs; each gets its own simulation)
- Optional **covariates**
- Optional **ID columns** to exclude from the model (e.g. person IDs)

Character columns can be read as factors. Continuous variables selected for standardization are modeled as `*_std` columns (same convention as `sensitivity_Analysis_V02.R`).

An optional synthetic example file, `test.csv`, is included for trying the app. Use your own study data for real analyses.

## Deploy to shinyapps.io

1. Create an account at [shinyapps.io](https://www.shinyapps.io) and copy your token.
2. In R:

```r
rsconnect::setAccountInfo(
  name   = "YOUR_ACCOUNT",
  token  = "YOUR_TOKEN",
  secret = "YOUR_SECRET"
)

setwd("path/to/Sensitivity-Analysis")
rsconnect::deployApp(
  appFiles = c("app.R", "test.csv"),  # omit test.csv if you prefer upload-only on the server
  appName  = "sensitivity-analysis",
  launch.browser = TRUE
)
```

## Files

| File | Role |
|------|------|
| `app.R` | Shiny app |
| `sensitivity_Analysis_V02.R` | Original scripted workflow |
| `test.csv` | Optional synthetic example data |
