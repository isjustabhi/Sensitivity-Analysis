############################
### Sensitivity Analysis Tool ###
### Lab mixed-model power / sensitivity app ###
############################

# ==============================================================================
# PACKAGES ----
# Do not call install.packages() here — shinyapps.io installs deps at deploy time.
# Locally, install once: install.packages(c("shiny","bslib","DT","tidyverse",
#   "simr","lme4","lmerTest","pwr","scales"))
# ==============================================================================

library(shiny)
library(bslib)
library(DT)
library(tidyverse)
library(simr)
library(lme4)
library(lmerTest)
library(pwr)
library(scales)

# Null-coalesce helper (works on older R / shinyapps images)
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

# ==============================================================================
# ANALYSIS HELPERS (from sensitivity_Analysis_V02.R) ----
# ==============================================================================

add_std_cols <- function(df, vars) {
  for (v in vars) {
    if (!v %in% names(df)) next
    if (!is.numeric(df[[v]])) {
      warning(sprintf("Skipping standardization of non-numeric column: %s", v))
      next
    }
    df[[paste0(v, "_std")]] <- as.numeric(scale(df[[v]]))
  }
  df
}

filter_complete_pairs <- function(df, vars, pair_col, label = "") {
  vars <- intersect(vars, names(df))
  df$..complete.. <- complete.cases(df[, vars, drop = FALSE])
  keep <- df %>%
    group_by(across(all_of(pair_col))) %>%
    summarise(ok = all(..complete..), .groups = "drop") %>%
    filter(ok) %>%
    pull(all_of(pair_col))
  out <- df %>%
    filter(.data[[pair_col]] %in% keep) %>%
    dplyr::select(-..complete..)
  n_pairs <- length(unique(out[[pair_col]]))
  list(data = out, n_pairs = n_pairs, label = label)
}

fit_lmer_model <- function(df, outcome_var, predictor, covariates, pair_col) {
  rhs_terms <- c(predictor, setdiff(covariates, predictor), paste0("(1 | ", pair_col, ")"))
  formula <- as.formula(paste(outcome_var, "~", paste(rhs_terms, collapse = " + ")))
  eval(bquote(lmer(
    .(formula),
    data = df,
    REML = FALSE,
    control = lmerControl(optimizer = "bobyqa")
  )))
}

power_curve_simr <- function(model, effect_sizes, nsim, param, test_var,
                             progress_cb = NULL) {
  map_dfr(seq_along(effect_sizes), function(i) {
    b <- effect_sizes[i]
    if (!is.null(progress_cb)) {
      progress_cb(i, length(effect_sizes), b)
    }
    m <- model
    fixef(m)[param] <- b
    sim <- powerSim(m, test = fixed(test_var), nsim = nsim, progress = FALSE)
    n_sig <- sum(sim$pval < sim$alpha, na.rm = TRUE)
    n_tot <- length(sim$pval)
    p_hat <- mean(sim$pval < sim$alpha, na.rm = TRUE)
    ci <- binom.test(x = n_sig, n = n_tot, conf.level = 0.95)$conf.int
    data.frame(
      effect_size = b,
      power = p_hat,
      lower = ci[1],
      upper = ci[2],
      n = sim$nrow,
      nsim = nsim
    )
  })
}

std_name <- function(v, standardize_vars) {
  if (v %in% standardize_vars) paste0(v, "_std") else v
}

guess_factor_cols <- function(df) {
  names(df)[vapply(df, function(x) {
    is.factor(x) || is.character(x) || is.logical(x)
  }, logical(1))]
}

guess_numeric_cols <- function(df) {
  names(df)[vapply(df, is.numeric, logical(1))]
}

# ==============================================================================
# UI ----
# ==============================================================================

app_theme <- bs_theme(
  version = 5,
  bootswatch = "flatly",
  primary = "#1B4F72",
  secondary = "#5D6D7E",
  font_scale = 0.95
)

ui <- page_navbar(
  title = "Sensitivity Analysis",
  theme = app_theme,
  fillable = FALSE,

  nav_panel(
    "Analysis",
    layout_sidebar(
      sidebar = sidebar(
        width = 360,
        open = "open",

        h5("1. Data"),
        fileInput("csv_file", "Upload CSV", accept = c(".csv", "text/csv")),
        helpText("Upload your own analysis dataset. Variable choices update after upload."),
        checkboxInput("strings_as_factors", "Read character columns as factors", TRUE),
        actionButton(
          "toggle_preview",
          "Show data preview",
          class = "btn-outline-secondary w-100",
          icon = icon("table")
        ),
        hr(),

        h5("2. Clustering & variables"),
        selectInput("pair_id", "Pair / cluster ID (L2)", choices = NULL),
        selectInput("exclude_ids", "Exclude ID columns (not modeled)",
                    choices = NULL, multiple = TRUE),
        selectInput("predictor", "Focal predictor (IV)", choices = NULL),
        selectInput("outcomes", "Focal outcomes (DVs)",
                    choices = NULL, multiple = TRUE),
        selectInput("covariates", "Covariates (plus focal IV is always included)",
                    choices = NULL, multiple = TRUE),
        selectInput(
          "factor_vars",
          "Treat as factors (categorical)",
          choices = NULL,
          multiple = TRUE
        ),
        selectInput(
          "standardize",
          "Standardize (creates _std columns for modeling)",
          choices = NULL,
          multiple = TRUE
        ),
        helpText(
          "Factors stay as factors in the model. Only numeric columns should be standardized."
        ),
        hr(),

        h5("3. Simulation settings"),
        textInput(
          "effect_sizes",
          "Effect sizes (comma-separated β)",
          value = "0.05, 0.10, 0.15, 0.20"
        ),
        numericInput("nsim", "Simulations per effect size (nsim)",
                     value = 10, min = 1, max = 5000, step = 1),
        helpText("Use ~10 to test the app; use ≥1000 for final reporting."),
        numericInput("seed", "Random seed", value = 1995, min = 1, step = 1),
        hr(),

        actionButton(
          "run_btn",
          "Run sensitivity analysis",
          class = "btn-primary w-100",
          icon = icon("play")
        ),
        br(), br(),
        downloadButton("download_results", "Download results CSV", class = "w-100"),
        downloadButton("download_plot", "Download plot PNG", class = "w-100")
      ),

      layout_columns(
        col_widths = 12,
        conditionalPanel(
          condition = "input.toggle_preview % 2 == 1",
          card(
            card_header("Data preview"),
            card_body(
              verbatimTextOutput("data_summary"),
              DTOutput("data_preview")
            )
          )
        ),
        card(
          full_screen = TRUE,
          card_header("Sensitivity curves"),
          card_body(
            div(
              style = "width:100%; min-height:520px;",
              plotOutput("power_plot", height = "520px", width = "100%")
            ),
            verbatimTextOutput("run_log")
          )
        ),
        card(
          card_header("Results table"),
          card_body(DTOutput("results_table"))
        )
      )
    )
  ),

  nav_panel(
    "Help",
    card(
      card_header("How to use this tool"),
      tags$ol(
        tags$li("Upload your twin/cluster CSV (any study dataset with a cluster ID)."),
        tags$li("Choose the cluster ID (e.g., pair / family ID). Exclude person IDs from modeling."),
        tags$li("Set the focal IV, one or more DVs (each gets its own simulation), and covariates."),
        tags$li(
          "Mark categorical columns as factors. Binary codes like sexF / zygF can stay numeric ",
          "or be converted to factors — factors enter the model as categorical predictors."
        ),
        tags$li(
          "Select continuous variables to standardize. Modeling then uses the ",
          code("_std"), " versions (matching the sensitivity analysis script convention)."
        ),
        tags$li(
          "Start with a small ", code("nsim"), " (e.g., 10). Increase to ≥1000 for manuscript figures."
        )
      ),
      tags$p(
        strong("Example setup:"),
        " IV = mv; DVs = connection, GrimAgeAA, DunedinPACE; ",
        "covariates = all other non-ID columns; standardize continuous focal vars."
      )
    )
  )
)

# ==============================================================================
# SERVER ----
# ==============================================================================

server <- function(input, output, session) {

  results_rv <- reactiveVal(NULL)
  plot_rv <- reactiveVal(NULL)
  log_rv <- reactiveVal("")

  append_log <- function(...) {
    msg <- paste0(...)
    log_rv(paste0(log_rv(), msg, "\n"))
  }

  raw_data <- reactive({
    validate(need(input$csv_file, "Upload a CSV file to begin."))
    read.csv(input$csv_file$datapath, stringsAsFactors = input$strings_as_factors)
  })

  # Keep selectors in sync when data changes
  observeEvent(raw_data(), {
    df <- raw_data()
    cols <- names(df)
    num_cols <- guess_numeric_cols(df)
    fac_guess <- guess_factor_cols(df)

    # Sensible guesses when common column names are present; otherwise leave open
    default_pair <- if ("pairid" %in% cols) {
      "pairid"
    } else if ("pair_id" %in% cols) {
      "pair_id"
    } else {
      cols[1]
    }

    default_pred <- if ("mv" %in% cols) {
      "mv"
    } else {
      setdiff(num_cols, default_pair)[1]
    }

    default_out <- intersect(c("connection", "GrimAgeAA", "DunedinPACE"), cols)

    id_guess <- intersect(
      c("studyid", "pair_id", "pairid", "fid", "id", "twinid", "person_id"),
      cols
    )
    # Covariates = everything except IDs and outcomes (IV stays included)
    default_cov <- setdiff(cols, c(default_pair, default_out, id_guess, "DNAmGrimAge2BasedOnRealAge"))
    default_std <- intersect(
      unique(c(default_pred, default_out, "age", "ave_age", "bmi", "mv", "connection",
               "DunedinPACE", "GrimAgeAA", "highested", "emp_hhinc")),
      num_cols
    )

    updateSelectInput(session, "pair_id", choices = cols, selected = default_pair)
    updateSelectInput(session, "exclude_ids", choices = cols, selected = setdiff(id_guess, default_pair))
    updateSelectInput(session, "predictor", choices = cols, selected = default_pred)
    updateSelectInput(session, "outcomes", choices = cols, selected = default_out)
    updateSelectInput(session, "covariates", choices = cols, selected = default_cov)
    updateSelectInput(session, "factor_vars", choices = cols, selected = fac_guess)
    updateSelectInput(session, "standardize", choices = num_cols, selected = default_std)
  }, ignoreNULL = TRUE)

  prepared_data <- reactive({
    df <- raw_data()
    facs <- intersect(input$factor_vars %||% character(0), names(df))
    for (v in facs) {
      df[[v]] <- as.factor(df[[v]])
    }
    df
  })

  observeEvent(input$toggle_preview, {
    showing <- (input$toggle_preview %% 2 == 1)
    updateActionButton(
      session,
      "toggle_preview",
      label = if (showing) "Hide data preview" else "Show data preview",
      icon = icon(if (showing) "eye-slash" else "table")
    )
  }, ignoreInit = TRUE)

  output$data_summary <- renderPrint({
    req(input$toggle_preview %% 2 == 1)
    df <- prepared_data()
    facs <- names(df)[vapply(df, is.factor, logical(1))]
    cat(sprintf("Rows: %d | Columns: %d\n", nrow(df), ncol(df)))
    cat(sprintf("Factor columns: %s\n",
                if (length(facs)) paste(facs, collapse = ", ") else "(none)"))
    str(df)
  })

  output$data_preview <- renderDT({
    req(input$toggle_preview %% 2 == 1)
    datatable(
      prepared_data(),
      options = list(pageLength = 8, scrollX = TRUE),
      rownames = FALSE
    )
  })

  output$run_log <- renderText(log_rv())

  observeEvent(input$run_btn, {
    log_rv("")
    results_rv(NULL)
    plot_rv(NULL)

    df0 <- prepared_data()
    pair_col <- input$pair_id
    predictor_raw <- input$predictor
    outcomes_raw <- input$outcomes
    cov_raw <- input$covariates
    std_vars <- input$standardize %||% character(0)
    exclude <- input$exclude_ids %||% character(0)

    if (is.null(pair_col) || !nzchar(pair_col)) {
      showNotification("Select a pair / cluster ID.", type = "error")
      return()
    }
    if (is.null(predictor_raw) || !nzchar(predictor_raw)) {
      showNotification("Select a focal predictor.", type = "error")
      return()
    }
    if (is.null(outcomes_raw) || length(outcomes_raw) < 1) {
      showNotification("Select at least one outcome.", type = "error")
      return()
    }

    effect_sizes <- suppressWarnings(
      as.numeric(trimws(unlist(strsplit(input$effect_sizes, ","))))
    )
    if (!length(effect_sizes) || any(!is.finite(effect_sizes))) {
      showNotification("Enter valid numeric effect sizes.", type = "error")
      return()
    }
    if (is.null(input$nsim) || input$nsim < 1) {
      showNotification("nsim must be ≥ 1.", type = "error")
      return()
    }

    # Drop excluded IDs from covariate set; ensure IV is in covariates
    cov_raw <- setdiff(unique(c(predictor_raw, cov_raw)), c(pair_col, exclude, outcomes_raw))
    std_vars <- setdiff(std_vars, c(pair_col, exclude))
    std_vars <- intersect(std_vars, names(df0))
    # Do not standardize factors
    std_vars <- std_vars[vapply(std_vars, function(v) is.numeric(df0[[v]]), logical(1))]

    withProgress(message = "Running sensitivity analysis", value = 0, {
      append_log("Standardizing: ", paste(std_vars, collapse = ", "))
      df <- add_std_cols(df0, std_vars)

      pred <- std_name(predictor_raw, std_vars)
      covs <- unname(vapply(cov_raw, std_name, character(1), standardize_vars = std_vars))
      outcome_specs <- lapply(outcomes_raw, function(o) {
        list(
          outcome_raw = o,
          outcome_var = std_name(o, std_vars),
          label = o
        )
      })

      append_log("Predictor in models: ", pred)
      append_log("Covariates in models: ", paste(covs, collapse = ", "))

      set.seed(input$seed)
      n_out <- length(outcome_specs)
      n_eff <- length(effect_sizes)
      total_steps <- max(1, n_out * n_eff)
      step_box <- new.env(parent = emptyenv())
      step_box$i <- 0

      power_results <- map_dfr(seq_along(outcome_specs), function(j) {
        spec <- outcome_specs[[j]]
        all_vars <- unique(c(spec$outcome_var, covs))
        filtered <- filter_complete_pairs(
          df = df, vars = all_vars, pair_col = pair_col, label = spec$label
        )
        append_log(sprintf("[%s] Complete pairs: %d", spec$label, filtered$n_pairs))

        if (filtered$n_pairs < 2) {
          append_log(sprintf("[%s] Skipped — not enough complete pairs.", spec$label))
          return(NULL)
        }

        model <- tryCatch(
          fit_lmer_model(
            df = filtered$data,
            outcome_var = spec$outcome_var,
            predictor = pred,
            covariates = covs,
            pair_col = pair_col
          ),
          error = function(e) {
            append_log(sprintf("[%s] Model error: %s", spec$label, e$message))
            NULL
          }
        )
        if (is.null(model)) return(NULL)
        if (!pred %in% names(fixef(model))) {
          append_log(sprintf(
            "[%s] Predictor '%s' not found in fixed effects (got: %s).",
            spec$label, pred, paste(names(fixef(model)), collapse = ", ")
          ))
          return(NULL)
        }

        power_curve_simr(
          model = model,
          effect_sizes = effect_sizes,
          nsim = input$nsim,
          param = pred,
          test_var = pred,
          progress_cb = function(i, n, b) {
            step_box$i <- step_box$i + 1
            setProgress(
              value = step_box$i / total_steps,
              detail = sprintf("%s | β = %s (%d/%d)", spec$label, b, i, n)
            )
          }
        ) %>%
          mutate(outcome = spec$label)
      })

      if (is.null(power_results) || nrow(power_results) == 0) {
        append_log("No results produced. Check variables and complete cases.")
        showNotification("No results produced.", type = "warning")
        return()
      }

      power_results <- power_results %>%
        mutate(outcome = factor(outcome, levels = outcomes_raw))

      n_panels <- length(unique(power_results$outcome))
      facet_nrow <- if (n_panels <= 3) 1L else 2L

      fig <- ggplot(power_results, aes(x = effect_size, y = power)) +
        geom_hline(yintercept = 0.80, linetype = "dashed",
                   color = "gray50", linewidth = 0.5) +
        geom_ribbon(aes(ymin = lower, ymax = upper),
                    alpha = 0.2, fill = "steelblue") +
        geom_line(color = "steelblue", linewidth = 0.8) +
        geom_point(color = "steelblue", size = 2) +
        facet_wrap(~ outcome, nrow = facet_nrow) +
        scale_y_continuous(
          limits = c(0, 1),
          breaks = seq(0, 1, by = 0.10),
          labels = scales::percent_format(accuracy = 1)
        ) +
        scale_x_continuous(breaks = unique(power_results$effect_size)) +
        labs(
          x = paste0("Standardized ", pred, " effect size (\u03b2)"),
          y = "Power",
          title = "Sensitivity analysis"
        ) +
        theme_bw(base_size = 12) +
        theme(
          plot.margin = margin(8, 12, 8, 8),
          strip.text = element_text(face = "bold", size = 11),
          axis.title = element_text(size = 11),
          axis.text = element_text(size = 9),
          panel.grid.minor = element_line(color = "gray90", linewidth = 0.3)
        )

      results_rv(power_results)
      plot_rv(fig)
      append_log("Done.")
      showNotification("Power analysis finished.", type = "message")
    })
  })

  output$power_plot <- renderPlot({
    req(plot_rv())
    # Guard against zero-size plot devices in fillable layouts
    op <- par(no.readonly = TRUE)
    on.exit(par(op), add = TRUE)
    par(mar = c(4, 4, 2, 1))
    print(plot_rv())
  }, width = function() {
    # Prefer session client width when available; fall back safely
    session$clientData$output_power_plot_width %||% 900
  }, height = 520, res = 96)

  output$results_table <- renderDT({
    req(results_rv())
    datatable(
      results_rv() %>%
        mutate(
          power = round(power, 3),
          lower = round(lower, 3),
          upper = round(upper, 3)
        ),
      options = list(pageLength = 12, scrollX = TRUE),
      rownames = FALSE
    )
  })

  output$download_results <- downloadHandler(
    filename = function() paste0("power_results_", Sys.Date(), ".csv"),
    content = function(file) {
      req(results_rv())
      write.csv(results_rv(), file, row.names = FALSE)
    }
  )

  output$download_plot <- downloadHandler(
    filename = function() paste0("power_plot_", Sys.Date(), ".png"),
    content = function(file) {
      req(plot_rv())
      ggsave(file, plot = plot_rv(), width = 11, height = 5.5, dpi = 150)
    }
  )
}

shinyApp(ui, server)
