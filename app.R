# =============================================================================
# Sensitivity Analysis — Shiny Application
# UI shell adapted from Associational Variable Analysis (AVA) app
# Analysis methodology: sensitivity_Analysis_V02.R (simr / lme4)
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Libraries
# -----------------------------------------------------------------------------
# Do not call install.packages() here — shinyapps.io installs deps at deploy time.
library(shiny)
library(DT)
library(tidyverse)
library(simr)
library(lme4)
library(lmerTest)
library(pwr)
library(scales)

`%||%` <- function(x, y) if (is.null(x)) y else x

# -----------------------------------------------------------------------------
# 2. Statistical helpers (methodology preserved)
# -----------------------------------------------------------------------------

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
  list(data = out, n_pairs = length(unique(out[[pair_col]])), label = label)
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
    if (!is.null(progress_cb)) progress_cb(i, length(effect_sizes), b)
    m <- model
    fixef(m)[param] <- b
    sim <- powerSim(m, test = fixed(test_var), nsim = nsim, progress = FALSE)
    n_sig <- sum(sim$pval < sim$alpha, na.rm = TRUE)
    n_tot <- length(sim$pval)
    ci <- binom.test(x = n_sig, n = n_tot, conf.level = 0.95)$conf.int
    data.frame(
      effect_size = b,
      power = mean(sim$pval < sim$alpha, na.rm = TRUE),
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

# -----------------------------------------------------------------------------
# 3. Theme / CSS / JS (AVA UI shell)
# -----------------------------------------------------------------------------

sa_css <- HTML("
@import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap');

:root {
  --bg: #F7F8FA;
  --surface: #FFFFFF;
  --surface-muted: #F3F4F6;
  --text: #111827;
  --muted: #667085;
  --border: #E5E7EB;
  --border-strong: #D0D5DD;
  --accent: #AB0520;
  --accent-soft: #FCE8EB;
  --accent-hover: #8B0419;
  --success: #047857;
  --success-soft: #ECFDF5;
  --info: #1D4ED8;
  --info-soft: #EFF6FF;
  --radius-sm: 8px;
  --radius-md: 12px;
  --radius-lg: 16px;
  --shadow-sm: 0 1px 2px rgba(16, 24, 40, 0.04);
  --shadow-md: 0 4px 12px rgba(16, 24, 40, 0.06);
  --sidebar-w: 300px;
  --font: Inter, system-ui, -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
}

*, *::before, *::after { box-sizing: border-box; }

html, body {
  margin: 0; padding: 0; height: 100%;
  font-family: var(--font); font-size: 15px; line-height: 1.55;
  color: var(--text); background: var(--bg);
  -webkit-font-smoothing: antialiased;
}

body.sa-app .container-fluid {
  padding: 0 !important; margin: 0 !important; max-width: none !important;
}

.app-shell {
  display: grid;
  grid-template-columns: var(--sidebar-w) 1fr;
  min-height: 100vh;
  background: var(--bg);
}

.sidebar {
  position: sticky; top: 0; height: 100vh; overflow-y: auto;
  background: var(--surface); border-right: 1px solid var(--border);
  padding: 20px 16px 28px;
  display: flex; flex-direction: column; gap: 18px; z-index: 40;
}

.workspace { min-width: 0; display: flex; flex-direction: column; min-height: 100vh; }

.topbar {
  position: sticky; top: 0; z-index: 30;
  min-height: 72px;
  background: rgba(247, 248, 250, 0.92);
  backdrop-filter: blur(10px);
  border-bottom: 1px solid var(--border);
  display: flex; align-items: center; justify-content: space-between; gap: 16px;
  padding: 14px 28px;
}

.content-area { padding: 24px 28px 48px; max-width: 1180px; width: 100%; }

.brand { display: flex; align-items: center; gap: 12px; padding: 4px 8px 8px; }
.brand-mark {
  width: 36px; height: 36px; border-radius: 10px;
  background: linear-gradient(145deg, #AB0520 0%, #8B0419 100%);
  color: white; font-weight: 700; font-size: 11px; letter-spacing: 0.02em;
  display: grid; place-items: center; flex-shrink: 0;
}
.brand-title { font-size: 15px; font-weight: 700; margin: 0; line-height: 1.2; }
.brand-sub { margin: 2px 0 0; font-size: 12px; color: var(--muted); }

.nav-section-label {
  font-size: 11px; font-weight: 600; letter-spacing: 0.06em;
  text-transform: uppercase; color: var(--muted); padding: 0 10px; margin: 0 0 8px;
}
.nav-list { display: flex; flex-direction: column; gap: 2px; }
.nav-item {
  appearance: none; border: none; background: transparent; width: 100%;
  text-align: left; display: flex; align-items: center; gap: 10px;
  padding: 10px 12px; border-radius: var(--radius-sm);
  color: var(--muted); font-size: 14px; font-weight: 500; cursor: pointer;
}
.nav-item:hover { background: var(--surface-muted); color: var(--text); }
.nav-item.active { background: var(--accent-soft); color: var(--accent); font-weight: 600; }

.control-group {
  background: var(--surface-muted); border: 1px solid var(--border);
  border-radius: var(--radius-md); padding: 14px;
}
.control-group + .control-group { margin-top: 10px; }
.control-group-title {
  font-size: 12px; font-weight: 600; letter-spacing: 0.04em;
  text-transform: uppercase; color: var(--muted); margin: 0 0 10px;
}

.sidebar .form-group { margin-bottom: 12px; }
.sidebar label.control-label, .sidebar .control-label {
  font-size: 13px !important; font-weight: 600 !important;
  color: var(--text) !important; margin-bottom: 4px !important;
}
.sidebar .help-block { font-size: 12px; color: var(--muted); margin-top: 4px; }
.sidebar .form-control, .sidebar select.form-control, .sidebar input.form-control {
  border-radius: var(--radius-sm) !important;
  border: 1px solid var(--border-strong) !important;
  box-shadow: none !important; font-size: 13px !important;
  min-height: 38px; background: var(--surface) !important; color: var(--text) !important;
}
.sidebar .form-control:focus {
  border-color: var(--accent) !important;
  box-shadow: 0 0 0 3px rgba(171, 5, 32, 0.15) !important;
}
.sidebar .checkbox { font-size: 13px; }

.sidebar-actions { display: grid; gap: 8px; margin-top: 4px; }

.btn-proxy-run, #run_btn {
  appearance: none; border: none; background: var(--accent); color: white !important;
  font-weight: 600; font-size: 14px; border-radius: var(--radius-sm);
  padding: 11px 16px; width: 100%; cursor: pointer;
  box-shadow: var(--shadow-sm);
}
.btn-proxy-run:hover, #run_btn:hover { background: var(--accent-hover); }

.btn-ava-secondary, .btn-outline-sa {
  appearance: none; display: inline-flex; align-items: center; justify-content: center;
  border: 1px solid var(--border-strong) !important;
  background: var(--surface) !important; color: var(--text) !important;
  font-weight: 600; font-size: 13px; border-radius: var(--radius-sm) !important;
  padding: 10px 14px !important; text-decoration: none !important; cursor: pointer;
}
.btn-ava-secondary:hover, .btn-outline-sa:hover { background: var(--surface-muted) !important; }

.topbar-title {
  margin: 0; font-size: 18px; font-weight: 700; letter-spacing: -0.02em;
}
.topbar-sub { margin: 2px 0 0; font-size: 13px; color: var(--muted); }
.topbar-right { display: flex; align-items: center; gap: 10px; flex-shrink: 0; }
.topbar-run #run_btn { width: auto; min-width: 160px; padding: 10px 18px; }

.meta-pill {
  display: inline-flex; align-items: center; gap: 4px;
  padding: 3px 8px; border-radius: 999px;
  background: var(--surface); border: 1px solid var(--border);
  font-size: 11px; color: var(--muted);
}
.meta-pill.is-ready {
  background: var(--success-soft); border-color: #A7F3D0; color: var(--success);
}

.card {
  background: var(--surface); border: 1px solid var(--border);
  border-radius: var(--radius-lg); box-shadow: var(--shadow-sm);
  padding: 18px 20px; margin-top: 0.5rem;
}
.card-kicker {
  font-size: 11px; font-weight: 700; letter-spacing: 0.06em;
  text-transform: uppercase; color: var(--muted); margin-bottom: 6px;
}
.card h3 { margin: 0 0 12px; font-size: 16px; font-weight: 700; }

.toolbar-row {
  display: flex; align-items: center; justify-content: space-between;
  flex-wrap: wrap; gap: 0; margin: 0; padding: 0;
}
.toolbar-left { display: flex; align-items: center; gap: 0; margin: 0; padding: 0; }
.toolbar-right { display: flex; align-items: center; gap: 0.5rem; margin: 0; }
.toolbar-left .btn-outline-sa:first-child {
  border-top-right-radius: 0 !important; border-bottom-right-radius: 0 !important; margin: 0;
}
.toolbar-left .btn-outline-sa:last-child {
  border-top-left-radius: 0 !important; border-bottom-left-radius: 0 !important;
  margin: 0; margin-left: -1px;
}

.empty-state {
  max-width: 720px; margin: 24px auto; text-align: left;
  background: var(--surface); border: 1px solid var(--border);
  border-radius: 20px; box-shadow: var(--shadow-md); padding: 36px;
}
.empty-state h1 { margin: 0 0 8px; font-size: 26px; font-weight: 700; }
.empty-state .lead { margin: 0; color: var(--muted); font-size: 15px; max-width: 48ch; }

.info-card {
  padding: 16px 18px; background: var(--info-soft); border: 1px solid #BFDBFE;
  border-radius: var(--radius-lg); margin-bottom: 16px;
}
.info-card p { margin: 0; color: #1E3A8A; font-size: 14px; }

#main_pages.tabbable > .nav-tabs { display: none !important; }
#main_pages .tab-content > .tab-pane { padding: 0; }

.shiny-notification {
  border-radius: var(--radius-md) !important;
  border: 1px solid var(--border) !important;
  box-shadow: var(--shadow-md) !important;
  font-family: var(--font) !important;
}

/* Progress overlay (AVA-style) */
.progress-overlay {
  position: fixed; inset: 0;
  background: rgba(17, 24, 39, 0.42);
  display: none; align-items: center; justify-content: center;
  z-index: 1000; padding: 24px;
}
.progress-overlay.is-visible { display: flex; }
.progress-panel {
  width: min(420px, 100%);
  background: var(--surface);
  border-radius: 20px;
  box-shadow: var(--shadow-md);
  padding: 28px;
}
.progress-panel h2 {
  margin: 0 0 6px; font-size: 20px; font-weight: 700;
}
.progress-panel .sub {
  margin: 0 0 20px; color: var(--muted); font-size: 13px;
}
.progress-stages { list-style: none; margin: 0; padding: 0; display: grid; gap: 10px; }
.progress-stages li {
  display: flex; align-items: center; gap: 10px;
  font-size: 14px; color: var(--muted);
}
.progress-stages li .mark {
  width: 18px; height: 18px; border-radius: 999px;
  border: 1.5px solid var(--border-strong);
  display: grid; place-items: center; font-size: 10px; flex-shrink: 0;
}
.progress-stages li.done { color: var(--text); }
.progress-stages li.done .mark {
  background: var(--success-soft); border-color: var(--success); color: var(--success);
}
.progress-stages li.active { color: var(--text); font-weight: 600; }
.progress-stages li.active .mark {
  border-color: var(--accent);
  box-shadow: 0 0 0 3px rgba(171, 5, 32, 0.15);
  animation: sa-pulse 1.2s ease infinite;
}
@keyframes sa-pulse {
  0%, 100% { opacity: 1; }
  50% { opacity: 0.4; }
}

.mobile-nav-toggle {
  display: none; appearance: none; border: 1px solid var(--border-strong);
  background: var(--surface); border-radius: var(--radius-sm);
  width: 40px; height: 40px; cursor: pointer; font-size: 18px;
}
.sidebar-backdrop {
  display: none; position: fixed; inset: 0;
  background: rgba(17,24,39,0.35); z-index: 35;
}

@media (max-width: 900px) {
  .app-shell { grid-template-columns: 1fr; }
  .sidebar {
    position: fixed; left: 0; top: 0; transform: translateX(-105%);
    transition: transform 0.22s ease; width: min(320px, 88vw);
    box-shadow: var(--shadow-md);
  }
  .sidebar.is-open { transform: translateX(0); }
  .sidebar-backdrop.is-open { display: block; }
  .mobile-nav-toggle { display: inline-grid; place-items: center; }
  .topbar, .content-area { padding-left: 16px; padding-right: 16px; }
}
")

sa_js <- HTML("
(function() {
  function setPage(page) {
    Shiny.setInputValue('nav_page', page, {priority: 'event'});
    document.querySelectorAll('.nav-item').forEach(function(el) {
      el.classList.toggle('active', el.getAttribute('data-page') === page);
    });
    var sidebar = document.querySelector('.sidebar');
    var backdrop = document.querySelector('.sidebar-backdrop');
    if (sidebar) sidebar.classList.remove('is-open');
    if (backdrop) backdrop.classList.remove('is-open');
  }

  document.addEventListener('click', function(e) {
    var nav = e.target.closest('.nav-item');
    if (nav) { e.preventDefault(); setPage(nav.getAttribute('data-page')); return; }
    var proxy = e.target.closest('#sidebar_run');
    if (proxy) {
      e.preventDefault();
      var runBtn = document.getElementById('run_btn');
      if (runBtn) runBtn.click();
      return;
    }
    var toggle = e.target.closest('#mobile_nav_toggle');
    if (toggle) {
      document.querySelector('.sidebar').classList.toggle('is-open');
      document.querySelector('.sidebar-backdrop').classList.toggle('is-open');
      return;
    }
    var backdrop = e.target.closest('.sidebar-backdrop');
    if (backdrop) {
      document.querySelector('.sidebar').classList.remove('is-open');
      backdrop.classList.remove('is-open');
    }
  });

  function setProgressStage(activeIdx) {
    var items = document.querySelectorAll('#sa_progress_stages li');
    items.forEach(function(li, i) {
      li.classList.remove('done', 'active');
      var mark = li.querySelector('.mark');
      if (i < activeIdx) {
        li.classList.add('done');
        if (mark) mark.textContent = '\u2713';
      } else if (i === activeIdx) {
        li.classList.add('active');
        if (mark) mark.textContent = '';
      } else if (mark) {
        mark.textContent = '';
      }
    });
  }

  Shiny.addCustomMessageHandler('saProgress', function(msg) {
    var overlay = document.getElementById('sa_progress_overlay');
    if (!overlay) return;
    var nStages = document.querySelectorAll('#sa_progress_stages li').length;
    if (msg.visible) {
      overlay.classList.add('is-visible');
      overlay.setAttribute('aria-hidden', 'false');
      var idx = typeof msg.stage === 'number' ? msg.stage : 0;
      setProgressStage(idx);
      var sub = document.getElementById('sa_progress_sub');
      if (sub && msg.detail) sub.textContent = msg.detail;
    } else {
      setProgressStage(nStages);
      setTimeout(function() {
        overlay.classList.remove('is-visible');
        overlay.setAttribute('aria-hidden', 'true');
        setProgressStage(0);
      }, 350);
    }
  });

  $(document).on('shiny:connected', function() { setPage('analysis'); });
})();
")

nav_button <- function(id, label) {
  tags$button(
    id = id, class = "nav-item", type = "button",
    `data-page` = gsub("^nav_", "", id),
    tags$span(class = "nav-label", label)
  )
}

# -----------------------------------------------------------------------------
# 4. UI
# -----------------------------------------------------------------------------

ui <- fluidPage(
  title = "Sensitivity Analysis",
  tags$head(
    tags$meta(name = "viewport", content = "width=device-width, initial-scale=1"),
    tags$style(sa_css),
    tags$script(sa_js)
  ),
  class = "sa-app",

  div(
    class = "app-shell",
    div(class = "sidebar-backdrop", `aria-hidden` = "true"),

    tags$aside(
      class = "sidebar",
      `aria-label` = "Sensitivity analysis controls",

      div(
        class = "brand",
        div(class = "brand-mark", "SA"),
        div(
          class = "brand-text",
          tags$p(class = "brand-title", "Sensitivity Analysis"),
          tags$p(class = "brand-sub", "Mixed-model power curves")
        )
      ),

      div(
        tags$p(class = "nav-section-label", "Navigate"),
        div(
          class = "nav-list",
          nav_button("nav_analysis", "Analysis"),
          nav_button("nav_help", "Help")
        )
      ),

      div(
        tags$p(class = "nav-section-label", "Configure"),

        div(
          class = "control-group",
          tags$p(class = "control-group-title", "Dataset"),
          fileInput("csv_file", label = NULL, accept = c(".csv", "text/csv"),
                    buttonLabel = "Upload CSV", placeholder = "No file selected"),
          checkboxInput("strings_as_factors", "Read character columns as factors", TRUE),
          actionButton("toggle_preview", "Show data preview",
                       class = "btn-ava-secondary", width = "100%"),
          uiOutput("dataset_status")
        ),

        div(
          class = "control-group",
          tags$p(class = "control-group-title", "Clustering & variables"),
          selectInput("pair_id", "Pair / cluster ID (L2)", choices = NULL),
          selectInput("exclude_ids", "Exclude ID columns", choices = NULL, multiple = TRUE),
          selectInput("predictor", "Focal predictor (IV)", choices = NULL),
          selectInput("outcomes", "Focal outcomes (DVs)", choices = NULL, multiple = TRUE),
          selectInput("covariates", "Covariates", choices = NULL, multiple = TRUE),
          selectInput("factor_vars", "Treat as factors", choices = NULL, multiple = TRUE),
          selectInput("standardize", "Standardize (creates _std)", choices = NULL, multiple = TRUE),
          helpText("Factors stay categorical. Only numeric columns should be standardized.")
        ),

        div(
          class = "control-group",
          tags$p(class = "control-group-title", "Simulation"),
          textInput("effect_sizes", "Effect sizes (comma-separated β)",
                    value = "0.05, 0.10, 0.15, 0.20"),
          numericInput("nsim", "Simulations per effect size (nsim)",
                       value = 10, min = 1, max = 5000, step = 1),
          helpText("Use ~10 to test; ≥1000 for reporting."),
          numericInput("seed", "Random seed", value = 1995, min = 1, step = 1)
        ),

        div(
          class = "sidebar-actions",
          tags$button(id = "sidebar_run", type = "button", class = "btn-proxy-run",
                      "Run sensitivity analysis")
        )
      )
    ),

    div(
      class = "workspace",
      tags$header(
        class = "topbar",
        div(
          style = "display:flex; align-items:center; gap:12px; min-width:0;",
          tags$button(id = "mobile_nav_toggle", class = "mobile-nav-toggle",
                      type = "button", `aria-label` = "Open navigation", "☰"),
          div(
            tags$h1(class = "topbar-title", "Sensitivity Analysis"),
            tags$p(class = "topbar-sub",
                   "Upload clustered/twin data, then simulate power across effect sizes.")
          )
        ),
        div(
          class = "topbar-right",
          uiOutput("topbar_status"),
          div(class = "topbar-run",
              actionButton("run_btn", "Run analysis", class = "btn-proxy-run"))
        )
      ),

      div(
        class = "content-area",
        tabsetPanel(
          id = "main_pages",
          type = "hidden",
          selected = "analysis",

          tabPanel(
            title = NULL, value = "analysis",
            conditionalPanel(
              condition = "!output.has_data_flag",
              div(
                class = "empty-state",
                tags$h1("Start a sensitivity analysis"),
                tags$p(
                  class = "lead",
                  "Upload a clustered or twin CSV, assign the pair ID, IV, outcomes, and covariates, then run power simulations across effect sizes."
                )
              )
            ),
            conditionalPanel(
              condition = "output.has_data_flag",
              conditionalPanel(
                condition = "input.toggle_preview % 2 == 1",
                div(
                  class = "card",
                  div(class = "card-kicker", "Data"),
                  tags$h3("Data preview"),
                  verbatimTextOutput("data_summary"),
                  DTOutput("data_preview")
                )
              ),
              div(
                class = "toolbar-row",
                div(
                  class = "toolbar-left",
                  actionButton("toggle_curves", "Show sensitivity curves",
                               class = "btn-outline-sa"),
                  actionButton("toggle_results", "Show results table",
                               class = "btn-outline-sa")
                ),
                div(
                  class = "toolbar-right",
                  downloadButton("download_results", "Download results CSV",
                                 class = "btn-ava-secondary"),
                  downloadButton("download_plot", "Download plot PNG",
                                 class = "btn-ava-secondary")
                )
              ),
              uiOutput("curves_panel"),
              uiOutput("results_panel")
            )
          ),
          tabPanel(
            title = NULL, value = "help",
            div(
              class = "card",
              div(class = "card-kicker", "Guide"),
              tags$h3("How to use this tool"),
              tags$ol(
                tags$li("Upload your twin/cluster CSV."),
                tags$li("Choose the cluster ID; exclude person IDs from modeling."),
                tags$li("Set the focal IV, one or more DVs, and covariates."),
                tags$li("Mark categorical columns as factors; standardize continuous variables as needed."),
                tags$li("Start with a small nsim (e.g., 10). Use ≥1000 for manuscript figures.")
              ),
              tags$p(
                tags$strong("Example setup: "),
                "IV = mv; DVs = connection, GrimAgeAA, DunedinPACE; ",
                "covariates = all other non-ID columns."
              )
            )
          )
        )
      )
    )
  ),

  # Analysis progress overlay (AVA-style)
  div(
    id = "sa_progress_overlay",
    class = "progress-overlay",
    `aria-hidden` = "true",
    role = "dialog",
    `aria-labelledby` = "sa_progress_title",
    div(
      class = "progress-panel",
      tags$h2(id = "sa_progress_title", "Running sensitivity analysis"),
      tags$p(id = "sa_progress_sub", class = "sub", "Preparing analysis…"),
      tags$ul(
        id = "sa_progress_stages",
        class = "progress-stages",
        tags$li(`data-stage` = "prep",
                tags$span(class = "mark"), tags$span("Preparing data")),
        tags$li(`data-stage` = "std",
                tags$span(class = "mark"), tags$span("Standardizing variables")),
        tags$li(`data-stage` = "fit",
                tags$span(class = "mark"), tags$span("Fitting mixed models")),
        tags$li(`data-stage` = "sim",
                tags$span(class = "mark"), tags$span("Running power simulations")),
        tags$li(`data-stage` = "results",
                tags$span(class = "mark"), tags$span("Preparing results"))
      )
    )
  )
)

# -----------------------------------------------------------------------------
# 5. Server
# -----------------------------------------------------------------------------

server <- function(input, output, session) {

  results_rv <- reactiveVal(NULL)
  plot_rv <- reactiveVal(NULL)
  log_rv <- reactiveVal("")
  curves_open <- reactiveVal(FALSE)
  results_open <- reactiveVal(FALSE)

  append_log <- function(...) {
    log_rv(paste0(log_rv(), paste0(...), "\n"))
  }

  observeEvent(input$nav_page, {
    page <- input$nav_page
    if (!is.null(page) && page %in% c("analysis", "help")) {
      updateTabsetPanel(session, "main_pages", selected = page)
    }
  }, ignoreInit = FALSE)

  has_data <- reactive(isTruthy(input$csv_file))
  output$has_data_flag <- reactive(has_data())
  outputOptions(output, "has_data_flag", suspendWhenHidden = FALSE)

  raw_data <- reactive({
    req(input$csv_file)
    read.csv(input$csv_file$datapath, stringsAsFactors = input$strings_as_factors)
  })

  observeEvent(raw_data(), {
    df <- raw_data()
    cols <- names(df)
    num_cols <- guess_numeric_cols(df)
    fac_guess <- guess_factor_cols(df)

    default_pair <- if ("pairid" %in% cols) "pairid" else if ("pair_id" %in% cols) "pair_id" else cols[1]
    default_pred <- if ("mv" %in% cols) "mv" else setdiff(num_cols, default_pair)[1]
    default_out <- intersect(c("connection", "GrimAgeAA", "DunedinPACE"), cols)
    id_guess <- intersect(
      c("studyid", "pair_id", "pairid", "fid", "id", "twinid", "person_id"), cols
    )
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
    for (v in facs) df[[v]] <- as.factor(df[[v]])
    df
  })

  output$dataset_status <- renderUI({
    if (!has_data()) {
      return(tags$p(style = "margin:8px 0 0; font-size:12px; color:var(--muted);",
                    "Upload a CSV to begin."))
    }
    df <- prepared_data()
    div(
      style = "margin-top:8px;",
      tags$span(class = "meta-pill is-ready", "Data ready"),
      tags$span(class = "meta-pill", sprintf("%s rows", format(nrow(df), big.mark = ","))),
      tags$span(class = "meta-pill", sprintf("%s vars", ncol(df)))
    )
  })

  output$topbar_status <- renderUI({
    if (!has_data()) {
      return(tags$span(class = "meta-pill", "No dataset"))
    }
    tagList(
      tags$span(class = "meta-pill is-ready", "Dataset loaded"),
      if (!is.null(results_rv())) tags$span(class = "meta-pill is-ready", "Results ready") else NULL
    )
  })

  observeEvent(input$toggle_preview, {
    showing <- (input$toggle_preview %% 2 == 1)
    updateActionButton(
      session, "toggle_preview",
      label = if (showing) "Hide data preview" else "Show data preview"
    )
  }, ignoreInit = TRUE)

  sync_toggle_labels <- function() {
    updateActionButton(
      session, "toggle_curves",
      label = if (curves_open()) "Hide sensitivity curves" else "Show sensitivity curves"
    )
    updateActionButton(
      session, "toggle_results",
      label = if (results_open()) "Hide results table" else "Show results table"
    )
  }

  observeEvent(input$toggle_curves, {
    if (is.null(plot_rv())) {
      showNotification("Run the analysis first to see curves.", type = "message")
      return()
    }
    curves_open(!curves_open())
    sync_toggle_labels()
  }, ignoreInit = TRUE)

  observeEvent(input$toggle_results, {
    if (is.null(results_rv())) {
      showNotification("Run the analysis first to see the results table.", type = "message")
      return()
    }
    results_open(!results_open())
    sync_toggle_labels()
  }, ignoreInit = TRUE)

  output$curves_panel <- renderUI({
    req(curves_open(), plot_rv())
    div(
      class = "card",
      div(class = "card-kicker", "Results"),
      tags$h3("Sensitivity curves"),
      div(style = "width:100%; min-height:520px;",
          plotOutput("power_plot", height = "520px", width = "100%")),
      verbatimTextOutput("run_log")
    )
  })

  output$results_panel <- renderUI({
    req(results_open(), results_rv())
    div(
      class = "card",
      div(class = "card-kicker", "Results"),
      tags$h3("Results table"),
      DTOutput("results_table")
    )
  })

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
    datatable(prepared_data(), options = list(pageLength = 8, scrollX = TRUE), rownames = FALSE)
  })

  output$run_log <- renderText(log_rv())

  sa_progress <- function(stage, detail) {
    session$sendCustomMessage(
      "saProgress",
      list(visible = TRUE, stage = as.integer(stage), detail = detail)
    )
  }

  observeEvent(input$run_btn, {
    log_rv("")
    results_rv(NULL)
    plot_rv(NULL)
    curves_open(FALSE)
    results_open(FALSE)
    sync_toggle_labels()

    df0 <- prepared_data()
    pair_col <- input$pair_id
    predictor_raw <- input$predictor
    outcomes_raw <- input$outcomes
    cov_raw <- input$covariates
    std_vars <- input$standardize %||% character(0)
    exclude <- input$exclude_ids %||% character(0)

    if (is.null(pair_col) || !nzchar(pair_col)) {
      showNotification("Select a pair / cluster ID.", type = "error"); return()
    }
    if (is.null(predictor_raw) || !nzchar(predictor_raw)) {
      showNotification("Select a focal predictor.", type = "error"); return()
    }
    if (is.null(outcomes_raw) || length(outcomes_raw) < 1) {
      showNotification("Select at least one outcome.", type = "error"); return()
    }

    effect_sizes <- suppressWarnings(
      as.numeric(trimws(unlist(strsplit(input$effect_sizes, ","))))
    )
    if (!length(effect_sizes) || any(!is.finite(effect_sizes))) {
      showNotification("Enter valid numeric effect sizes.", type = "error"); return()
    }
    if (is.null(input$nsim) || input$nsim < 1) {
      showNotification("nsim must be ≥ 1.", type = "error"); return()
    }

    cov_raw <- setdiff(unique(c(predictor_raw, cov_raw)), c(pair_col, exclude, outcomes_raw))
    std_vars <- setdiff(std_vars, c(pair_col, exclude))
    std_vars <- intersect(std_vars, names(df0))
    std_vars <- std_vars[vapply(std_vars, function(v) is.numeric(df0[[v]]), logical(1))]

    on.exit(session$sendCustomMessage("saProgress", list(visible = FALSE)), add = TRUE)
    sa_progress(0L, "Preparing data…")

    withProgress(message = "Running sensitivity analysis", value = 0, {
      sa_progress(1L, "Standardizing variables…")
      append_log("Standardizing: ", paste(std_vars, collapse = ", "))
      df <- add_std_cols(df0, std_vars)

      pred <- std_name(predictor_raw, std_vars)
      covs <- unname(vapply(cov_raw, std_name, character(1), standardize_vars = std_vars))
      outcome_specs <- lapply(outcomes_raw, function(o) {
        list(outcome_raw = o, outcome_var = std_name(o, std_vars), label = o)
      })

      append_log("Predictor in models: ", pred)
      append_log("Covariates in models: ", paste(covs, collapse = ", "))

      set.seed(input$seed)
      total_steps <- max(1, length(outcome_specs) * length(effect_sizes))
      step_box <- new.env(parent = emptyenv())
      step_box$i <- 0

      sa_progress(2L, "Fitting mixed models…")

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
            df = filtered$data, outcome_var = spec$outcome_var,
            predictor = pred, covariates = covs, pair_col = pair_col
          ),
          error = function(e) {
            append_log(sprintf("[%s] Model error: %s", spec$label, e$message))
            NULL
          }
        )
        if (is.null(model)) return(NULL)
        if (!pred %in% names(fixef(model))) {
          append_log(sprintf("[%s] Predictor '%s' not in fixed effects.", spec$label, pred))
          return(NULL)
        }

        sa_progress(
          3L,
          sprintf("Running power simulations… %s", spec$label)
        )

        power_curve_simr(
          model = model, effect_sizes = effect_sizes, nsim = input$nsim,
          param = pred, test_var = pred,
          progress_cb = function(i, n, b) {
            step_box$i <- step_box$i + 1
            setProgress(
              value = step_box$i / total_steps,
              detail = sprintf("%s | β = %s (%d/%d)", spec$label, b, i, n)
            )
            sa_progress(
              3L,
              sprintf("Running power simulations… %s | β = %s (%d/%d)",
                      spec$label, b, i, n)
            )
          }
        ) %>% mutate(outcome = spec$label)
      })

      if (is.null(power_results) || nrow(power_results) == 0) {
        append_log("No results produced.")
        showNotification("No results produced.", type = "warning")
        return()
      }

      sa_progress(4L, "Preparing results…")

      power_results <- power_results %>%
        mutate(outcome = factor(outcome, levels = outcomes_raw))

      n_panels <- length(unique(power_results$outcome))
      facet_nrow <- if (n_panels <= 3) 1L else 2L

      fig <- ggplot(power_results, aes(x = effect_size, y = power)) +
        geom_hline(yintercept = 0.80, linetype = "dashed",
                   color = "gray50", linewidth = 0.5) +
        geom_ribbon(aes(ymin = lower, ymax = upper),
                    alpha = 0.2, fill = "#AB0520") +
        geom_line(color = "#AB0520", linewidth = 0.8) +
        geom_point(color = "#AB0520", size = 2) +
        facet_wrap(~ outcome, nrow = facet_nrow) +
        scale_y_continuous(
          limits = c(0, 1), breaks = seq(0, 1, by = 0.10),
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
          panel.grid.minor = element_line(color = "gray90", linewidth = 0.3)
        )

      results_rv(power_results)
      plot_rv(fig)
      curves_open(TRUE)
      results_open(TRUE)
      sync_toggle_labels()
      append_log("Done.")
      showNotification("Sensitivity analysis finished.", type = "message")
    })
  })

  output$power_plot <- renderPlot({
    req(plot_rv())
    print(plot_rv())
  }, height = 520, res = 96)

  output$results_table <- renderDT({
    req(results_rv())
    datatable(
      results_rv() %>%
        mutate(power = round(power, 3), lower = round(lower, 3), upper = round(upper, 3)),
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

# -----------------------------------------------------------------------------
# 6. App launch
# -----------------------------------------------------------------------------
shinyApp(ui = ui, server = server)
