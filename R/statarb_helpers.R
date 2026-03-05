# statarb_helpers.R
# Shared utility functions for the stat-arb research notebooks
# Source via: source("https://raw.githubusercontent.com/RWLab/rwlab-helpers/main/R/statarb_helpers.R")
#
# Contains only generic, non-proprietary functions:
#   - Data loading and processing
#   - Portfolio math helpers
#   - Commission presets
#
# Signal construction, weight calculation, and strategy logic stay in the notebooks.

# ── Data Processing ──────────────────────────────────────────────────────────

#' Flatten raw spread data to bidirectional ticker-level observations
#'
#' Each pair appears twice with opposite z-scores and displacements,
#' so that downstream aggregation naturally considers both legs.
process_raw_spreads <- function(raw_df) {
  df <- raw_df %>%
    mutate(date = as.Date(date)) %>%
    arrange(date, ticker) %>%
    rename("stock1" = ticker)

  # legacy alpha field
  if ("alpha" %in% colnames(df)) {
    df <- rename(df, "displacement" = alpha)
  }

  df_rev <- df %>%
    mutate(.s1 = stock1, .s2 = stock2) %>%
    mutate(
      stock1 = .s2,
      stock2 = .s1,
      zscore = -zscore,
      displacement = -displacement
    ) %>%
    select(-.s1, -.s2)

  df <- bind_rows(df, df_rev) %>%
    filter(stock1 != stock2) %>%
    arrange(date, stock1) %>%
    rename("ticker" = stock1)

  df
}

#' Load all stat-arb datasets from GCS
#'
#' Returns a named list with spreads, prices (with forward returns),
#' acquisitions, and short borrow data.
load_statarb_data <- function(pod = "EquityFactors",
                              short_borrow_years = 2020:2026,
                              auth_method = c("auto", "colab", "local"),
                              local_json_path = NULL) {
  auth_method <- match.arg(auth_method)

  # Auth
 if (auth_method == "colab" || (auth_method == "auto" && exists("rwlab_data_auth"))) {
    rwRtools::rwlab_data_auth()
  } else if (!is.null(local_json_path)) {
    googleCloudStorageR::gcs_auth(json_file = local_json_path)
  }
  # else: assume gcloud ADC is already configured

  # Load core datasets
  spreads <- rwRtools::load_lab_object(pod = pod, object = "statarb/spreads.feather")
  if ("alpha" %in% colnames(spreads)) {
    spreads <- rename(spreads, "displacement" = alpha)
  }

  prices <- rwRtools::load_lab_object(pod = pod, object = "statarb/prices.feather")
  acquisitions <- rwRtools::load_lab_object(pod = pod, object = "statarb/acquisitions.feather")

  # Short borrow data
  short_borrow <- short_borrow_years %>%
    purrr::map_dfr(~ rwRtools::macro_get_historical_short_sale(year = .x)) %>%
    arrange(date, ticker)

  # Add returns to prices
  prices <- prices %>%
    arrange(date) %>%
    group_by(ticker) %>%
    mutate(
      log_return = log(closeadj / lag(closeadj)),
      fwd_log_return = lead(log_return),
      fwd2_log_return = lead(log_return, 2)
    ) %>%
    ungroup()

  list(
    spreads = spreads,
    prices = prices,
    acquisitions = acquisitions,
    short_borrow = short_borrow
  )
}

# ── Short Borrow Aggregation ────────────────────────────────────────────────

#' Aggregate historical short borrow rates by ticker
#'
#' Returns named vectors for mean, p75, and max borrow rates.
#' Tickers in the universe with no borrow data are filled with the
#' universe mean for each metric.
aggregate_short_borrows <- function(short_borrow, universe_tickers) {
  agg <- short_borrow %>%
    filter(ticker %in% universe_tickers) %>%
    group_by(ticker) %>%
    summarise(
      mean_borrow = mean(-rebate_rate / 100, na.rm = TRUE),
      max_borrow = max(-rebate_rate / 100, na.rm = TRUE),
      p75_borrow = quantile(-rebate_rate / 100, 0.75, na.rm = TRUE),
      p90_borrow = quantile(-rebate_rate / 100, 0.90, na.rm = TRUE),
      .groups = "drop"
    )

  # Fill missing tickers with universe mean
  missing <- setdiff(universe_tickers, agg$ticker)
  if (length(missing) > 0) {
    defaults <- agg %>%
      summarise(
        mean_borrow = mean(mean_borrow, na.rm = TRUE),
        max_borrow = mean(max_borrow, na.rm = TRUE),
        p75_borrow = mean(p75_borrow, na.rm = TRUE),
        p90_borrow = mean(p90_borrow, na.rm = TRUE)
      )
    fill_rows <- tibble(ticker = missing) %>% bind_cols(defaults)
    agg <- bind_rows(agg, fill_rows)

    warning(
      glue::glue(
        "{length(missing)} tickers missing from borrow data, filled with universe mean: ",
        "{paste(missing, collapse = ', ')}"
      )
    )
  }

  list(
    mean = deframe(agg %>% select(ticker, mean_borrow)),
    max = deframe(agg %>% select(ticker, max_borrow)),
    p75 = deframe(agg %>% select(ticker, p75_borrow)),
    p90 = deframe(agg %>% select(ticker, p90_borrow))
  )
}

# ── Simulation Helpers ───────────────────────────────────────────────────────

#' Iterative weight capping preserving long/short balance
#'
#' Caps absolute weights at max_weight while maintaining the original
#' gross exposure on each side (long/short) via iterative rescaling.
apply_max_weight_ls <- function(df, max_weight = 0.05, max_iter = 20) {
  for (i in seq_len(max_iter)) {
    if (!any(abs(df$weight) > max_weight + 1e-9, na.rm = TRUE)) break

    df <- df %>%
      mutate(side = sign(weight)) %>%
      filter(side != 0) %>%
      group_by(date, side) %>%
      mutate(
        original_gross = sum(abs(weight)),
        capped = sign(weight) * pmin(abs(weight), max_weight),
        weight = capped * original_gross / sum(abs(capped))
      ) %>%
      ungroup() %>%
      bind_rows(df %>% filter(sign(weight) == 0))
  }

  if (any(abs(df$weight) > max_weight + 1e-9, na.rm = TRUE)) {
    warning("Max iterations reached, some weights still exceed cap")
  }

  df
}

#' Calculate annualised Sharpe ratio from rsims backtest results
calc_sharpe <- function(backtest_results, capitalise_profits = TRUE, initial_cash = 0) {
  if(capitalise_profits == FALSE & initial_cash == 0) {
    stop("If capitalise_profits is FALSE, initial_cash must be > 0.")
  }

  if(capitalise_profits == TRUE) {
    sharpe <- backtest_results %>%
      rsims::calc_port_returns() %>%
      summarise(
        ann_return = 252 * mean(returns, na.rm = TRUE),
        ann_vol = sqrt(252) * sd(returns, na.rm = TRUE),
        sharpe = ann_return / ann_vol
      ) %>%
      pull(sharpe)
    
    return(sharpe)
  } else {
    sharpe <- backtest_results %>%
      group_by(date) %>%
      summarise(port_value = sum(exposure)) %>%
      mutate(port_return = (port_value - lag(port_value))/initial_cash) %>%
      summarise(
        port_annualised_return = mean(port_return, na.rm = TRUE)*trading_days_per_year,
        port_vol = sd(port_return, na.rm = TRUE)*sqrt(trading_days_per_year),
        sharpe = sqrt(trading_days_per_year)* mean(port_return, na.rm = TRUE)/sd(port_return, na.rm = TRUE)
      ) %>% 
      pull(sharpe)

      return(sharpe)
  }
    
}

#' Extract NAV time series from rsims backtest results
append_nav_to_bt_results <- function(bt_results_df) {
  port_nav <- bt_results_df %>%
    select(date, exposure, ticker) %>%
    group_by(date) %>%
    summarise(exposure = sum(exposure, na.rm = TRUE), .groups = "drop") %>%
    mutate(ticker = "NAV")

  bt_results_df %>%
    select(date, exposure, ticker) %>%
    bind_rows(port_nav) %>%
    arrange(date)
}

# ── Commission Presets ───────────────────────────────────────────────────────

commissions_cost_free <- list(
  fun = rsims::us_tiered_commission,
  max_pct_per_order = 0.0,
  min_dollars_per_order = 0.0,
  dollars_per_share = 0.0
)

commissions_ib_tiered_high <- list(
  fun = rsims::us_tiered_commission,
  max_pct_per_order = 0.01,
  min_dollars_per_order = 0.35,
  dollars_per_share = 0.0035
)

commissions_ib_tiered_mid <- list(
  fun = rsims::us_tiered_commission,
  max_pct_per_order = 0.01,
  min_dollars_per_order = 0.35,
  dollars_per_share = 0.0015
)

commissions_ib_tiered_low <- list(
  fun = rsims::us_tiered_commission,
  max_pct_per_order = 0.01,
  min_dollars_per_order = 0.35,
  dollars_per_share = 0.0005
)

commissions_ib_fixed <- list(
  fun = rsims::us_tiered_commission,
  max_pct_per_order = 0.01,
  min_dollars_per_order = 1.0,
  dollars_per_share = 0.005
)
