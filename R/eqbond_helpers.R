# eqbond_helpers.R
# Shared functions for the equity-bond relative performance reversion work.
# Source via:
#   source("https://raw.githubusercontent.com/RWLab/rwlab-helpers/master/R/eqbond_helpers.R")
#
# Two notebooks use these and they must agree:
#
#   research/equity_bond_flow_effects/   does the effect exist, and how strong
#   trading/equity_bond_eom_performance_reversion/   can you keep any of it
#
# If the signal were written out in both, they would drift, and the drift would
# be silent: the research chart and the simulation would gradually describe
# different strategies while both continued to run.
#
# Signal construction lives here. Cost assumptions, trade buffers and anything
# else about how a particular person implements it stay in the notebooks.

# ── Data preparation ────────────────────────────────────────────────────────

#' Filter to the traded instruments and add trading day of month
#'
#' Drops months with fewer than `min_trading_days` observations, since a
#' partial month has no meaningful "first half" to measure.
#'
#' @param prices Price data with date, ticker, close, closeadjusted
#' @param tickers Instruments to keep
#' @param min_trading_days Months shorter than this are dropped
#' @return date, ticker, close, closeadjusted, tdm, log_return
prepare_eqbond_prices <- function(prices,
                                  tickers = c("VTI", "TLT"),
                                  min_trading_days = 15) {
  p <- prices %>%
    dplyr::filter(ticker %in% tickers) %>%
    dplyr::mutate(year = lubridate::year(date), month = lubridate::month(date)) %>%
    dplyr::arrange(date)

  incomplete <- p %>%
    dplyr::group_by(ticker, year, month) %>%
    dplyr::summarise(trading_days = dplyr::n(), .groups = "keep") %>%
    dplyr::filter(trading_days < min_trading_days)

  p %>%
    dplyr::anti_join(incomplete, by = c("ticker", "year", "month")) %>%
    dplyr::group_by(ticker, year, month) %>%
    dplyr::mutate(tdm = dplyr::row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::select(date, ticker, close, closeadjusted, tdm) %>%
    dplyr::group_by(ticker) %>%
    dplyr::arrange(date) %>%
    dplyr::mutate(log_return = log(closeadjusted / dplyr::lag(closeadjusted, n = 1))) %>%
    dplyr::ungroup() %>%
    stats::na.omit()
}

# ── Relative performance ────────────────────────────────────────────────────

#' Equity minus bond return over each half of each month
#'
#' The original notebook bucketed `tdm` into "weeks of the month" and then
#' collapsed those buckets into a two-way split at day 15. The intermediate
#' variable never affected anything, so the split point is stated directly and
#' is a parameter.
#'
#' @param prices Output of prepare_eqbond_prices()
#' @param split_day Trading day of month dividing first half from second
#' @return year, month, part_mnth_1, part_mnth_2
eqbond_relative_performance <- function(prices, split_day = 15) {
  prices %>%
    tidyr::pivot_wider(id_cols = c(date, tdm), names_from = ticker,
                       values_from = log_return) %>%
    stats::na.omit() %>%
    tidyr::pivot_longer(cols = c(VTI, TLT), names_to = "ticker",
                        values_to = "log_return") %>%
    dplyr::mutate(
      month_split = dplyr::if_else(tdm <= split_day, 1, 2),
      month = lubridate::month(date),
      year  = lubridate::year(date)
    ) %>%
    dplyr::group_by(ticker, year, month, month_split) %>%
    dplyr::summarise(partial_return = sum(log_return), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = ticker, values_from = partial_return) %>%
    dplyr::mutate(eq_bond_outperf = log(1 + (exp(VTI) - 1) - (exp(TLT) - 1))) %>%
    tidyr::pivot_wider(id_cols = c(year, month), names_from = month_split,
                       names_prefix = "part_mnth_", values_from = eq_bond_outperf)
}

# ── Signal ──────────────────────────────────────────────────────────────────

#' Positions from the reversion rule
#'
#' Reversion: whichever instrument underperformed over the first half of the
#' month is bought for the second half. Long only, one leg at a time, flat for
#' the first half.
#'
#' @param prices Output of prepare_eqbond_prices()
#' @param split_day Trading day of month dividing the halves
#' @param signal_threshold Relative return must exceed this to take a side
#' @param legs "both", "equity_only" or "bond_only"
#' @return prices plus part_mnth_1 and position (0 or 1)
eqbond_signal <- function(prices,
                          split_day = 15,
                          signal_threshold = 0,
                          legs = c("both", "equity_only", "bond_only")) {
  legs <- match.arg(legs)
  take_equity <- legs %in% c("both", "equity_only")
  take_bond   <- legs %in% c("both", "bond_only")

  rel <- eqbond_relative_performance(prices, split_day = split_day)

  prices %>%
    dplyr::mutate(year = lubridate::year(date), month = lubridate::month(date)) %>%
    dplyr::left_join(rel, by = c("month", "year")) %>%
    dplyr::mutate(position = dplyr::case_when(
      take_equity & ticker == "VTI" & tdm > split_day & part_mnth_1 < -signal_threshold ~ 1,
      take_bond   & ticker == "TLT" & tdm > split_day & part_mnth_1 >  signal_threshold ~ 1,
      TRUE ~ 0
    )) %>%
    dplyr::select(date, ticker, close, closeadjusted, tdm, log_return,
                  part_mnth_1, position)
}

# ── Frictionless returns ────────────────────────────────────────────────────

#' What the signal would have earned with no costs and no friction
#'
#' The research question is whether the effect exists at all. If it does not
#' survive zero costs there is nothing to implement, and this answers that
#' without a simulator.
#'
#' Deliberately naive: positions are taken at the close and earn the next day's
#' return, with no commissions, no spread, no trade buffer, no interest on idle
#' cash and no position sizing beyond 0 or 1.
#'
#' **This is a ceiling, not a forecast.** How much of it survives a real cost
#' model is the implementation notebook's question, not this one's.
#'
#' @param signal_df Output of eqbond_signal()
#' @return date, strategy_return, cumulative_return
eqbond_frictionless_returns <- function(signal_df) {
  signal_df %>%
    dplyr::group_by(ticker) %>%
    dplyr::arrange(date) %>%
    dplyr::mutate(fwd_log_return = dplyr::lead(log_return)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(contribution = position * fwd_log_return) %>%
    dplyr::group_by(date) %>%
    dplyr::summarise(strategy_return = sum(contribution, na.rm = TRUE),
                     .groups = "drop") %>%
    dplyr::arrange(date) %>%
    dplyr::mutate(cumulative_return = cumsum(strategy_return))
}
