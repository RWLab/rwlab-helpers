# vix_helpers.R
# Shared functions for the UVXY-VXZ spread.
# Source via:
#   source("https://raw.githubusercontent.com/RWLab/rwlab-helpers/master/R/vix_helpers.R")
#
# Three consumers have to agree about these, and they are in three places:
#
#   research/vxz_uvxy_spread_hedge_engineering.ipynb   where the hedge came from
#   trading/uvxy_vxz_spread/    the simulation, and the script that sizes today
#   your own operational layer  whatever computes the position you actually hold
#
# The hedge ratio in particular was written out twice already, once in the
# notebook as an inline sigmoid and once in uvxy-vxz-dynamic-hedge.R. Two copies
# of a formula drift, and the drift is silent: the simulation and the live
# position gradually describe different strategies while both keep running.

# ── The hedge ratio ─────────────────────────────────────────────────────────

#' VIX-aware dynamic hedge ratio
#'
#' Blends a hedge ratio smoothly across three volatility regimes using sigmoid
#' transitions centred on two VIX thresholds. Continuous and differentiable
#' rather than a hard regime switch, so the position does not jump when VIX
#' crosses a boundary.
#'
#' @param vix_level Numeric. VIX level, scalar or vector.
#' @param b1 Numeric. Boundary between the low and mid regimes (e.g. 18).
#' @param b2 Numeric. Boundary between the mid and high regimes (e.g. 28).
#'   Must be greater than `b1`.
#' @param h_low,h_mid,h_high Numeric. Hedge ratio in each regime.
#' @param k Numeric. Steepness of the transitions. Larger is sharper.
#' @return Numeric vector of blended hedge ratios, same length as `vix_level`.
#'   Each value is a convex combination of the three regime ratios.
#' @examples
#' vix_aware_hedge(22.3, b1 = 18, b2 = 28, h_low = 3.7, h_mid = 2.7,
#'                 h_high = 3.2, k = 1)
vix_aware_hedge <- function(vix_level, b1, b2, h_low, h_mid, h_high, k) {
  if (b1 >= b2) stop("`b1` must be less than `b2`")

  sigmoid <- function(x) 1 / (1 + exp(-x))

  s1 <- sigmoid(k * (vix_level - b1))
  s2 <- sigmoid(k * (vix_level - b2))

  (1 - s1) * h_low + (s1 - s2) * h_mid + s2 * h_high
}

# ── Spread weights ──────────────────────────────────────────────────────────

#' Target weights for the UVXY-VXZ spread
#'
#' Short one unit of UVXY against `hedge_ratio` units of VXZ, then normalise so
#' gross exposure is exactly 1. The sleeve's budget scales it from there, so the
#' weights published here are at the strategy's own native scale.
#'
#' Two different things produce a missing hedge ratio, and they need different
#' treatment. **Before the signal starts** -- while a rolling regression window
#' is still filling, or before the VIX series begins -- there is genuinely no
#' position, and those dates are dropped so they do not pad the NAV series with
#' flat days that dilute the volatility estimate. **After it starts**, a gap
#' means something upstream failed, usually a join. Zeroing that quietly takes
#' the book flat for a day and pays a round trip to do it, and it looks exactly
#' like a deliberate flat, so it errors instead.
#'
#' @param df Long data frame with date, ticker, and a hedge_ratio column
#' @param short_ticker Instrument held short (one unit before normalisation)
#' @param hedge_col Name of the column holding the hedge ratio
#' @return `df` from the first live date onward, with a `weight` column summing
#'   to gross 1 on every date
uvxy_vxz_weights <- function(df, short_ticker = "UVXY", hedge_col = "hedge_ratio") {
  df <- dplyr::arrange(df, date)

  is_hedged <- function(d) d$ticker != short_ticker

  live <- df$date[is_hedged(df) & !is.na(df[[hedge_col]])]
  if (length(live) == 0) {
    stop("hedge ratio is NA on every date. Check the join that supplies it.")
  }

  df <- df[df$date >= min(live), ]

  # Recomputed against the filtered frame. Carrying a logical vector across the
  # filter and indexing it with the new frame's rows silently misaligns them.
  gaps <- unique(df$date[is_hedged(df) & is.na(df[[hedge_col]])])
  if (length(gaps) > 0) {
    stop(sprintf(
      "hedge ratio missing on %d date(s) after the signal starts, first %s. Zeroing those would take the book flat for a day, pay a round trip to do it, and look like a deliberate flat. Fix the upstream join.",
      length(gaps), min(gaps)
    ))
  }

  df %>%
    dplyr::mutate(weight = dplyr::if_else(ticker == short_ticker, -1, .data[[hedge_col]])) %>%
    dplyr::group_by(date) %>%
    dplyr::mutate(weight = weight / sum(abs(weight))) %>%
    dplyr::ungroup()
}

# ── Splits ──────────────────────────────────────────────────────────────────

#' UVXY reverse splits
#'
#' `ratio` is old shares per new share, so a 1-for-5 reverse split is 5.
#'
#' This table has to be maintained by hand, and UVXY reverse-splits every year
#' or two. It lives here rather than in the notebook so that updating it once
#' reaches every member on their next run, rather than each of them having to
#' notice a split and edit their own copy.
#'
#' `uvxy_unadjusted_close()` checks whether it has gone stale.
UVXY_SPLITS <- data.frame(
  date  = as.Date(c("2018-09-18", "2021-05-21", "2023-06-03",
                    "2024-04-11", "2025-11-20")),
  ratio = c(5, 10, 10, 5, 5)
)

#' Reconstruct UVXY's traded (non split-adjusted) closing price
#'
#' `rsims` needs split-adjusted prices for P&L and raw prices for commissions,
#' because commissions are charged per share and the share count is the one you
#' actually traded. For UVXY the two differ by a factor of over a thousand, so
#' using the wrong one does not produce a small error.
#'
#' **Why this is checked rather than trusted.** If a split is missing from the
#' table, the whole history before it is divided by too small a factor and comes
#' out uniformly inflated. Nothing about that looks wrong: the series is still
#' smooth, still the right shape, and the simulation still runs. It just charges
#' the wrong commissions, which is precisely what the trade buffer is calibrated
#' against.
#'
#' The detector uses the fund's own behaviour. ProShares reverse-splits UVXY to
#' keep its nominal price in a tradeable band, so a reconstructed price far above
#' that band early in the sample means a split is missing, and one far below it at
#' the end of the sample means a split has happened or is due.
#'
#' @param prices_df UVXY rows with date and close (split-adjusted)
#' @param splits Split table, defaults to `UVXY_SPLITS`
#' @param plausible Range the raw price is expected to stay within
#' @param on_stale "warn" or "stop"
#' @return `prices_df` with `cum_split_factor` and `unadjustedclose`
uvxy_unadjusted_close <- function(prices_df,
                                  splits = UVXY_SPLITS,
                                  plausible = c(4, 250),
                                  on_stale = c("warn", "stop")) {
  on_stale <- match.arg(on_stale)
  prices_df <- dplyr::arrange(prices_df, date)

  cum_split_factor <- rep(1, nrow(prices_df))
  for (i in seq_len(nrow(splits))) {
    before <- prices_df$date < splits$date[i]
    cum_split_factor[before] <- cum_split_factor[before] * splits$ratio[i]
  }

  out <- prices_df %>%
    dplyr::mutate(cum_split_factor = cum_split_factor,
                  unadjustedclose  = close / cum_split_factor)

  complain <- if (on_stale == "stop") stop else warning
  px <- out$unadjustedclose[!is.na(out$unadjustedclose)]
  n  <- length(px)

  # The primary detector. A reverse-splitting fund is held in a price band on
  # purpose, so the start and the end of the sample should sit at comparable
  # levels. A missing split inflates everything before it by the split ratio,
  # and the smallest ratio UVXY has ever done is 5, so a step of that size
  # stands well clear of ordinary drift.
  if (n >= 50) {
    fifth <- max(1, floor(n / 5))
    early <- stats::median(utils::head(px, fifth))
    late  <- stats::median(utils::tail(px, fifth))
    if (is.finite(early) && is.finite(late) && late > 0 && early / late > 4) {
      complain(sprintf(
        "Reconstructed UVXY price starts the sample around $%.2f and ends it around $%.2f, a factor of %.1f. A reverse-splitting fund is kept in a band, so a split is probably missing from the table. Last one recorded: %s.",
        early, late, early / late, max(splits$date)))
    }
  }

  # Backstops, for the cases the ratio test cannot see: an inflation that is
  # already present at the start of the sample, and a split that happened after
  # the last recorded one.
  hi <- max(px)
  if (hi > plausible[2]) {
    complain(sprintf(
      "Reconstructed UVXY price reaches $%.0f, above the $%.0f a reverse-splitting fund would be expected to stay under. Check the split table against the fund's history.",
      hi, plausible[2]))
  }
  recent <- utils::tail(px, 20)
  if (length(recent) > 0 && mean(recent) < plausible[1]) {
    complain(sprintf(
      "Reconstructed UVXY price is averaging $%.2f at the end of the sample. UVXY reverse-splits before it gets that low, so a split has probably happened since %s and is not in the table.",
      mean(recent), max(splits$date)))
  }

  out
}
