# rp_helpers.R
# Shared estimators for the risk premia strategies.
# Source via:
#   source("https://raw.githubusercontent.com/RWLab/rwlab-helpers/master/R/rp_helpers.R")
#
# These were written out twice: once in the RP7 simulation notebook and once in
# the RP7 dashboard. Two copies of an estimator drift, and the drift is silent
# -- the weights you simulate and the weights that get published gradually
# describe different strategies while both keep running.

# ── EWMA volatility ─────────────────────────────────────────────────────────

#' Exponentially weighted moving average of a vector
#'
#' Uses the CURRENT observation in each update, seeded with the first value.
#' Note the contrast with `ewma_cor()` below, which does not.
#'
#' @param x Numeric vector
#' @param lambda Decay. Higher means slower, more weight on history.
#' @return Numeric vector, same length as `x`
ewma <- function(x, lambda) {
  out <- vector(mode = "double", length = length(x))
  out[1] <- x[1]
  for (i in 2:length(x)) {
    out[i] <- (1 - lambda) * x[i] + lambda * out[i - 1]
  }
  out
}

# ── EWMA correlation ────────────────────────────────────────────────────────

#' EWMA pairwise correlation estimate
#'
#' Covariance and both variances are tracked with the same decay, and the
#' correlation is their ratio. `lambda` follows the RiskMetrics convention, so
#' higher values weight history more.
#'
#' **This estimate is one period stale, and deliberately left that way.** The
#' update uses the previous observation rather than the current one, so the
#' value at index `i` is what the correct recursion would produce at `i - 1`.
#' Verified exactly: the two series have a correlation of 1.0000 at a one-day
#' offset, and both recover the true correlation on simulated data, so it is a
#' lag rather than a bias, and it errs towards being stale rather than
#' peeking.
#'
#' It is preserved because the published RP7 dashboard computes it the same
#' way. Changing it here alone would make members' simulations disagree with
#' the weights they are given. Worth noting that `ewma()` above does use the
#' current observation, so a strategy combining the two is mixing a volatility
#' estimate through `t` with a correlation estimate through `t-1`.
#'
#' @param x,y Numeric vectors of returns, equal length, NAs removed
#' @param lambda Decay
#' @param initialisation_wdw Observations used to seed, returned as NA
#' @return Numeric vector, same length as `x`
ewma_cor <- function(x, y, lambda, initialisation_wdw = 100) {
  if (length(x) != length(y)) stop("x and y must be the same length")
  if (length(x) <= initialisation_wdw) {
    stop(sprintf("need more than initialisation_wdw (%d) observations, got %d",
                 initialisation_wdw, length(x)))
  }
  if (lambda <= 0 || lambda >= 1) stop("lambda must be between 0 and 1")
  if (anyNA(x) || anyNA(y)) stop("x and y must not contain NA")

  init_x <- x[1:initialisation_wdw]
  init_y <- y[1:initialisation_wdw]
  num_obs <- length(x)

  old_cov   <- stats::cov(init_x, init_y)
  old_var_x <- stats::var(init_x)
  old_var_y <- stats::var(init_y)
  old_x     <- mean(init_x)
  old_y     <- mean(init_y)

  out <- vector(mode = "numeric", length = num_obs)
  out[1:initialisation_wdw] <- NA

  for (i in c((initialisation_wdw + 1):num_obs)) {
    this_cov   <- lambda * old_cov   + (1 - lambda) * (old_x * old_y)
    this_var_x <- lambda * old_var_x + (1 - lambda) * old_x^2
    this_var_y <- lambda * old_var_y + (1 - lambda) * old_y^2

    out[i] <- this_cov / (sqrt(this_var_x) * sqrt(this_var_y))

    old_cov   <- this_cov
    old_var_x <- this_var_x
    old_var_y <- this_var_y
    old_x     <- x[i]      # previous observation on the next pass: see above
    old_y     <- y[i]
  }

  out
}

# ── Correlation matrix recovery ─────────────────────────────────────────────

#' Rebuild a correlation matrix from a long frame of pairwise correlations
#'
#' Expects the lower triangle including the diagonal, ordered as
#' `lower.tri(..., diag = TRUE)` fills it.
#'
#' @param long_cors Data frame with an `ewma_cor` column
#' @param tickers Character vector naming the assets
#' @param num_assets Length of `tickers`
#' @return Symmetric numeric matrix
recover_cormat <- function(long_cors, tickers, num_assets) {
  cor_mat <- matrix(rep(0, num_assets * num_assets), num_assets)
  dimnames(cor_mat) <- list(tickers, tickers)
  cor_mat[lower.tri(cor_mat, diag = TRUE)] <- long_cors$ewma_cor
  cor_mat[upper.tri(cor_mat)] <- t(cor_mat)[upper.tri(cor_mat)]
  cor_mat
}

# ── Mean pairwise correlation, per asset per day ────────────────────────────

#' Each asset's mean pairwise correlation, and its distance from the average
#'
#' The original built a correlation matrix per day inside a loop that filtered
#' the whole frame on every pass, which is quadratic in the number of days.
#' This computes the same numbers by summing the long frame directly. Identical
#' output, and it turns minutes into seconds.
#'
#' **The diagonal is included**, matching the original's `rowMeans()` over a
#' matrix whose diagonal is 1. That does not change the ranking, because every
#' asset gets the same constant, but it does scale the resulting delta by
#' `(n-1)/n`. So `cor_multiplier` means slightly different things at different
#' universe sizes: at seven assets the effective multiplier is 6/7 of the
#' number you set, and at four assets it would be 3/4. Set
#' `include_diagonal = FALSE` to remove that dependence, which changes results.
#'
#' @param ewma_cors Long frame: date, tickers ("A, B"), ewma_cor
#' @param include_diagonal Keep the original's treatment (TRUE) or drop self-pairs
#' @return date, ticker, mean_pw_cor, port_mean_pw_cor, pw_cor_delta
mean_pairwise_correlation <- function(ewma_cors, include_diagonal = TRUE) {
  long <- ewma_cors %>%
    tidyr::separate(tickers, into = c("ticker_x", "ticker_y"), sep = ",\\s*")

  self_pairs  <- long %>% dplyr::filter(ticker_x == ticker_y)
  cross_pairs <- long %>% dplyr::filter(ticker_x != ticker_y)

  # Each cross pair contributes to both of its assets.
  contributions <- dplyr::bind_rows(
    cross_pairs %>% dplyr::select(date, ticker = ticker_x, ewma_cor),
    cross_pairs %>% dplyr::select(date, ticker = ticker_y, ewma_cor)
  )

  if (include_diagonal) {
    contributions <- dplyr::bind_rows(
      contributions,
      self_pairs %>% dplyr::select(date, ticker = ticker_x, ewma_cor)
    )
  }

  contributions %>%
    dplyr::group_by(date, ticker) %>%
    dplyr::summarise(total = sum(ewma_cor), n_terms = dplyr::n(), .groups = "drop") %>%
    dplyr::group_by(date) %>%
    dplyr::mutate(
      mean_pw_cor      = total / n_terms,
      port_mean_pw_cor = mean(mean_pw_cor),
      pw_cor_delta     = port_mean_pw_cor - mean_pw_cor
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(date, ticker, mean_pw_cor, port_mean_pw_cor, pw_cor_delta)
}

# ── Hard weight constraints ─────────────────────────────────────────────────

#' Apply weight caps and redistribute the freed weight
#'
#' Two kinds of cap, in any combination:
#'
#' **Per-asset**, via `max_single`. A single number applies to every asset; a
#' named vector sets different limits for different assets, and anything
#' unnamed is uncapped.
#'
#' **Per-group**, via `groups`. Each entry names a set of tickers and a limit on
#' their combined absolute weight. Any number of groups, and an asset may
#' belong to more than one.
#'
#' When a cap binds, the freed weight is redistributed across the names still
#' free to receive it, in proportion to their current weight, so the book stays
#' fully invested. **This has to iterate**, and that is the part a hand-rolled
#' version usually gets wrong. Two ways to get it wrong:
#'
#' Renormalising everything at the end by dividing through by the total scales
#' the capped group back up along with everything else, undoing the very cap
#' you just applied. A group held to 35% of a book summing to 0.75 comes back
#' as 46.7%.
#'
#' Giving the whole deficit to the names outside the group, pro rata, respects
#' the group cap but can push one of those names past its own per-asset cap.
#'
#' The fix for both: only give weight to names that are free to take it, then
#' look again, because some of them will no longer be free.
#'
#' Names already at zero stay at zero. An asset dropped upstream is not
#' resurrected here.
#'
#' @param tickers Character vector, same order as `w`
#' @param w Numeric weights for a single date
#' @param max_single Scalar, or named vector, or NULL for no per-asset cap
#' @param groups List of `list(tickers = , max = )`, or NULL
#' @param max_group,group_tickers Shorthand for a single group. Ignored if
#'   `groups` is supplied.
#' @param max_iter Safety stop on the redistribution loop
#' @param tol Convergence tolerance
#' @return Numeric vector of capped weights
#' @examples
#' tk <- c("EMB", "GLD", "IEF", "TLT", "VEA", "VTI", "VWO")
#' w  <- c(0.18, 0.06, 0.26, 0.16, 0.06, 0.22, 0.06)
#'
#' # RW's published configuration
#' apply_caps_and_normalise(tk, w, max_single = 0.25,
#'                          max_group = 0.35, group_tickers = c("EMB","IEF","TLT"))
#'
#' # Your own: cap bonds harder, cap equities too, and let gold run
#' apply_caps_and_normalise(tk, w,
#'   max_single = c(EMB = 0.15, IEF = 0.15, TLT = 0.15, VTI = 0.30),
#'   groups = list(
#'     bonds    = list(tickers = c("EMB", "IEF", "TLT"), max = 0.30),
#'     equities = list(tickers = c("VEA", "VTI", "VWO"), max = 0.50)
#'   ))
apply_caps_and_normalise <- function(tickers, w,
                                     max_single = NULL,
                                     groups = NULL,
                                     max_group = NULL,
                                     group_tickers = c("EMB", "IEF", "TLT"),
                                     max_iter = 50, tol = 1e-10) {

  # Shorthand for the single-group case.
  if (is.null(groups) && !is.null(max_group)) {
    groups <- list(list(tickers = group_tickers, max = max_group))
  }
  if (is.null(max_single) && length(groups) == 0) return(w)

  n <- length(w)

  # Per-asset limits, as a vector aligned to `tickers`.
  single <- rep(Inf, n)
  if (!is.null(max_single)) {
    if (is.null(names(max_single))) {
      if (length(max_single) != 1 && length(max_single) != n) {
        stop("max_single must be a single number, a named vector, or one value per ticker")
      }
      single <- rep(max_single, length.out = n)
    } else {
      unknown <- setdiff(names(max_single), tickers)
      if (length(unknown) > 0) {
        stop(sprintf("max_single names not in tickers: %s", paste(unknown, collapse = ", ")))
      }
      single[match(names(max_single), tickers)] <- unname(max_single)
    }
  }

  masks <- lapply(groups, function(g) tickers %in% g$tickers)
  limits <- vapply(groups, function(g) g$max, numeric(1))

  for (iter in seq_len(max_iter)) {
    w_prev <- w

    w <- sign(w) * pmin(abs(w), single)

    for (k in seq_along(masks)) {
      s <- sum(abs(w[masks[[k]]]))
      if (s > limits[k]) w[masks[[k]]] <- w[masks[[k]]] * (limits[k] / s)
    }

    total <- sum(abs(w))
    if (total < tol) return(w)
    if (abs(total - 1) < tol) return(w)
    deficit <- 1 - total

    if (deficit < 0) {
      w <- w / total
      next
    }

    abs_w <- abs(w)

    # Free to receive weight: below its own cap, in no group that is at its
    # cap, and not already zero.
    at_single <- abs_w >= single - tol
    in_full_group <- rep(FALSE, n)
    for (k in seq_along(masks)) {
      if (sum(abs_w[masks[[k]]]) >= limits[k] - tol) {
        in_full_group <- in_full_group | masks[[k]]
      }
    }
    flexible <- !at_single & !in_full_group & (abs_w > tol)

    if (!any(flexible)) return(w)

    flex_sum <- sum(abs_w[flexible])
    w[flexible] <- w[flexible] * ((flex_sum + deficit) / flex_sum)

    if (max(abs(w - w_prev)) < tol) return(w)
  }

  warning("apply_caps_and_normalise: max_iter reached without converging")
  w
}
