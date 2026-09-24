# Tests for R/statarb_helpers.R
# Run from the repo root:  Rscript tests/test_statarb_helpers.R
suppressPackageStartupMessages({ library(dplyr); library(tibble) })

this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
root <- normalizePath(file.path(dirname(this_file), ".."))
source(file.path(root, "R", "statarb_helpers.R"))

n_pass <- 0; n_fail <- 0
check <- function(name, expr) {
  res <- tryCatch({ force(expr); TRUE }, error = function(e) conditionMessage(e))
  if (isTRUE(res)) { n_pass <<- n_pass + 1; cat("PASS ", name, "\n") }
  else { n_fail <<- n_fail + 1; cat("FAIL ", name, "\n      ", res, "\n") }
}
expect_equal <- function(actual, expected) {
  if (!isTRUE(all.equal(actual, expected))) stop(sprintf("expected %s, got %s",
    paste(format(expected), collapse = ","), paste(format(actual), collapse = ",")))
}

# ── hold_state(): a two-gate latch ─────────────────────────────────────────
cal <- tibble(trading_date = as.Date("2024-01-01") + 0:5, t_idx = 1:6)

gates <- function(absz, entry, exit, dates = cal$trading_date[seq_along(absz)], ticker = "A") {
  tibble(date = dates, ticker = ticker, absz = absz,
         .in_ok = absz >= entry, .out_ok = absz >= exit)
}
held_of <- function(df) df %>% arrange(ticker, date) %>% pull(held)

check("hold_state: exit == entry reproduces the hard gate exactly", {
  absz <- c(0.6, 0.4, 0.7, 0.2, 0.5)
  out <- hold_state(gates(absz, entry = 0.5, exit = 0.5), cal)
  expect_equal(held_of(out), absz >= 0.5)
})

check("hold_state: a name that cleared the entry bar is carried through the band and dropped below the exit bar", {
  out <- hold_state(gates(c(0.6, 0.4, 0.3, 0.2), entry = 0.5, exit = 0.25), cal)
  expect_equal(held_of(out), c(TRUE, TRUE, TRUE, FALSE))
})

check("hold_state: being inside the band is not enough, the entry bar must have been cleared first", {
  out <- hold_state(gates(c(0.3, 0.4, 0.6, 0.3), entry = 0.5, exit = 0.25), cal)
  expect_equal(held_of(out), c(FALSE, FALSE, TRUE, TRUE))
})

check("hold_state: leaving the universe resets the latch", {
  dates <- cal$trading_date[c(1, 2, 4, 5)]          # absent on day 3
  out <- hold_state(gates(c(0.6, 0.4, 0.3, 0.6), entry = 0.5, exit = 0.25, dates = dates), cal)
  expect_equal(held_of(out), c(TRUE, TRUE, FALSE, TRUE))
})

check("hold_state: an NA gate is not held and breaks the run", {
  out <- hold_state(gates(c(0.6, NA, 0.4, 0.6), entry = 0.5, exit = 0.25), cal)
  expect_equal(held_of(out), c(TRUE, FALSE, FALSE, TRUE))
})

check("hold_state: tickers are latched independently", {
  df <- bind_rows(gates(c(0.6, 0.4), entry = 0.5, exit = 0.25, ticker = "A"),
                  gates(c(0.3, 0.3), entry = 0.5, exit = 0.25, ticker = "B"))
  out <- hold_state(df, cal)
  expect_equal(out %>% filter(ticker == "A") %>% arrange(date) %>% pull(held), c(TRUE, TRUE))
  expect_equal(out %>% filter(ticker == "B") %>% arrange(date) %>% pull(held), c(FALSE, FALSE))
})

check("hold_state: returns the input columns plus `held`, gate and scratch columns removed", {
  out <- hold_state(gates(c(0.6, 0.4), entry = 0.5, exit = 0.25), cal)
  expect_equal(sort(names(out)), sort(c("date", "ticker", "absz", "held")))
  expect_equal(nrow(out), 2L)
})

check("hold_state: a weekend (consecutive trading days, non-consecutive dates) does not reset the latch", {
  wk <- tibble(trading_date = as.Date(c("2024-01-05", "2024-01-08", "2024-01-09")), t_idx = 1:3)  # Fri, Mon, Tue
  out <- hold_state(gates(c(0.6, 0.4, 0.3), entry = 0.5, exit = 0.25, dates = wk$trading_date), wk)
  expect_equal(held_of(out), c(TRUE, TRUE, TRUE))
})

check("hold_state: an exit bar of zero holds a name until it leaves the panel", {
  out <- hold_state(gates(c(0.6, 0.1, 0.0, 0.05), entry = 0.5, exit = 0), cal)
  expect_equal(held_of(out), c(TRUE, TRUE, TRUE, TRUE))
})

check("hold_state: a date missing from the calendar is isolated, never carried into", {
  dates <- c(cal$trading_date[1:2], as.Date("2030-01-01"))
  out <- hold_state(gates(c(0.6, 0.4, 0.4), entry = 0.5, exit = 0.25, dates = dates), cal)
  expect_equal(held_of(out), c(TRUE, TRUE, FALSE))
})

check("hold_state: returns rows in the caller's order, not sorted", {
  df <- gates(c(0.6, 0.4, 0.3), entry = 0.5, exit = 0.25)[c(3, 1, 2), ]
  out <- hold_state(df, cal)
  expect_equal(out$date, df$date)
  expect_equal(out$held, c(TRUE, TRUE, TRUE))
})

check("hold_state: a duplicated calendar date is an error, not silently duplicated rows", {
  bad_cal <- bind_rows(cal, cal[2, ])
  res <- tryCatch({ hold_state(gates(c(0.6, 0.4), entry = 0.5, exit = 0.25), bad_cal); "no error" },
                  error = function(e) "error")
  expect_equal(res, "error")
})

cat(sprintf("\n%d passed, %d failed\n", n_pass, n_fail))
if (n_fail > 0) quit(status = 1)
