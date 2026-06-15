## Phase 2 compute_pulse_score: structure, the direction= argument, and the
## consistency relationships that must hold by construction.

skip_if_not_installed("ape")
skip_if_not_installed("phytools")

set.seed(42)
tree <- ape::rcoal(60)
ages <- c(0.2, 0.4, 0.6, 0.8)

test_that("compute_pulse_score validates inputs", {
  expect_error(compute_pulse_score("nope", ages), "phylo")
  expect_error(compute_pulse_score(tree, numeric(0)), "non-empty")
  expect_error(compute_pulse_score(tree, c(0.2, NA)), "positive and non-NA")
  expect_error(compute_pulse_score(tree, c(0.2, -1)), "positive and non-NA")
  expect_error(compute_pulse_score(tree, ages, direction = "sideways"))
})

test_that("output has one row per tip in original order with expected columns", {
  ps <- compute_pulse_score(tree, ages)
  expect_equal(nrow(ps), ape::Ntip(tree))
  expect_equal(ps$species, tree$tip.label)         # original tip order preserved
  expect_true(all(c("species", "n_slices", "n_flagged",
                    "pulse_fraction", "pulse_score") %in% names(ps)))
  ## one clade_ + one flagged_ column per (sorted, unique) age
  for (a in ages) {
    expect_true(paste0("clade_", a, "Ma")   %in% names(ps))
    expect_true(paste0("flagged_", a, "Ma") %in% names(ps))
  }
  expect_identical(attr(ps, "direction"), "oversized")
})

test_that("score bookkeeping is internally consistent", {
  ps <- compute_pulse_score(tree, ages)
  expect_true(all(ps$n_slices == length(ages)))
  ## pulse_fraction == n_flagged / n_slices
  expect_equal(ps$pulse_fraction, ps$n_flagged / ps$n_slices)
  expect_true(all(ps$pulse_fraction >= 0 & ps$pulse_fraction <= 1))
  ## score is 0 exactly when never flagged, positive otherwise; always finite
  expect_true(all(is.finite(ps$pulse_score)))
  expect_true(all((ps$pulse_score == 0) == (ps$n_flagged == 0)))
  expect_true(all(ps$pulse_score >= 0))
})

test_that("duplicate / unsorted ages collapse to the sorted unique set", {
  ps <- compute_pulse_score(tree, c(0.4, 0.2, 0.4, 0.6))
  expect_true(all(ps$n_slices == 3))
  expect_true(all(c("clade_0.2Ma", "clade_0.4Ma", "clade_0.6Ma") %in% names(ps)))
})

test_that("direction = 'undersized' scores the lower tail instead", {
  pu <- compute_pulse_score(tree, ages, direction = "undersized")
  expect_identical(attr(pu, "direction"), "undersized")
  expect_equal(nrow(pu), ape::Ntip(tree))
  expect_true(all(is.finite(pu$pulse_score)))
  expect_true(all((pu$pulse_score == 0) == (pu$n_flagged == 0)))
})

test_that("adjust = 'BH' is accepted and never increases the flag count", {
  raw <- compute_pulse_score(tree, ages, adjust = "none")
  bh  <- compute_pulse_score(tree, ages, adjust = "BH")
  expect_equal(nrow(bh), ape::Ntip(tree))
  ## BH correction can only make flags more conservative, never less
  expect_true(sum(bh$n_flagged) <= sum(raw$n_flagged))
})
