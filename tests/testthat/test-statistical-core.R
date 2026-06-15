## Phase 1 statistical core, validated against the pencil-and-paper example
## vector (mean = 6, n = 20, S = 120) used in the original development script.

cs <- c(1, 1, 1, 1, 1, 2, 2, 2, 3, 3,
        4, 4, 5, 6, 7, 8, 10, 12, 20, 27)

test_that("check_clade_sizes rejects bad input via the public functions", {
  expect_error(geom_expectation(c(1, 2, NA)), "NA")
  expect_error(geom_expectation(c(0, 1, 2)), "positive")
  expect_error(geom_expectation(c(1, 2.5, 3)), "whole")
  expect_silent(invisible(geom_expectation(cs)))
})

test_that("geom_expectation matches hand-computed values", {
  g <- geom_expectation(cs)
  expect_equal(g$n_clades, 20L)
  expect_equal(g$total_species, 120)
  expect_equal(g$mean_size, 6)
  expect_equal(g$geom_p, 1 / 6)
  expect_equal(g$expected_sd, sqrt(30))               # sqrt(mean*(mean-1))
  expect_equal(g$expected_prop_singleton, 1 / 6)
  expect_equal(g$prop_singleton, 0.25)                # 5 singletons / 20
  expect_equal(g$simpson_div, 1 - sum((cs / 120)^2))  # Gini-Simpson, inline (no vegan)
  expect_equal(g$observed_sd, stats::sd(cs))
})

test_that("clade_tests geometric tails and flags are correct", {
  r <- clade_tests(cs)
  expect_equal(nrow(r), 20L)
  expect_true(all(c("p_lower", "p_upper", "too_few", "too_many") %in% names(r)))
  expect_false(any(c("p_lower_adj", "p_upper_adj") %in% names(r))) # none requested
  expect_equal(unique(r$mean_clade_size), 6)
  ## sorted smallest to largest
  expect_equal(r$clade_sizes, sort(cs))
  ## raw p_upper = P(K >= k) = (1 - 1/6)^(k - 1)
  expect_equal(r$p_upper, (5 / 6)^(r$clade_sizes - 1))
  ## only the two largest clades (20, 27) exceed the geometric upper tail
  expect_equal(sum(r$too_many), 2L)
  expect_equal(sort(r$clade_sizes[r$too_many]), c(20, 27))
  ## no clade is significantly too small under the marginal geometric
  expect_equal(sum(r$too_few), 0L)
})

test_that("clade_tests conditional tails match the exact broken-stick formula", {
  ## Small hand-computable case: S = 4, n = 2 -> compositions (1,3),(2,2),(3,1).
  ## Marginal of a single part is uniform on {1, 2, 3}, so
  ##   P(K >= k) = choose(S - k, n - 1) / choose(S - 1, n - 1) = (4 - k) / 3.
  r <- clade_tests(c(1, 3), method = "conditional")
  expect_equal(attr(r, "method"), "conditional")
  ## sorted smallest -> largest: size 1 then size 3
  expect_equal(r$clade_sizes, c(1, 3))
  expect_equal(r$p_upper, c(1, 1 / 3))      # P(K >= 1) = 1 ; P(K >= 3) = 1/3
  expect_equal(r$p_lower, c(1 / 3, 1))      # P(K <= 1) = 1/3 ; P(K <= 3) = 1
})

test_that("clade_tests conditional behaves on the worked example", {
  r <- clade_tests(cs, method = "conditional")
  expect_equal(nrow(r), 20L)
  ## all tail probabilities are valid
  expect_true(all(r$p_lower >= 0 & r$p_lower <= 1))
  expect_true(all(r$p_upper >= 0 & r$p_upper <= 1))
  ## upper tail is monotone non-increasing in clade size (sorted ascending)
  expect_false(is.unsorted(rev(r$p_upper)))
  ## the two largest clades are still flagged as too_many, as under geometric
  expect_true(all(c(20, 27) %in% r$clade_sizes[r$too_many]))
  expect_equal(sum(r$too_few), 0L)
})

test_that("clade_tests method is validated and geometric is the default", {
  expect_equal(attr(clade_tests(cs), "method"), "geometric")
  expect_error(clade_tests(cs, method = "bogus"), "should be one of")
})

test_that("clade_tests conditional is degenerate for a single clade", {
  r <- clade_tests(7, method = "conditional")
  expect_equal(r$p_lower, 1)
  expect_equal(r$p_upper, 1)
  expect_false(r$too_few)
  expect_false(r$too_many)
})

test_that("clade_tests adjust= adds adjusted columns and is more conservative", {
  r_raw <- clade_tests(cs)
  r_bh  <- clade_tests(cs, adjust = "BH")
  expect_true(all(c("p_lower_adj", "p_upper_adj") %in% names(r_bh)))
  expect_true(all(r_bh$p_upper_adj >= r_bh$p_upper))          # adjusted >= raw
  expect_lte(sum(r_bh$too_many), sum(r_raw$too_many))         # no more discoveries
  expect_error(clade_tests(cs, adjust = "not_a_method"), "adjust must be one of")
})

test_that("sd_sim_test is reproducible and reports the observed SD", {
  s1 <- sd_sim_test(cs, nsim = 2000, seed = 42)
  s2 <- sd_sim_test(cs, nsim = 2000, seed = 42)
  expect_equal(s1, s2)                                  # seed -> reproducible
  expect_equal(s1$n_clades, 20L)
  expect_equal(s1$total_species, 120)
  expect_equal(s1$observed_sd, stats::sd(cs))
  expect_true(s1$p_sd_too_large > 0 && s1$p_sd_too_large <= 1)
  expect_true(s1$p_sd_too_small > 0 && s1$p_sd_too_small <= 1)
  expect_true(s1$p_two_sided <= 1)
})

test_that("sd_sim_test warns on a degenerate null", {
  expect_warning(sd_sim_test(c(1, 1, 1), nsim = 10), "degenerate")
})

test_that("clade_rank_data ranks largest-first with the geometric prediction", {
  rd <- clade_rank_data(cs)
  expect_equal(nrow(rd), 20L)
  expect_equal(rd$clade_sizes, sort(cs, decreasing = TRUE))
  expect_equal(rd$clade_rank, seq_len(20))
  expect_equal(rd$clade_sizes[1], 27)                   # rank 1 = largest
  expect_equal(rd$ln_clade_rank, log(rd$clade_rank))
  expect_equal(rd$ln_rank_predicted,
               log(20) + (rd$clade_sizes - 1) * log(5 / 6))
})
