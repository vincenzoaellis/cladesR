## Phase 2 tree-slicing core (extract_clade_sizes / build_clade_table /
## ancestor_descendant_stats). Uses a small, fixed, hand-checkable ultrametric
## tree so the expected clades are known exactly and no RNG is involved.

skip_if_not_installed("ape")
skip_if_not_installed("phytools")

## A balanced, ultrametric 4-tip tree of total height 3:
##   ((a:1,b:1):2,(c:1,d:1):2);
## Heights from the root: the two cherries (ab, cd) join at height 2; tips at 3.
## So slicing at age_before_present = 1.5 (root-height 1.5) cuts ABOVE the
## cherry nodes -> two 2-species clades {a,b} and {c,d}.
## Slicing at age_before_present = 0.5 (root-height 2.5) cuts BELOW the cherries
## -> four singleton clades.
fixed_tree <- ape::read.tree(text = "((a:1,b:1):2,(c:1,d:1):2);")

test_that("extract_clade_sizes input validation works", {
  expect_error(extract_clade_sizes("not a tree", 1), "phylo")
  expect_error(extract_clade_sizes(fixed_tree, -1), "non-negative")
  expect_error(extract_clade_sizes(fixed_tree, NA_real_), "NA")
  expect_error(extract_clade_sizes(fixed_tree, c(1, 2)), "single numeric")
  expect_error(extract_clade_sizes(fixed_tree, 100), "outside the age range")
  expect_error(extract_clade_sizes(fixed_tree, 1, ultrametric_option = 3), "either 1 or 2")
})

test_that("extract_clade_sizes recovers the two expected cherries", {
  s <- extract_clade_sizes(fixed_tree, age_before_present = 1.5)
  expect_equal(length(s$clade_sizes), 2L)
  expect_equal(sort(unname(s$clade_sizes)), c(2L, 2L))
  expect_equal(sum(s$clade_sizes), ape::Ntip(fixed_tree))
  ## membership is {a,b} and {c,d} (order-independent check)
  members <- lapply(s$species, sort)
  expect_true(any(vapply(members, identical, logical(1), c("a", "b"))))
  expect_true(any(vapply(members, identical, logical(1), c("c", "d"))))
  ## nodes rooting 2-species clades are internal (> Ntip)
  expect_true(all(s$nodes > ape::Ntip(fixed_tree)))
  ## returned clades object is a multiPhylo
  expect_s3_class(s$clades, "multiPhylo")
})

test_that("extract_clade_sizes age == 0 gives singletons for every tip", {
  s0 <- extract_clade_sizes(fixed_tree, age_before_present = 0)
  expect_equal(length(s0$clade_sizes), ape::Ntip(fixed_tree))
  expect_true(all(s0$clade_sizes == 1L))
  ## singleton nodes are tip nodes (1:Ntip)
  expect_true(all(s0$nodes %in% seq_len(ape::Ntip(fixed_tree))))
})

test_that("slicing near the present yields four singletons", {
  s <- extract_clade_sizes(fixed_tree, age_before_present = 0.5)
  expect_equal(length(s$clade_sizes), 4L)
  expect_true(all(s$clade_sizes == 1L))
})

test_that("build_clade_table assembles a correct species-by-slice table", {
  ages   <- c(1.5, 0.5)              # supplied youngest-last to test re-ordering
  slices <- lapply(ages, function(a) extract_clade_sizes(fixed_tree, a))
  tab    <- build_clade_table(slices, ages)

  ## one row per species, species col + one col per slice
  expect_equal(nrow(tab), ape::Ntip(fixed_tree))
  expect_equal(names(tab), c("species", "clade_1.5Ma", "clade_0.5Ma"))
  expect_setequal(tab$species, fixed_tree$tip.label)

  ## at the 1.5 Ma slice, a & b share a clade; c & d share a (different) clade
  ab <- tab$clade_1.5Ma[match(c("a", "b"), tab$species)]
  cd <- tab$clade_1.5Ma[match(c("c", "d"), tab$species)]
  expect_equal(ab[1], ab[2])
  expect_equal(cd[1], cd[2])
  expect_false(ab[1] == cd[1])

  ## at the 0.5 Ma slice every species is in its own clade
  expect_equal(length(unique(tab$clade_0.5Ma)), ape::Ntip(fixed_tree))
})

test_that("build_clade_table validates inputs", {
  expect_error(build_clade_table(list(1), c(1, 2)), "same length")
  expect_error(build_clade_table(list(list(foo = 1)), 1), "extract_clade_sizes")
})

test_that("ancestor_descendant_stats summarizes adjacent slice pairs", {
  ages   <- c(1.5, 0.5)
  slices <- lapply(ages, function(a) extract_clade_sizes(fixed_tree, a))
  tab    <- build_clade_table(slices, ages)
  ad     <- ancestor_descendant_stats(tab)

  expect_true(all(c("old_slice", "new_slice", "parent_clade", "parent_size",
                    "n_desc", "mean_desc_size", "sd_desc_size",
                    "max_desc_size") %in% names(ad)))
  ## two parents at the 1.5 Ma slice (the two cherries), one slice pair
  expect_equal(nrow(ad), 2L)
  ## each cherry has 2 species and splits into 2 singleton descendants
  expect_true(all(ad$parent_size == 2))
  expect_true(all(ad$n_desc == 2))
  expect_true(all(ad$max_desc_size == 1))
  ## parent sizes across the pair account for all species
  expect_equal(sum(ad$parent_size), ape::Ntip(fixed_tree))
})

test_that("ancestor_descendant_stats validates inputs", {
  expect_error(ancestor_descendant_stats(list()), "data frame")
  expect_error(ancestor_descendant_stats(data.frame(x = 1)), "species")
  expect_error(
    ancestor_descendant_stats(data.frame(species = "a", clade_1Ma = "c1")),
    "at least two slice columns"
  )
})
