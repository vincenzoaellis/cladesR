## Phase 3 visualization. The three data builders (make_rect_tree_coords,
## build_truncated_tree_segments, build_collapsed_clade_data) are pure data and
## get exact, deterministic unit tests on a fixed tree. plot_cladetracker only
## gets smoke tests: it must return a ggplot that builds without error (pixel
## comparisons via vdiffr are intentionally avoided as too fragile).

skip_if_not_installed("ape")
skip_if_not_installed("phytools")

## Balanced ultrametric 4-tip tree, total height 3.
##   tips a,b,c,d at x = 3, y = 1,2,3,4
##   cherries (nodes 6,7) at x = 2; root (node 5) at x = 0
fixed_tree <- ape::read.tree(text = "((a:1,b:1):2,(c:1,d:1):2);")

test_that("make_rect_tree_coords returns correct phylogram coordinates", {
  co <- make_rect_tree_coords(fixed_tree)
  expect_equal(co$tree_height, 3)
  expect_equal(nrow(co$nodes), 7L)          # 4 tips + 3 internal
  expect_equal(nrow(co$tips), 4L)
  expect_equal(nrow(co$edges), 6L)          # one per edge
  ## tips sit at the tree height; root sits at x = 0
  expect_true(all(co$tips$x == 3))
  expect_equal(co$nodes$x[co$nodes$node == 5], 0)
  ## tip y-coordinates follow tip.label order a,b,c,d -> 1,2,3,4
  ty <- stats::setNames(co$tips$y, co$tips$species)
  expect_equal(unname(ty[c("a", "b", "c", "d")]), c(1, 2, 3, 4))
  ## edge endpoints reference parent/child node coordinates consistently
  expect_true(all(co$edges$x_child >= co$edges$x_parent))
})

test_that("build_truncated_tree_segments validates the age range", {
  expect_error(build_truncated_tree_segments(fixed_tree, c(1.5, 3)),
               "strictly between 0 and the tree height")
  expect_error(build_truncated_tree_segments(fixed_tree, c(1.5, 0)),
               "strictly between 0 and the tree height")
  expect_silent(build_truncated_tree_segments(fixed_tree, c(1.5, 0.5)))
})

test_that("build_truncated_tree_segments returns one panel per slice", {
  seg <- build_truncated_tree_segments(fixed_tree, c(1.5, 0.5))
  expect_equal(nrow(seg$panels), 2L)
  expect_equal(seg$tree_height, 3)
  expect_true(all(c("horizontal", "vertical", "panels", "tree_height") %in% names(seg)))
  ## panels are laid out left-to-right (oldest first => smaller tree_xmin)
  expect_true(seg$panels$tree_xmin[seg$panels$age == 1.5] <
              seg$panels$tree_xmin[seg$panels$age == 0.5])
})

test_that("build_collapsed_clade_data builds triangles and connectors", {
  ages   <- c(1.5, 0.5)
  slices <- lapply(ages, function(a) extract_clade_sizes(fixed_tree, a))
  cd     <- build_collapsed_clade_data(fixed_tree, slices, ages)

  expect_true(all(c("tree_horizontal", "tree_vertical", "panels",
                    "triangles", "triangle_meta", "connectors") %in% names(cd)))
  ## 2 clades at the 1.5 Ma slice + 4 singletons at 0.5 Ma = 6 triangles
  expect_equal(nrow(cd$triangle_meta), 6L)
  expect_equal(nrow(cd$triangles), 6L * 3L)          # 3 vertices per triangle
  ## each young singleton connects to its parent cherry -> 4 connectors
  expect_equal(nrow(cd$connectors), 4L)
  ## triangle area scales with clade size: the 2-species cherries are wider
  ## than the singletons
  meta <- cd$triangle_meta
  w    <- meta$x_tip - meta$x_base
  expect_true(all(w[meta$n_species == 2] > w[meta$n_species == 1]))
})

test_that("build_collapsed_clade_data validates inputs", {
  ages   <- c(1.5, 0.5)
  slices <- lapply(ages, function(a) extract_clade_sizes(fixed_tree, a))
  expect_error(build_collapsed_clade_data(fixed_tree, slices, 1.5), "same length")
  expect_error(build_collapsed_clade_data(fixed_tree, list(list(foo = 1), slices[[2]]), ages),
               "extract_clade_sizes")
  ## species from a different tree are rejected
  other <- ape::read.tree(text = "((w:1,x:1):2,(y:1,z:1):2);")
  other_slices <- lapply(ages, function(a) extract_clade_sizes(other, a))
  expect_error(build_collapsed_clade_data(fixed_tree, other_slices, ages),
               "same tree")
})

test_that("plot_cladetracker returns a ggplot that builds (smoke tests)", {
  skip_if_not_installed("ggplot2")

  ages   <- c(1.5, 0.5)
  slices <- lapply(ages, function(a) extract_clade_sizes(fixed_tree, a))

  ## default plus each major option combination should build without error
  variants <- list(
    default        = list(),
    color_clades   = list(color_clades = TRUE),
    heat_splits    = list(heat_splits = TRUE),
    heat_scaled    = list(heat_splits = TRUE, scale_by_slice = TRUE),
    flipped        = list(flip_triangles = TRUE),
    labels         = list(show_clade_labels = TRUE, show_split_counts = TRUE),
    no_extras      = list(show_slice_ages = FALSE, show_legend = FALSE)
  )

  for (nm in names(variants)) {
    p <- do.call(plot_cladetracker,
                 c(list(tree = fixed_tree, slice_list = slices, slice_ages = ages),
                   variants[[nm]]))
    expect_s3_class(p, "ggplot")
    expect_silent(ggplot2::ggplot_build(p))
  }
})

test_that("plot_cladetracker warns when both color modes are requested", {
  skip_if_not_installed("ggplot2")
  ages   <- c(1.5, 0.5)
  slices <- lapply(ages, function(a) extract_clade_sizes(fixed_tree, a))
  expect_warning(
    plot_cladetracker(fixed_tree, slices, ages,
                      color_clades = TRUE, heat_splits = TRUE),
    "cannot both be TRUE"
  )
})
