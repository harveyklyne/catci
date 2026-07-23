# Differential-test fixture export (CODE_REVIEW.md step 2).
#
# Freezes the *deterministic seams* of the R method as JSON so the future Python
# port can be tested against them. RNG streams will not match across languages,
# so we deliberately avoid the bootstrap/tie-break paths and pin only the
# input -> output maps that are pure functions of their arguments:
#
#   1. form_T_Sigma   : (x, y, f, g, normalise) -> (T_vector, Sigma)
#   2. search_paths   : (T, Sigma, dx, dy, method) -> criterion values + partitions
#   3. scalar_methods : (T, Sigma) -> mGCM / max / euclid
#   4. approx_chi     : (normsq, tr, tr2) -> Box (1954) chi-square CDF
#   5. rank_one_updates: fast update formulae (24)-(27) vs dense recomputation
#   6. tree_structure : make_binary_tree + permitted-merge maps (the structure object)
#
# Run from the package root:  Rscript data-raw/fixtures/export_fixtures.R
# Regenerate whenever the frozen R method changes (and re-tag).

suppressMessages(devtools::load_all("."))
library(jsonlite)

set.seed(20260723L)

out_dir <- "fixtures"
dir.create(out_dir, showWarnings = FALSE)

git_commit <- tryCatch(
  trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE)),
  error = function(e) NA_character_
)

# ---- serialization helpers ------------------------------------------------
# Ambiguity we must remove for numpy: matrix orientation and length-1 arrays.
# Matrices are stored row-major with explicit dims; every vector that could be
# length 1 is wrapped in I() so jsonlite always emits a JSON array.

mat_to_json <- function(M) {
  M <- as.matrix(M)
  list(
    nrow = jsonlite::unbox(nrow(M)),
    ncol = jsonlite::unbox(ncol(M)),
    data = lapply(seq_len(nrow(M)), function(i) I(as.numeric(M[i, ])))
  )
}

vec_to_json <- function(v) I(as.numeric(v))
ivec_to_json <- function(v) I(as.integer(v))

# A partition is list(x = <groups>, y = <groups>); each group is an integer
# vector of original labels. Force every group to a JSON array.
partition_to_json <- function(categories) {
  grp <- function(gs) lapply(gs, function(g) I(as.integer(unlist(g))))
  list(x = grp(categories$x), y = grp(categories$y))
}

fixtures <- list(
  metadata = list(
    description = jsonlite::unbox("catci differential-test fixtures: deterministic seams of the R method, frozen for the Python port."),
    generated   = jsonlite::unbox(as.character(Sys.Date())),
    git_commit  = jsonlite::unbox(git_commit),
    r_version   = jsonlite::unbox(R.version.string),
    seed        = jsonlite::unbox(20260723L),
    tolerance   = jsonlite::unbox("Compare floats with atol=1e-9, rtol=1e-7. Matrices are row-major with explicit nrow/ncol. All arrays are JSON arrays even when length 1.")
  )
)

# ---------------------------------------------------------------------------
# 1. form_T_Sigma : fixed (x, y, f, g) -> (T_vector, Sigma)
# ---------------------------------------------------------------------------
# Build realistic fixed inputs: row-stochastic propensities f, g and labels
# drawn from them, so Sigma is well-conditioned. The *arrays* are the fixture
# input (frozen into JSON) -- Python reads them verbatim, so the RNG used to
# create them here is irrelevant to the cross-language comparison.

n  <- 200L
dx <- 8L
dy <- 8L

rand_probs <- function(n, d) {
  M <- matrix(stats::runif(n * d, 0.2, 1.0), n, d)
  M / rowSums(M)
}
f <- rand_probs(n, dx)
g <- rand_probs(n, dy)
x <- apply(f, 1, function(p) sample.int(dx, 1, prob = p))
y <- apply(g, 1, function(p) sample.int(dy, 1, prob = p))

ts_norm  <- form_T_Sigma(x, y, f, g, normalise = TRUE)
ts_unorm <- form_T_Sigma(x, y, f, g, normalise = FALSE)

fixtures$form_T_Sigma <- list(
  description = jsonlite::unbox("Pure map (x,y,f,g,normalise) -> (T_vector, Sigma). No RNG in form_T_Sigma itself."),
  inputs = list(
    n = jsonlite::unbox(n), dx = jsonlite::unbox(dx), dy = jsonlite::unbox(dy),
    x = ivec_to_json(x), y = ivec_to_json(y),
    f = mat_to_json(f), g = mat_to_json(g)
  ),
  cases = list(
    list(normalise = jsonlite::unbox(TRUE),
         T_vector = vec_to_json(ts_norm$T_vector),  Sigma = mat_to_json(ts_norm$Sigma)),
    list(normalise = jsonlite::unbox(FALSE),
         T_vector = vec_to_json(ts_unorm$T_vector), Sigma = mat_to_json(ts_unorm$Sigma))
  )
)

# The (T, Sigma) reused by the search / scalar / update fixtures below.
T_vector <- ts_norm$T_vector
Sigma    <- ts_norm$Sigma

# ---------------------------------------------------------------------------
# 2. search_paths : greedy_query criterion values + partition sequence
# ---------------------------------------------------------------------------
# colsample_bylevel = 1 makes `values` and the partition sequence deterministic
# functions of (T, Sigma). We pin ordinal, tree, and greedy. This is the
# fixture that catches finding #1 (tree != ordinal path).

run_search <- function(xsearch, ysearch, trees) {
  res <- greedy_query(T_vector, Sigma, dx = dx, dy = dy,
                      metric = "approx_chi",
                      xsearch = xsearch, ysearch = ysearch,
                      colsample_bylevel = 1, trees = trees)
  list(
    xsearch = jsonlite::unbox(xsearch),
    ysearch = jsonlite::unbox(ysearch),
    values = vec_to_json(res$values),
    partitions = lapply(res$categories, partition_to_json)
  )
}

fixtures$search_paths <- list(
  description = jsonlite::unbox("greedy_query on the shared (T, Sigma) [normalise=TRUE] at dx=dy=8. `values`[0] is the criterion before any merge; each subsequent entry is the max criterion at that level. `partitions`[k] is the label partition after k merges."),
  dx = jsonlite::unbox(dx), dy = jsonlite::unbox(dy),
  T_vector = vec_to_json(T_vector),
  Sigma = mat_to_json(Sigma),
  cases = list(
    ordinal = run_search("ordinal", "ordinal", list(NULL, NULL)),
    tree    = run_search("tree", "tree", list(make_binary_tree(dx), make_binary_tree(dy))),
    greedy  = run_search("greedy", "greedy", list(NULL, NULL))
  )
)

# ---------------------------------------------------------------------------
# 3. scalar_methods : non-adaptive comparators on the same (T, Sigma)
# ---------------------------------------------------------------------------
fixtures$scalar_methods <- list(
  description = jsonlite::unbox("Non-adaptive query_lookup methods (depth-1 searches) on the shared (T, Sigma)."),
  mGCM   = jsonlite::unbox(query_lookup("mGCM")(T_vector, Sigma, dx, dy)),
  max    = jsonlite::unbox(query_lookup("max")(T_vector, Sigma, dx, dy)),
  euclid = jsonlite::unbox(query_lookup("euclid")(T_vector, Sigma, dx, dy))
)

# ---------------------------------------------------------------------------
# 4. approx_chi : Box (1954) chi-square CDF (must match scipy pchisq)
# ---------------------------------------------------------------------------
chi_cases <- list(
  c(normsq = 5.0,   tr = 10.0, tr2 = 30.0),
  c(normsq = 12.5,  tr = 8.0,  tr2 = 16.0),
  c(normsq = 64.0,  tr = 64.0, tr2 = 64.0),
  c(normsq = 100.0, tr = 20.0, tr2 = 55.0),
  c(normsq = sum(T_vector^2), tr = sum(diag(Sigma)), tr2 = sum(Sigma^2))
)
fixtures$approx_chi <- list(
  description = jsonlite::unbox("approx_chi_metric(normsq, tr, tr2) = pchisq(normsq / (tr2/tr), df = tr^2/tr2)."),
  cases = lapply(chi_cases, function(cc) list(
    normsq = jsonlite::unbox(cc[["normsq"]]),
    tr     = jsonlite::unbox(cc[["tr"]]),
    tr2    = jsonlite::unbox(cc[["tr2"]]),
    value  = jsonlite::unbox(approx_chi_metric(cc[["normsq"]], cc[["tr"]], cc[["tr2"]]))
  ))
)

# ---------------------------------------------------------------------------
# 5. rank_one_updates : fast update formulae vs dense recomputation
# ---------------------------------------------------------------------------
# Pins the mathematically load-bearing part (updates (24)-(27)). For each merge
# we store both the fast-update output AND the dense recomputation, so Python
# can assert fast == dense == R.

update_case <- function(dimension, j1, j2, dxc, dyc, Tv, Sig) {
  index1 <- get_index(dimension, j1, dxc, dyc)
  index2 <- get_index(dimension, j2, dxc, dyc)

  new_T    <- update_T(Tv, index1, index2)
  new_Sig  <- update_Sigma(Sig, index1, index2)
  normsq_f <- update_normsq(sum(Tv^2), Tv, index1, index2)
  tr_f     <- update_tr(sum(diag(Sig)), Sig, index1, index2)
  tr2_f    <- update_tr2(sum(Sig^2), Sig, index1, index2)

  # dense oracle
  normsq_d <- sum(new_T^2)
  tr_d     <- sum(diag(new_Sig))
  tr2_d    <- sum(new_Sig^2)

  list(
    dimension = jsonlite::unbox(dimension),
    j1 = jsonlite::unbox(j1), j2 = jsonlite::unbox(j2),
    index1 = I(as.integer(index1)), index2 = I(as.integer(index2)),
    new_T = vec_to_json(new_T),
    new_Sigma = mat_to_json(new_Sig),
    fast   = list(normsq = jsonlite::unbox(normsq_f), tr = jsonlite::unbox(tr_f), tr2 = jsonlite::unbox(tr2_f)),
    dense  = list(normsq = jsonlite::unbox(normsq_d), tr = jsonlite::unbox(tr_d), tr2 = jsonlite::unbox(tr2_d))
  )
}

fixtures$rank_one_updates <- list(
  description = jsonlite::unbox("Merge labels j1,j2 in dimension (1=X, 2=Y) on the shared (T, Sigma). index1/index2 are 0/1 masks (row-order (1,1),(2,1),...,(dx,1),(1,2),...). `fast` uses update formulae (24)-(27); `dense` recomputes from new_T/new_Sigma; they must agree."),
  dx = jsonlite::unbox(dx), dy = jsonlite::unbox(dy),
  cases = list(
    update_case(1L, 1L, 2L, dx, dy, T_vector, Sigma),
    update_case(1L, 3L, 7L, dx, dy, T_vector, Sigma),
    update_case(2L, 2L, 5L, dx, dy, T_vector, Sigma),
    update_case(2L, 1L, 8L, dx, dy, T_vector, Sigma)
  )
)

# ---------------------------------------------------------------------------
# 6. tree_structure : make_binary_tree + permitted-merge maps
# ---------------------------------------------------------------------------
# Pins the structure object the Python redesign centres on (§4.1): the map
# `permitted_merges(structure, partition) -> index pairs`. For each search we
# record, on the initial (unmerged) partition, get_ind2 for every ind1 and the
# total get_num_levels. Also freeze the raw tree shape.

tree_shape <- function(node) {
  # recursively serialize the tree: root labels + children
  list(
    root = I(as.integer(node$root)),
    children = if (length(node$children) == 0) list() else lapply(node$children, tree_shape)
  )
}

permitted_map <- function(d, search, tree) {
  category <- as.list(seq_len(d))
  pairs <- list()
  for (j1 in seq_len(d - 1)) {
    j2s <- get_ind2(ind1 = j1, d = d, search = search, tree = tree, category = category)
    for (j2 in j2s) pairs <- c(pairs, list(I(c(as.integer(j1), as.integer(j2)))))
  }
  list(
    search = jsonlite::unbox(search),
    d = jsonlite::unbox(as.integer(d)),
    num_levels = jsonlite::unbox(get_num_levels(search, d, tree, category)),
    pairs = pairs
  )
}

tree_dims <- c(2L, 4L, 8L)
fixtures$tree_structure <- list(
  description = jsonlite::unbox("make_binary_tree(d) shape, and permitted-merge maps on the INITIAL partition {1},{2},...,{d}. `pairs` are [ind1, ind2] 1-based label indices permitted for merging; `num_levels` is their count (guarded to 0 at d=2)."),
  trees = lapply(tree_dims, function(d) list(d = jsonlite::unbox(d), shape = tree_shape(make_binary_tree(d)))),
  permitted = lapply(tree_dims, function(d) {
    tr <- make_binary_tree(d)
    list(
      d = jsonlite::unbox(d),
      ordinal = permitted_map(d, "ordinal", NULL),
      greedy  = permitted_map(d, "greedy", NULL),
      tree    = permitted_map(d, "tree", tr)
    )
  })
)

# ---------------------------------------------------------------------------
# Write it out
# ---------------------------------------------------------------------------
path <- file.path(out_dir, "catci_fixtures.json")
write_json(fixtures, path, pretty = TRUE, digits = 17, na = "string", null = "null")
cat("Wrote", path, "\n")
cat("Sections:", paste(setdiff(names(fixtures), "metadata"), collapse = ", "), "\n")
