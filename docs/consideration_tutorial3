# Tutorial: Same Total Cell Count and Proportions, Different Spatial Patterns

> **Goal:** Keep **cell-type proportions** and **total cell count (N)** constant, while varying **the spatial arrangement of cell types**.

> **Purpose:** Demonstrate that CRAWDAD z-scores change meaningfully when **the spatial organization** of cell types differs, if the overall composition and total number of cells remain constant, reflecting true biological differences in spatial relationships rather than sampling effects.
---

## Overview

In this tutorial we will:
1. Define several distinct spatial geometries (e.g., four clusters, center, diagonal).
2. Assign cell types A–E to different spatial roles (core, ring, background, far) in each geometry.
3. Keep the total cell count (N = 5000) and cell-type proportions constant across all simulations.
4. Run CRAWDAD analysis on each dataset.
5. Compare and interpret how changing spatial patterns affects the resulting z-scores, while keeping total cell count and proportions constant.

---

## 1. Setup
<details>
<summary><b>Explanation</b></summary>
We load the required packages and define basic sampling functions for building different spatial arrangements.
These functions generate points for each cell type under distinct patterns (e.g., four clusters, central, diagonal).
Each pattern is generated deterministically to make the results reproducible.
</details>

```r
suppressPackageStartupMessages({
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install ggplot2")
})

.subseed <- function(base_seed, name) {
  salt <- sum(utf8ToInt(name))
  as.integer((as.numeric(base_seed) + salt * 9973) %% .Machine$integer.max)
}
within_bounds <- function(x, y, size) x >= 0 & x <= size & y >= 0 & y <= size
sample_uniform <- function(n, size) data.frame(x = runif(n, 0, size), y = runif(n, 0, size))


# ---- stream generators (core/ring/bg/far) ----
stream_gaussian_around_centers <- function(n_stream, centers, sigma, size, base_seed, tag) {
  set.seed(.subseed(base_seed, tag))
  per_center <- ceiling(n_stream / nrow(centers))
  xs <- ys <- numeric(0)
  for (i in seq_len(nrow(centers))) {
    need <- per_center; xi <- yi <- numeric(0)
    while (length(xi) < need) {
      batch <- max(need - length(xi), 512)
      bx <- rnorm(batch, centers$x[i], sigma)
      by <- rnorm(batch, centers$y[i], sigma)
      ok <- within_bounds(bx, by, size)
      xi <- c(xi, bx[ok]); yi <- c(yi, by[ok])
    }
    xs <- c(xs, xi[seq_len(need)]); ys <- c(ys, yi[seq_len(need)])
  }
  data.frame(x = xs[seq_len(min(n_stream, length(xs)))],
             y = ys[seq_len(min(n_stream, length(ys)))])
}

stream_ring_around_centers <- function(n_stream, centers, r_mean, r_sd, size, base_seed, tag) {
  set.seed(.subseed(base_seed, tag))
  per_center <- ceiling(n_stream / nrow(centers))
  xs <- ys <- numeric(0)
  for (i in seq_len(nrow(centers))) {
    need <- per_center; xi <- yi <- numeric(0)
    while (length(xi) < need) {
      batch <- max(need - length(xi), 512)
      theta <- runif(batch, 0, 2*pi)
      r <- pmax(0, rnorm(batch, r_mean, r_sd))
      bx <- centers$x[i] + r * cos(theta)
      by <- centers$y[i] + r * sin(theta)
      ok <- within_bounds(bx, by, size)
      xi <- c(xi, bx[ok]); yi <- c(yi, by[ok])
    }
    xs <- c(xs, xi[seq_len(need)]); ys <- c(ys, yi[seq_len(need)])
  }
  data.frame(x = xs[seq_len(min(n_stream, length(xs)))],
             y = ys[seq_len(min(n_stream, length(ys)))])
}

stream_uniform_excluding_centers <- function(n_stream, size, centers, radius, base_seed, tag) {
  set.seed(.subseed(base_seed, tag))
  xs <- ys <- numeric(0); r2 <- radius^2
  while (length(xs) < n_stream) {
    batch <- max(n_stream - length(xs), 2048)
    bx <- runif(batch, 0, size); by <- runif(batch, 0, size)
    keep <- rep(TRUE, batch)
    for (i in seq_len(nrow(centers))) {
      dx <- bx - centers$x[i]; dy <- by - centers$y[i]
      keep <- keep & (dx*dx + dy*dy > r2)
    }
    if (any(keep)) { xs <- c(xs, bx[keep]); ys <- c(ys, by[keep]) }
  }
  data.frame(x = xs[seq_len(n_stream)], y = ys[seq_len(n_stream)])
}

# ---- build shuffled streams for a chosen geometry ----
build_streams_by_pattern <- function(
  pattern = c("four_corners","center","two_diagonal","line_horizontal","grid3x3"),
  N = 6000, tissue_size = 1000, core_sigma = 35,
  ring_radius = 120, ring_sd = 20, exclusion_radius = 140,
  base_seed = 42, stream_multiplier = 1.2, centers_frac = 0.25
){
  pattern <- match.arg(pattern)
  centers <- switch(
    pattern,
    four_corners    = fixed_four_centers(tissue_size, centers_frac),
    center          = center_middle(tissue_size),
    two_diagonal    = centers_two_diagonal(tissue_size, frac = centers_frac),
    line_horizontal = centers_line_horizontal(tissue_size, xs = c(0.2,0.4,0.6,0.8)),
    grid3x3         = centers_grid(tissue_size, kx = 3, ky = 3, margin_frac = 0.18)
  )
  n_stream <- ceiling(N * stream_multiplier)

  core_stream <- stream_gaussian_around_centers(
    n_stream, centers, core_sigma, tissue_size, base_seed, paste0("core_", pattern)
  )
  ring_stream <- stream_ring_around_centers(
    n_stream, centers, ring_radius, ring_sd, tissue_size, base_seed, paste0("ring_", pattern)
  )
  set.seed(.subseed(base_seed, paste0("bg_", pattern)))
  bg_stream <- sample_uniform(n_stream, tissue_size)
  far_stream <- stream_uniform_excluding_centers(
    n_stream, tissue_size, centers, exclusion_radius, base_seed, paste0("far_", pattern)
  )

  shuffle_idx <- function(n, base_seed, tag) {
    set.seed(.subseed(base_seed, paste0("shuffle_", tag))); sample.int(n, n)
  }
  core_stream <- core_stream[shuffle_idx(nrow(core_stream), base_seed, paste0("core_", pattern)), , drop = FALSE]
  ring_stream <- ring_stream[shuffle_idx(nrow(ring_stream), base_seed, paste0("ring_", pattern)), , drop = FALSE]
  bg_stream   <- bg_stream  [shuffle_idx(nrow(bg_stream),   base_seed, paste0("bg_", pattern)),   , drop = FALSE]
  far_stream  <- far_stream [shuffle_idx(nrow(far_stream),  base_seed, paste0("far_", pattern)),  , drop = FALSE]

  list(
    pattern = pattern, centers = centers,
    core_stream = core_stream, ring_stream = ring_stream,
    bg_stream = bg_stream, far_stream = far_stream,
    params = list(
      N = N, tissue_size = tissue_size, core_sigma = core_sigma,
      ring_radius = ring_radius, ring_sd = ring_sd,
      exclusion_radius = exclusion_radius, base_seed = base_seed,
      stream_multiplier = stream_multiplier, centers_frac = centers_frac
    )
  )
}

# ---- role-based slicing with fixed proportions ----
simulate_cells_roles <- function(
  proportions, N, streams, role_map, plot = TRUE
){
  wanted <- c("A","B","C","D","E")
  if (is.null(names(proportions))) names(proportions) <- wanted
  proportions <- proportions[wanted] / sum(proportions)

  counts <- as.integer(round(proportions * N))
  diffN  <- N - sum(counts)
  if (diffN != 0) {
    fr <- proportions * N - round(proportions * N)
    ord <- order(fr, decreasing = TRUE)
    for (i in seq_len(abs(diffN))) {
      idx <- ord[((i - 1) %% length(ord)) + 1]
      counts[idx] <- counts[idx] + sign(diffN)
    }
  }
  names(counts) <- wanted

  groups <- c("core","ring","bg","far")
  stopifnot(all(groups %in% names(role_map)))
  all_labels <- unlist(role_map[groups], use.names = FALSE)
  stopifnot(setequal(all_labels, wanted))  # each label exactly once

  take_from_stream <- function(df_stream, labels_vec) {
    out <- list(); start <- 1L
    for (lab in labels_vec) {
      n <- counts[[lab]]
      if (n > 0) {
        stopifnot(start + n - 1 <= nrow(df_stream))
        chunk <- df_stream[seq.int(start, start + n - 1), , drop = FALSE]
        chunk$celltype <- lab
        out[[lab]] <- chunk
        start <- start + n
      }
    }
    do.call(rbind, out)
  }

  parts <- list(
    core = take_from_stream(streams$core_stream, role_map$core),
    ring = take_from_stream(streams$ring_stream, role_map$ring),
    bg   = take_from_stream(streams$bg_stream,   role_map$bg),
    far  = take_from_stream(streams$far_stream,  role_map$far)
  )
  cells <- do.call(rbind, parts)
  cells$celltype <- factor(cells$celltype, levels = wanted)

  plt <- NULL
  if (plot && requireNamespace("ggplot2", quietly = TRUE)) {
    ts <- streams$params$tissue_size
    subtitle_line1 <- sprintf(
      "Pattern=%s | Roles: core=%s; ring=%s; bg=%s; far=%s",
      streams$pattern,
      paste(role_map$core, collapse = ","),
      paste(role_map$ring, collapse = ","),
      paste(role_map$bg,   collapse = ","),
      paste(role_map$far,  collapse = ",")
    )
    subtitle_line2 <- paste0("N=", nrow(cells), "   |   ", .format_props(cells))
    plt <- ggplot2::ggplot(cells, ggplot2::aes(x, y, color = celltype)) +
      ggplot2::geom_point(size = 0.6, alpha = 0.8) +
      ggplot2::coord_fixed(xlim = c(0, ts), ylim = c(0, ts), expand = FALSE) +
      ggplot2::labs(
        title = "Fixed proportions, different relationships",
        subtitle = paste(subtitle_line1, subtitle_line2, sep = "\n"),
        color = "type"
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        panel.grid       = ggplot2::element_blank(),
        plot.subtitle    = ggplot2::element_text(lineheight = 1.1),
        plot.background  = ggplot2::element_rect(fill = "white", colour = NA),
        panel.background = ggplot2::element_rect(fill = "white", colour = NA)
      )
  }

  list(
    cells = cells,
    plot = plt,
    pattern = streams$pattern,
    role_map = role_map
  )
}

# ---- suite wrapper returning items compatible with run_crawdad_all ----
simulate_suite_roles <- function(
  scenarios, proportions,
  N = 6000, tissue_size = 1000,
  core_sigma = 35, ring_radius = 120, ring_sd = 20,
  exclusion_radius = 140, base_seed = 42,
  stream_multiplier = 1.2, centers_frac = 0.25,
  plot_each = TRUE
){
  out <- lapply(names(scenarios), function(nm) {
    sc  <- scenarios[[nm]]
    str <- build_streams_by_pattern(
      pattern = sc$pattern, N = N, tissue_size = tissue_size,
      core_sigma = core_sigma, ring_radius = ring_radius, ring_sd = ring_sd,
      exclusion_radius = exclusion_radius, base_seed = base_seed,
      stream_multiplier = stream_multiplier, centers_frac = centers_frac
    )
    sim <- simulate_cells_roles(
      proportions = proportions, N = N, streams = str,
      role_map = sc$role_map, plot = plot_each
    )
    list(
      name = nm,
      data = list(data = sim$cells),
      plot = sim$plot,
      pattern = sim$pattern,
      role_map = sim$role_map
    )
  })
  names(out) <- names(scenarios)
  out
}

```

## 2. Define geometric patterns
<details>
<summary><b>Explanation</b></summary>
We define several alternative geometric layouts representing different tissue architectures:
four_corners - Four clusters placed at each corner.
center - One large cluster in the middle.
two_diagonal - Two clusters along the diagonal.
Each geometry is later combined with a role map defining which cell types occupy core, ring, background, and far regions.
</details>

```r
fixed_four_centers <- function(size, frac = 0.25) {
  xs <- size * c(frac, 1 - frac, frac, 1 - frac)
  ys <- size * c(frac, frac, 1 - frac, 1 - frac)
  data.frame(x = xs, y = ys)
}
center_middle <- function(size) data.frame(x = size/2, y = size/2)
centers_two_diagonal <- function(size, frac = 0.25) {
  data.frame(x = size * c(frac, 1 - frac), y = size * c(frac, 1 - frac))
}
centers_line_horizontal <- function(size, xs = c(0.2, 0.4, 0.6, 0.8)) {
  data.frame(x = size * xs, y = rep(size * 0.5, length(xs)))
}
centers_grid <- function(size, kx = 3, ky = 3, margin_frac = 0.18) {
  xs <- seq(margin_frac, 1 - margin_frac, length.out = kx) * size
  ys <- seq(margin_frac, 1 - margin_frac, length.out = ky) * size
  expand.grid(x = xs, y = ys)
}


```

## 3. Generate simulated datasets (same proportions, different spatial roles)
<details> <summary><b>Explanation</b></summary>

Here, we use fixed cell-type proportions (A=0.05, B=0.10, C=0.20, D=0.30, E=0.35) and fixed total N = 5000,
but vary the spatial role assignment and geometry.

Each scenario reuses the same total number of cells and proportions,
but cell types take on different spatial positions (e.g., A and B as cores in one pattern, D as core in another).

This allows us to observe how CRAWDAD responds to genuine spatial rearrangements.
</details>

```r
# Helper for summarizing proportions
.format_props <- function(cells) {
  lvls <- c("A","B","C","D","E")
  cnt  <- table(factor(cells$celltype, levels = lvls))
  N    <- sum(cnt)
  pct  <- round(100 * as.numeric(cnt) / N, 1)
  paste(sprintf("%s=%s%%", lvls, format(pct, trim = TRUE)), collapse = "   ")
}

# Define fixed proportions and scenarios
fixed_props <- c(A=0.05, B=0.10, C=0.20, D=0.30, E=0.35)

scenarios <- list(
  S1 = list(
    pattern  = "four_corners",
    role_map = list(core = c("A","B"), ring = c("C"), bg = c("E"), far = c("D"))
  ),
  S2 = list(
    pattern  = "center",
    role_map = list(core = c("D"), ring = c("E"), bg = c("B","C"), far = c("A"))
  ),
  S3 = list(
    pattern  = "two_diagonal",
    role_map = list(core = c("C","E"), ring = c("A"), bg = c("B"), far = c("D"))
  )
)

# Run the suite
sims <- simulate_suite_roles(
  scenarios     = scenarios,
  proportions   = fixed_props,
  N             = 5000,
  tissue_size   = 1000,
  base_seed     = 42,
  plot_each     = TRUE
)

```

<table align="center">
  <thead>
    <tr>
      <th align="center">S1</th>
      <th align="center">S2</th>
      <th align="center">S3</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td align="center"><img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/considerations_tutorials_files/tutorial3/coloc_dotplot_S1.png" width="380"></td>
      <td align="center"><img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/considerations_tutorials_files/tutorial3/coloc_dotplot_S2.png" width="380"></td>
      <td align="center"><img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/considerations_tutorials_files/tutorial3/coloc_dotplot_S3.png" width="380"></td>
    </tr>
  </tbody>
</table>

<p align="center"><b>Figure 1.</b> Simulated tissues with identical cell-type proportions and total N = 5000, but distinct spatial organizations and role assignments.</p>

## 4. CRAWDAD analysis
<details> <summary><b>Explanation</b></summary>
We analyze each simulated tissue with the same CRAWDAD parameters
(neighDist = 30, scales = 30–500).
All shuffles and distances are identical across simulations, ensuring differences in z-scores reflect only spatial patterns.
</details>

```r
run_crawdad_all <- function(
  sims,
  neighDist = 30,
  scales = seq(30, 500, by = 30),
  perms = 3,
  seed = 123,
  verbose = TRUE
) {
  if (min(scales) < neighDist)
    scales <- sort(unique(c(neighDist, scales[scales >= neighDist])))

  lapply(sims, function(sim) {
    nm <- sim$name; df <- sim$data$data
    message(">>> CRAWDAD for ", nm)
    cells_sf <- crawdad::toSF(pos = df[, c("x","y")], cellTypes = df$celltype)
    shuffleList <- crawdad:::makeShuffledCells(
      cells_sf, scales = scales, perms = perms,
      seed = seed, verbose = verbose
    )
    results <- crawdad::findTrends(
      cells_sf, neighDist = neighDist, shuffleList = shuffleList,
      verbose = verbose, returnMeans = FALSE
    )
    dat  <- crawdad::meltResultsList(results, withPerms = TRUE)
    zthr <- crawdad::correctZBonferroni(dat)

    p <- crawdad::vizColocDotplot(
      dat, zSigThresh = zthr, zScoreLimit = 6, symmetrical = TRUE,
      dotSizes = c(3, 12)
    ) +
      ggplot2::theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
      ggplot2::labs(title = paste0("CRAWDAD analysis — ", nm)) +
      ggplot2::scale_colour_gradient2(
        low = "blue", mid = "white", high = "red",
        midpoint = 0, limits = c(-6, 6),
        oob = scales::squish, na.value = "lightgray"
      )

    list(name = nm, dat = dat, zsig = zthr, plot = p)
  })
}

crawdad_out <- run_crawdad_all(sims)

```

<table>
  <thead>
    <tr>
      <th align="left">Scenario</th>
      <th align="center">Visualization</th>
      <th align="center">Dot-plot</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><b>S1</b></td>
      <td align="center"><img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/considerations_tutorials_files/tutorial3/sim_plot_S1.png" width="400"></td>
      <td align="center"><img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/considerations_tutorials_files/tutorial3/coloc_dotplot_S1.png" width="400"></td>
    </tr>
    <tr>
      <td><b>S2</b></td>
      <td align="center"><img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/considerations_tutorials_files/tutorial3/sim_plot_S2.png" width="400"></td>
      <td align="center"><img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/considerations_tutorials_files/tutorial3/coloc_dotplot_S2.png" width="400"></td>
    </tr>
    <tr>
      <td><b>S3</b></td>
      <td align="center"><img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/considerations_tutorials_files/tutorial3/sim_plot_S3.png" width="400"></td>
      <td align="center"><img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/considerations_tutorials_files/tutorial3/coloc_dotplot_S3.png" width="400"></td>
    </tr>
  </tbody>
</table>

<p align="center"><b>Figure 2.</b> Simulated tissues with identical cell-type proportions and total N = 5000, but distinct spatial organizations and role assignments.</p>

## 5. How to interpret the results
<details> <summary><b>Explanation</b></summary>
Each dot-plot shows pairwise spatial associations between cell types.
Rows represent the reference cell type, and columns represent its neighbor.
Red indicates co-localization (positive z-scores), while blue indicates spatial depletion (negative z-scores).

In this experiment, both the total number of cells (N) and the cell-type proportions remain constant across simulations.

However, the spatial organization of those cell types is deliberately altered - each scenario assigns A–E to different spatial roles (core, ring, background, far) and uses distinct geometric patterns (e.g., four clusters, center, diagonal).
As a result, the nature of the spatial relationships between cell types changes completely.

Cell pairs that were co-localized in one pattern may become separated or even depleted in another, leading to qualitatively different z-score structures.

Unlike previous tutorials, where only the intensity of associations changed, here the direction and identity of associations change as well - reflecting true biological differences in spatial arrangement rather than statistical artifacts.

These simulations highlight that CRAWDAD captures biologically meaningful spatial variation.When z-score patterns differ across tissues with the same composition and density, it signals genuine differences in cell-type organization rather than sampling bias.
</details>