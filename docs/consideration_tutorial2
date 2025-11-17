# Tutorial: Same Spatial Pattern and Proportions, Different Total Cell Count (N)

> **Goal:** Keep both the **spatial pattern** and **cell-type proportions** fixed, while varying the **total cell count (N)**.

> **Purpose:** Demonstrate that even when tissues share identical geometry and composition, CRAWDAD z-scores may still vary slightly due to differences in cell density and the resulting expected neighborhood counts.
---

## Overview

In this tutorial we will:

1. Build a deterministic 4-cluster spatial pattern (used across all scenarios). 
2. Generate simulated tissues with identical cell-type proportions but different total N (e.g., 2 000, 5 000, 10 000 cells).
3. Run CRAWDAD analysis on each dataset.
4. Compare and interpret how increasing or decreasing overall cell density affects z-scores, even when the spatial structure and proportions remain constant.

---

## 1. Setup
<details>
<summary><b>Explanation</b></summary>

We load the required packages and set up basic plotting options.  
`crawdad` performs spatial relationship analysis, and `ggplot2` is used for visualization.
</details>

```r
suppressPackageStartupMessages({
  library(ggplot2)
  library(crawdad)
})
```

## 2. Build deterministic base pattern
<details>
<summary><b>Explanation</b></summary>
We create a fixed 4-cluster layout representing generic tissue organization:
A / B: occupy the same Gaussian “core” clusters around four corner centers.
C: forms a ring around those cores.
D: fills uniform background.
E: fills uniform “far” regions excluding the cluster disks.
These streams are built once and reused so all simulations share the same spatial geometry.
</details>

```r
# ---- helper functions ----
.subseed <- function(base_seed, name) {
  salt <- sum(utf8ToInt(name))
  as.integer((as.numeric(base_seed) + salt * 9973) %% .Machine$integer.max)
}
.within_bounds <- function(x, y, size) x >= 0 & x <= size & y >= 0 & y <= size
.sample_uniform <- function(n, size) data.frame(x = runif(n, 0, size), y = runif(n, 0, size))

# ---- fixed four centers (corners inset by 'frac') ----
.fixed_four_centers <- function(size, frac = 0.25) {
  stopifnot(frac > 0, frac < 0.5)
  xs <- size * c(frac, 1 - frac,  frac,      1 - frac)
  ys <- size * c(frac, frac,      1 - frac,  1 - frac)
  data.frame(x = xs, y = ys)
}

# ---- Gaussian “core” samples around each center ----
.stream_core <- function(n_stream, centers, sigma, size, base_seed, tag) {
  set.seed(.subseed(base_seed, tag))
  per_center <- ceiling(n_stream / nrow(centers))
  xs <- ys <- numeric(0)
  for (i in seq_len(nrow(centers))) {
    need <- per_center; xi <- yi <- numeric(0)
    while (length(xi) < need) {
      batch <- max(need - length(xi), 512)
      bx <- rnorm(batch, centers$x[i], sigma)
      by <- rnorm(batch, centers$y[i], sigma)
      ok <- .within_bounds(bx, by, size)
      xi <- c(xi, bx[ok]); yi <- c(yi, by[ok])
    }
    xs <- c(xs, xi[seq_len(need)]); ys <- c(ys, yi[seq_len(need)])
  }
  data.frame(
    x = xs[seq_len(min(n_stream, length(xs)))],
    y = ys[seq_len(min(n_stream, length(ys)))]
  )
}

# ---- Annular “ring” samples around each center ----
.stream_ring <- function(n_stream, centers, r_mean, r_sd, size, base_seed, tag) {
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
      ok <- .within_bounds(bx, by, size)
      xi <- c(xi, bx[ok]); yi <- c(yi, by[ok])
    }
    xs <- c(xs, xi[seq_len(need)]); ys <- c(ys, yi[seq_len(need)])
  }
  data.frame(
    x = xs[seq_len(min(n_stream, length(xs)))],
    y = ys[seq_len(min(n_stream, length(ys)))]
  )
}

# ---- Uniform background excluding disks around centers ----
.stream_uniform_excluding <- function(n_stream, size, centers, radius, base_seed, tag) {
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

# ---- Build deterministic streams once ----
build_streams <- function(
    N = 5000,
    tissue_size = 1000,
    core_sigma = 35,
    ring_radius = 120,
    ring_sd = 20,
    exclusion_radius = 140,
    base_seed = 42,
    stream_multiplier = 1.2,
    centers_frac = 0.25
) {
  centers  <- .fixed_four_centers(tissue_size, centers_frac)
  n_stream <- ceiling(N * stream_multiplier)
  
  core_stream <- .stream_core(n_stream, centers, core_sigma, tissue_size, base_seed, "core")
  ring_stream <- .stream_ring(n_stream, centers, ring_radius, ring_sd, tissue_size, base_seed, "ring")
  set.seed(.subseed(base_seed, "bg"))
  bg_stream   <- .sample_uniform(n_stream, tissue_size)
  far_stream  <- .stream_uniform_excluding(n_stream, tissue_size, centers, exclusion_radius, base_seed, "far")
  
  # Stable shuffles per stream 
  .shuffle <- function(df, base_seed, tag) {
    set.seed(.subseed(base_seed, paste0("shuffle_", tag)))
    df[sample.int(nrow(df), nrow(df)), , drop = FALSE]
  }
  
  list(
    centers      = centers,
    core_stream  = .shuffle(core_stream, base_seed, "core"),
    ring_stream  = .shuffle(ring_stream, base_seed, "ring"),
    bg_stream    = .shuffle(bg_stream,   base_seed, "bg"),
    far_stream   = .shuffle(far_stream,  base_seed, "far"),
    params = list(
      N=N, tissue_size=tissue_size, core_sigma=core_sigma,
      ring_radius=ring_radius, ring_sd=ring_sd, exclusion_radius=exclusion_radius,
      base_seed=base_seed, stream_multiplier=stream_multiplier, centers_frac=centers_frac
    )
  )
}

```

## 3. Generate simulated datasets (same pattern, same proportions, different N)
<details> <summary><b>Explanation</b></summary>

In this experiment, all simulated tissues share the exact same spatial pattern and the same cell-type proportions
(A = 0.10, B = 0.10, C = 0.20, D = 0.30, E = 0.30).

The only variable we change is the total number of cells (N) - for example, 2 000, 5 000, and 10 000.
By holding composition and structure constant, we can isolate the effect of density on CRAWDAD’s z-scores.

As cell density increases, each cell has more neighbors within the same distance threshold,
slightly shifting the expected counts used in z-score normalization.

This helps demonstrate that minor z-score fluctuations can arise even when biological organization is unchanged.

</details>

```r
.format_props <- function(celltypes, lvls = c("A","B","C","D","E")) {
  ct  <- table(factor(celltypes, levels = lvls))
  N   <- sum(ct)
  pct <- round(100 * as.numeric(ct) / N, 1)
  paste(sprintf("%s=%s%%", lvls, format(pct, trim = TRUE)), collapse = "   ")
}

simulate_cells <- function(
  proportions = c(A=.20, B=.20, C=.15, D=.25, E=.20),
  N = 5000,
  streams,
  plot = TRUE
) {
  wanted <- c("A","B","C","D","E")
  if (is.null(names(proportions))) names(proportions) <- wanted
  proportions <- proportions[wanted] / sum(proportions)

  counts <- as.integer(round(proportions * N))
  diffN <- N - sum(counts)
  if (diffN != 0) {
    fr  <- proportions * N - round(proportions * N)
    ord <- order(fr, decreasing = TRUE)
    for (i in seq_len(abs(diffN))) {
      idx <- ord[((i - 1) %% length(ord)) + 1]
      counts[idx] <- counts[idx] + sign(diffN)
    }
  }
  names(counts) <- wanted
  nA <- counts["A"]; nB <- counts["B"]; nC <- counts["C"]; nD <- counts["D"]; nE <- counts["E"]

  A_core <- if (nA>0) { out <- streams$core_stream[seq_len(nA), , drop=FALSE]; out$celltype <- "A"; out } else NULL
  B_core <- if (nB>0) { out <- streams$core_stream[(nA+1):(nA+nB), , drop=FALSE]; out$celltype <- "B"; out } else NULL
  C_ring <- if (nC>0) { out <- streams$ring_stream[seq_len(nC), , drop=FALSE]; out$celltype <- "C"; out } else NULL
  D_bg   <- if (nD>0) { out <- streams$bg_stream [seq_len(nD), , drop=FALSE]; out$celltype <- "D"; out } else NULL
  E_far  <- if (nE>0) { out <- streams$far_stream[seq_len(nE), , drop=FALSE]; out$celltype <- "E"; out } else NULL

  cells <- do.call(rbind, Filter(Negate(is.null), list(A_core,B_core,C_ring,D_bg,E_far)))
  cells$celltype <- factor(cells$celltype, levels = wanted)

  if (plot) {
    ts <- streams$params$tissue_size
    subtitle <- paste0("N=", nrow(cells), " | ", .format_props(cells$celltype))
    ggplot(cells, aes(x, y, color = celltype)) +
      geom_point(size = 0.6, alpha = 0.8) +
      coord_fixed(xlim = c(0, ts), ylim = c(0, ts), expand = FALSE) +
      labs(title = "Fixed 4-cluster spatial pattern",
           subtitle = subtitle, color = "Type") +
      theme_minimal(base_size = 12) +
      theme(panel.grid = element_blank(),
            plot.subtitle = element_text(lineheight = 1.1))
  } else {
    invisible(cells)
  }
}

baseline_props <- c(A = 0.10, B = 0.10, C = 0.20, D = 0.30, E = 0.30)
Ns <- c(2000, 5000, 10000)

# Build streams ONCE at max N so spatial geometry is identical across densities
streams_fixed <- build_streams(
  N = max(Ns),
  tissue_size = 1000,
  core_sigma = 35,
  ring_radius = 120,
  ring_sd = 20,
  exclusion_radius = 140,
  base_seed = 42,
  stream_multiplier = 1.2,
  centers_frac = 0.25
)

# Helper to slice the same streams to a given N with the same proportions
.simulate_fixed_props_varyN <- function(N) {
  sim <- simulate_cells(
    proportions = baseline_props,
    N = N,
    streams = streams_fixed,
    plot = TRUE
  )
  sim$name <- paste0("N", N)
  sim
}

sims <- setNames(lapply(Ns, .simulate_fixed_props_varyN), paste0("N", Ns))

```

<table align="center">
  <thead>
    <tr>
      <th align="center">N=2000</th>
      <th align="center">N=5000</th>
      <th align="center">N=10000</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td align="center"><img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/considerations_tutorials_files/tutorial2/simulation_N2000.png" width="380"></td>
      <td align="center"><img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/considerations_tutorials_files/tutorial2/simulation_N5000.png" width="380"></td>
      <td align="center"><img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/considerations_tutorials_files/tutorial2/simulation_N10000.png" width="380"></td>
    </tr>
  </tbody>
</table>

<p align="center"><b>Figure 1.</b> Simulation patterns showing identical spatial layout and cell-type proportions but different total numbers of cells.</p>

## 4. CRAWDAD analysis
<details> <summary><b>Explanation</b></summary>
Now we analyze each simulated tissue.
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
      <td><b>N=2000</b></td>
      <td align="center"><img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/considerations_tutorials_files/tutorial2/simulation_N2000.png" width="400"></td>
      <td align="center"><img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/considerations_tutorials_files/tutorial2/coloc_dotplot_N2000.png" width="400"></td>
    </tr>
    <tr>
      <td><b>N=5000</b></td>
      <td align="center"><img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/considerations_tutorials_files/tutorial2/simulation_N5000.png" width="400"></td>
      <td align="center"><img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/considerations_tutorials_files/tutorial2/coloc_dotplot_N5000.png" width="400"></td>
    </tr>
        <tr>
      <td><b>N=10000</b></td>
      <td align="center"><img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/considerations_tutorials_files/tutorial2/simulation_N10000.png" width="400"></td>
      <td align="center"><img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/considerations_tutorials_files/tutorial2/coloc_dotplot_N10000.png" width="400"></td>
    </tr>
  </tbody>
</table>


<p align="center"><b>Figure 2.</b> Simulation patterns showing identical spatial layout and the same cell type proportions but different total number of cells.</p>


## 5. How to interpret the results
<details> <summary><b>Explanation</b></summary>
Each dot-plot shows pairwise spatial associations between cell types.
Rows represent the reference cell type, and columns represent its neighbor.
Red indicates co-localization (positive z-scores), while blue indicates spatial depletion (negative z-scores).


When comparing the dot-plots across simulations, we see that both the spatial pattern and cell-type proportions remain perfectly constant.
The only factor that changes is the total number of cells (N), which affects the overall cell density in the tissue.

In our simulations, the direction of the relationships (which cell types tend to co-localize or avoid each other) stays the same, but the magnitude of the z-scores varies slightly.
This happens because CRAWDAD’s z-score depends on expected neighborhood counts - when density increases, more cells fall within a fixed distance, subtly altering the expected and observed values used in z-score normalization.
These differences reflect statistical scaling, not biological change.

To make results more comparable across tissues with different densities, one can normalize by expected neighborhood counts or apply density matching (for example, by downsampling) before analysis.
This ensures that z-score differences reflect true biological organization, not sampling density.
</details>
