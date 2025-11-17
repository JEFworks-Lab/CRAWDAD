# Tutorial: Same Spatial Pattern, Same Total Cell Count, Different Cell-Type Proportions

> **Goal:** Hold the **spatial pattern** and **total cell count (N)** fixed, but vary **cell-type proportions**.  

> **Purpose:** Show that CRAWDAD z-scores can differ even when the spatial arrangement is identical, because overall cell composition changes expected neighborhood frequencies.

---

## Overview

In this tutorial we will:

1. Build a deterministic 4-cluster spatial pattern (used across all scenarios).  
2. Generate multiple simulated tissues with **different cell-type proportions** but the same pattern and same total N.  
3. Run CRAWDAD analysis for each simulation at fixed parameters (`neighDist`, `scales`).  
4. Visualize and compare results.

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

## 3. Generate simulated datasets (same pattern, different proportions)
<details> <summary><b>Explanation</b></summary>
Each scenario will have the same N = 5000 cells and identical spatial pattern,
but different proportions for cell types A–E.
This lets us isolate the effect of composition on CRAWDAD results.
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

  A_core <- if (nA>0) { out <- streams$core_stream[seq_len(nA), ]; out$celltype <- "A"; out } else NULL
  B_core <- if (nB>0) { out <- streams$core_stream[(nA+1):(nA+nB), ]; out$celltype <- "B"; out } else NULL
  C_ring <- if (nC>0) { out <- streams$ring_stream[seq_len(nC), ]; out$celltype <- "C"; out } else NULL
  D_bg   <- if (nD>0) { out <- streams$bg_stream[seq_len(nD), ]; out$celltype <- "D"; out } else NULL
  E_far  <- if (nE>0) { out <- streams$far_stream[seq_len(nE), ]; out$celltype <- "E"; out } else NULL

  cells <- do.call(rbind, Filter(Negate(is.null), list(A_core,B_core,C_ring,D_bg,E_far)))
  cells$celltype <- factor(cells$celltype, levels = wanted)

  if (plot) {
    ts <- 1000
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

simulate_suite <- function(scenarios, N = 5000, tissue_size = 1000, base_seed = 42) {
  streams <- build_streams(N = N, tissue_size = tissue_size, base_seed = base_seed)
  lapply(names(scenarios), function(nm) {
    list(name = nm, data = simulate_cells(scenarios[[nm]], N, streams))
  })
}

scenarios <- list(
  S0  = c(A=0.10, B=0.10, C=0.20, D=0.30, E=0.30),
  S1  = c(A=0.25, B=0.05, C=0.20, D=0.30, E=0.20),
  S2  = c(A=0.25, B=0.10, C=0.05, D=0.30, E=0.30),
  S3  = c(A=0.25, B=0.10, C=0.20, D=0.15, E=0.30),
  S4  = c(A=0.25, B=0.10, C=0.20, D=0.30, E=0.15),
  S5  = c(A=0.10, B=0.25, C=0.05, D=0.30, E=0.30),
  S6  = c(A=0.10, B=0.25, C=0.20, D=0.15, E=0.30),
  S7  = c(A=0.10, B=0.25, C=0.20, D=0.30, E=0.15),
  S8  = c(A=0.10, B=0.10, C=0.35, D=0.15, E=0.30),
  S9  = c(A=0.10, B=0.10, C=0.35, D=0.30, E=0.15),
  S10 = c(A=0.10, B=0.10, C=0.20, D=0.45, E=0.15)
)

sims <- simulate_suite(scenarios, N = 5000)

```

<table align="center">
  <thead>
    <tr>
      <th align="center">S0</th>
      <th align="center">S1</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td align="center"><img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/considerations_tutorials_files/tutorial1/plot_S0.png" width="380"></td>
      <td align="center"><img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/considerations_tutorials_files/tutorial1/plot_S1.png" width="380"></td>
    </tr>
    <tr>
      <th align="center">S2</th>
      <th align="center">S3</th>
    </tr>
    <tr>
      <td align="center"><img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/considerations_tutorials_files/tutorial1/plot_S2.png" width="380"></td>
      <td align="center"><img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/considerations_tutorials_files/tutorial1/plot_S3.png" width="380"></td>
    </tr>
    <tr>
      <th align="center">S4</th>
      <th align="center">S5</th>
    </tr>
    <tr>
      <td align="center"><img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/considerations_tutorials_files/tutorial1/plot_S4.png" width="380"></td>
      <td align="center"><img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/considerations_tutorials_files/tutorial1/plot_S5.png" width="380"></td>
    </tr>
    <tr>
      <th align="center">S6</th>
      <th align="center">S7</th>
    </tr>
    <tr>
      <td align="center"><img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/considerations_tutorials_files/tutorial1/plot_S6.png" width="380"></td>
      <td align="center"><img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/considerations_tutorials_files/tutorial1/plot_S7.png" width="380"></td>
    </tr>
    <tr>
      <th align="center">S8</th>
      <th align="center">S9</th>
    </tr>
    <tr>
      <td align="center"><img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/considerations_tutorials_files/tutorial1/plot_S8.png" width="380"></td>
      <td align="center"><img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/considerations_tutorials_files/tutorial1/plot_S9.png" width="380"></td>
    </tr>
  </tbody>
</table>

<p align="center"><b>Figure 1.</b> Simulation patterns showing identical spatial layout but different cell-type proportions.</p>

## 4. CRAWDAD analysis
<details> <summary><b>Explanation</b></summary>
We analyze each simulated tissue with the same CRAWDAD parameters
(neighDist = 30, scales = 30–500).
All shuffles and distances are identical across simulations, ensuring differences in z-scores reflect only composition.
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
      <th align="center">Histogram</th>
      <th align="center">Dot-plot</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><b>S0</b></td>
      <td align="center"><img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/considerations_tutorials_files/tutorial1/plot_S0.png" width="300"></td>
      <td align="center"><img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/considerations_tutorials_files/tutorial1/hist_S0.png" width="300"></td>
      <td align="center"><img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/considerations_tutorials_files/tutorial1/coloc_dotplot_S0.png" width="300"></td>
    </tr>
    <tr>
      <td><b>S1</b></td>
      <td align="center"><img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/considerations_tutorials_files/tutorial1/plot_S1.png" width="300"></td>
      <td align="center"><img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/considerations_tutorials_files/tutorial1/hist_S1.png" width="300"></td>
      <td align="center"><img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/considerations_tutorials_files/tutorial1/coloc_dotplot_S1.png" width="300"></td>
    </tr>
    <tr>
      <td><b>S2</b></td>
      <td align="center"><img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/considerations_tutorials_files/tutorial1/plot_S2.png" width="300"></td>
      <td align="center"><img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/considerations_tutorials_files/tutorial1/hist_S2.png" width="300"></td>
      <td align="center"><img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/considerations_tutorials_files/tutorial1/coloc_dotplot_S2.png" width="300"></td>
    </tr>
    <tr>
      <td><b>S3</b></td>
      <td align="center"><img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/considerations_tutorials_files/tutorial1/plot_S3.png" width="300"></td>
      <td align="center"><img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/considerations_tutorials_files/tutorial1/hist_S3.png" width="300"></td>
      <td align="center"><img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/considerations_tutorials_filess/tutorial1/coloc_dotplot_S3.png" width="300"></td>
    </tr>
    <tr>
      <td><b>S4</b></td>
      <td align="center"><img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/considerations_tutorials_files/tutorial1/plot_S4.png" width="300"></td>
      <td align="center"><img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/considerations_tutorials_files/tutorial1/hist_S4.png" width="300"></td>
      <td align="center"><img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/considerations_tutorials_files/tutorial1/coloc_dotplot_S4.png" width="300"></td>
    </tr>
    <tr>
      <td><b>S5</b></td>
      <td align="center"><img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/considerations_tutorials_files/tutorial1/plot_S5.png" width="300"></td>
      <td align="center"><img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/considerations_tutorials_filess/tutorial1/hist_S5.png" width="300"></td>
      <td align="center"><img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/considerations_tutorials_files/tutorial1/coloc_dotplot_S5.png" width="300"></td>
    </tr>
    <tr>
      <td><b>S6</b></td>
      <td align="center"><img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/considerations_tutorials_files/tutorial1/plot_S6.png" width="300"></td>
      <td align="center"><img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/considerations_tutorials_files/tutorial1/hist_S6.png" width="300"></td>
      <td align="center"><img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/considerations_tutorials_files/tutorial1/coloc_dotplot_S6.png" width="300"></td>
    </tr>
    <tr>
      <td><b>S7</b></td>
      <td align="center"><img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/considerations_tutorials_files/tutorial1/plot_S7.png" width="300"></td>
      <td align="center"><img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/considerations_tutorials_files/tutorial1/hist_S7.png" width="300"></td>
      <td align="center"><img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/considerations_tutorials_files/tutorial1/coloc_dotplot_S7.png" width="300"></td>
    </tr>
    <tr>
      <td><b>S8</b></td>
      <td align="center"><img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/considerations_tutorials_files/tutorial1/plot_S8.png" width="300"></td>
      <td align="center"><img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/considerations_tutorials_files/tutorial1/hist_S8.png" width="300"></td>
      <td align="center"><img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/considerations_tutorials_files/tutorial1/coloc_dotplot_S8.png" width="300"></td>
    </tr>
    <tr>
      <td><b>S9</b></td>
      <td align="center"><img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/considerations_tutorials_files/tutorial1/plot_S9.png" width="300"></td>
      <td align="center"><img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/considerations_tutorials_files/tutorial1/hist_S9.png" width="300"></td>
      <td align="center"><img src="https://github.com/JEFworks/CRAWDAD/blob/main/docs/considerations_tutorials_files/tutorial1/coloc_dotplot_S9.png" width="300"></td>
    </tr>
  </tbody>
</table>


<p align="center"><b>Figure 2.</b> Simulation patterns showing identical spatial layout but different cell type proportions.</p>


## 5. How to interpret the results
<details> <summary><b>Explanation</b></summary>
Each dot-plot shows pairwise spatial associations between cell types.
Rows represent the reference cell type, and columns represent its neighbor.
Red indicates co-localization (positive z-scores), while blue indicates spatial depletion (negative z-scores).


When comparing the dot-plots across simulations, we see that the overall spatial pattern and total number of cells (N) remain constant.
Therefore, any differences in z-scores come purely from changes in cell-type proportions rather than real differences in spatial relationships.

In our simulations, the nature of the relationships (which cell types tend to co-localize or avoid each other) did not change, but the intensity of these signals did.
This happens because z-scores are influenced by the expected frequencies of cell-type interactions - when a cell type becomes more or less abundant, the denominator in the z-score calculation changes, inflating or deflating its apparent association strength.

To make fair comparisons in real datasets, one potential solution is to apply downsampling so that all samples have similar cell-type proportions before running CRAWDAD analysis.
This ensures that observed z-score differences reflect true biological variation, not differences in composition.
</details>
