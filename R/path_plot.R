#' Plot shortest paths from a source gene to core genes
#'
#' Finds and visualises all shortest paths (by hop count) from
#' \code{source_gene} to each gene in \code{core_genes} through the T1 causal
#' network. Paths sharing intermediate nodes are merged into a DAG view.
#' If no T1 path exists within \code{max_hops} and
#' \code{show_t2_fallback = TRUE}, the single shortest T2 dense-network path
#' is overlaid in dashed lines.
#'
#' Requires the packages \pkg{ggplot2} and \pkg{arrow}
#' (listed under \code{Suggests} in DESCRIPTION).
#'
#' @param source_gene     Character scalar. Starting gene symbol.
#' @param core_genes      Character vector. Target core gene symbols.
#' @param data_dir        Path to the directory containing
#'   \file{combined_edges_v5.parquet} (T1) and
#'   \file{step1_betamat_raw_v5.parquet} (T2).
#' @param max_hops        Maximum path length to consider (default \code{5}).
#' @param show_t2_fallback Logical. When \code{TRUE} (default), draw the
#'   shortest T2 path (dashed) for any core gene unreachable via T1 within
#'   \code{max_hops}.
#' @param output_file     PDF file path, or \code{NULL} to return the
#'   \code{ggplot} object without saving.
#'
#' @return A \code{ggplot} object (invisibly). When \code{output_file} is not
#'   \code{NULL} a PDF is written and a message is printed.
#'
#' @examples
#' \dontrun{
#' udon_path_plot(
#'   source_gene = "CELSR2",
#'   core_genes  = c("PCSK9", "LDLR", "ASGR1"),
#'   data_dir    = "/path/to/GxG_v5/data",
#'   output_file = "celsr2_paths.pdf"
#' )
#' }
#'
#' @export
udon_path_plot <- function(source_gene,
                           core_genes,
                           data_dir,
                           max_hops         = 5L,
                           show_t2_fallback = TRUE,
                           output_file      = NULL) {

  for (pkg in c("ggplot2", "arrow")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop(sprintf(
        "Package '%s' required. Install with: install.packages('%s')", pkg, pkg
      ))
  }

  data_dir   <- as.character(data_dir)
  max_hops   <- as.integer(max_hops)
  core_genes <- unique(core_genes)

  # ── Internal graph helpers (pure R, no igraph) ────────────────────────────

  # Build adjacency list: list named by "from", each element a char vector of "to"
  .adj <- function(df) split(df$to, df$from)

  # BFS distance from source; returns environment (hash) of distances
  .bfs_dist <- function(adj, src, max_h) {
    dist <- new.env(hash = TRUE, parent = emptyenv())
    assign(src, 0L, envir = dist)
    queue <- src
    while (length(queue)) {
      nd <- queue[1L]; queue <- queue[-1L]
      d  <- get(nd, envir = dist)
      if (d >= max_h) next
      for (nb in adj[[nd]])
        if (!exists(nb, envir = dist, inherits = FALSE)) {
          assign(nb, d + 1L, envir = dist)
          queue <- c(queue, nb)
        }
    }
    dist
  }

  # All simple paths of exactly min_d hops from src to tgt (stack-based DFS)
  .all_paths <- function(adj, src, tgt, min_d) {
    results <- list()
    stack   <- list(list(path = src, vis = src))
    while (length(stack)) {
      top   <- stack[[length(stack)]]; stack[[length(stack)]] <- NULL
      node  <- top$path[length(top$path)]
      depth <- length(top$path) - 1L
      if (node == tgt) {
        if (depth == min_d) results <- c(results, list(top$path))
        next
      }
      if (depth >= min_d) next
      for (nb in setdiff(adj[[node]], top$vis))
        stack <- c(stack, list(list(path = c(top$path, nb),
                                    vis  = c(top$vis,  nb))))
    }
    results
  }

  # Single shortest path via BFS (for T2 fallback); returns char vector or NULL
  .shortest_path <- function(adj, src, tgt, max_h) {
    prev  <- new.env(hash = TRUE, parent = emptyenv())
    dist  <- new.env(hash = TRUE, parent = emptyenv())
    assign(src, 0L, envir = dist)
    queue <- src; found <- FALSE
    while (length(queue) && !found) {
      nd <- queue[1L]; queue <- queue[-1L]
      d  <- get(nd, envir = dist)
      if (d >= max_h) next
      for (nb in adj[[nd]])
        if (!exists(nb, envir = dist, inherits = FALSE)) {
          assign(nb, d + 1L, envir = dist)
          assign(nb, nd,     envir = prev)
          queue <- c(queue, nb)
          if (nb == tgt) { found <- TRUE; break }
        }
    }
    if (!found) return(NULL)
    path <- tgt
    while (path[1L] != src) path <- c(get(path[1L], envir = prev), path)
    path
  }

  # ── Load and deduplicate T1 edges ─────────────────────────────────────────
  t1_raw <- arrow::read_parquet(file.path(data_dir, "combined_edges_v5.parquet"))
  t1_df  <- data.frame(
    from  = as.character(t1_raw$source),
    to    = as.character(t1_raw$target),
    beta  = as.numeric(t1_raw$beta),
    layer = as.character(t1_raw$layer),
    stringsAsFactors = FALSE
  )
  t1_df <- t1_df[t1_df$from != t1_df$to, ]
  t1_df <- t1_df[order(-abs(t1_df$beta)), ]
  t1_df <- t1_df[!duplicated(paste(t1_df$from, t1_df$to)), ]

  t1_beta_lup  <- setNames(t1_df$beta,  paste(t1_df$from, t1_df$to))
  t1_layer_lup <- setNames(t1_df$layer, paste(t1_df$from, t1_df$to))
  t1_nodes     <- unique(c(t1_df$from, t1_df$to))
  adj1         <- .adj(t1_df)

  if (!source_gene %in% t1_nodes)
    stop("'", source_gene, "' not found in T1 network.")

  # ── Find all shortest T1 paths per core gene ──────────────────────────────
  collected <- list()
  t2_needed <- character(0)

  t1_dist <- .bfs_dist(adj1, source_gene, max_hops)

  for (cg in core_genes) {
    if (cg == source_gene) {
      message("Skipping '", cg, "': same as source gene."); next
    }
    if (!exists(cg, envir = t1_dist, inherits = FALSE)) {
      msg <- if (!cg %in% t1_nodes)
        paste0("'", cg, "' not in T1 network")
      else
        paste0("No T1 path from '", source_gene, "' to '", cg,
               "' within ", max_hops, " hops")
      message(msg, if (show_t2_fallback) "; will try T2 fallback." else ".")
      if (show_t2_fallback) t2_needed <- c(t2_needed, cg)
      next
    }
    min_d <- get(cg, envir = t1_dist)
    paths <- .all_paths(adj1, source_gene, cg, min_d)

    if (length(paths) > 20L)
      message(length(paths), " shortest paths of length ", min_d,
              " to '", cg, "' — all will be drawn.")

    for (p in paths) {
      for (i in seq_len(length(p) - 1L)) {
        key <- paste(p[i], p[i + 1L])
        collected[[length(collected) + 1L]] <- data.frame(
          from   = p[i],
          to     = p[i + 1L],
          beta   = t1_beta_lup[[key]],
          layer  = t1_layer_lup[[key]],
          via_t2 = FALSE,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  # ── T2 fallback paths ──────────────────────────────────────────────────────
  if (show_t2_fallback && length(t2_needed) > 0L) {
    t2_raw <- arrow::read_parquet(
      file.path(data_dir, "step1_betamat_raw_v5.parquet")
    )
    t2_df <- data.frame(
      from  = as.character(t2_raw$source),
      to    = as.character(t2_raw$target),
      beta  = as.numeric(t2_raw$beta_raw),
      layer = "dense",
      stringsAsFactors = FALSE
    )
    t2_df  <- t2_df[t2_df$from != t2_df$to, ]
    t2_df  <- t2_df[order(-abs(t2_df$beta)), ]
    t2_df  <- t2_df[!duplicated(paste(t2_df$from, t2_df$to)), ]
    t2_beta_lup <- setNames(t2_df$beta, paste(t2_df$from, t2_df$to))
    adj2        <- .adj(t2_df)

    for (cg in t2_needed) {
      path <- .shortest_path(adj2, source_gene, cg, max_hops)
      if (is.null(path)) {
        message("No T2 path to '", cg, "' within ", max_hops, " hops; skipping.")
        next
      }
      for (i in seq_len(length(path) - 1L)) {
        key <- paste(path[i], path[i + 1L])
        collected[[length(collected) + 1L]] <- data.frame(
          from   = path[i],
          to     = path[i + 1L],
          beta   = t2_beta_lup[[key]],
          layer  = "dense",
          via_t2 = TRUE,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(collected) == 0L)
    stop("No paths found. Increase max_hops or check gene name spelling.")

  edge_df <- unique(do.call(rbind, collected))

  # ── Node column positions (BFS distance from source over path subgraph) ───
  all_nodes <- unique(c(edge_df$from, edge_df$to))
  col_of    <- setNames(rep(NA_integer_, length(all_nodes)), all_nodes)
  col_of[source_gene] <- 0L
  queue <- source_gene
  while (length(queue)) {
    nd   <- queue[1L]; queue <- queue[-1L]
    nbrs <- edge_df$to[edge_df$from == nd]
    for (nb in nbrs) if (is.na(col_of[nb])) {
      col_of[nb] <- col_of[nd] + 1L; queue <- c(queue, nb)
    }
  }

  # ── Node table with layout positions ──────────────────────────────────────
  node_df <- data.frame(
    gene      = all_nodes,
    col       = col_of[all_nodes],
    is_source = all_nodes == source_gene,
    is_core   = all_nodes %in% core_genes,
    stringsAsFactors = FALSE
  )
  node_df <- node_df[order(node_df$col, !node_df$is_core, node_df$gene), ]

  node_df$row_in_col <- unlist(lapply(
    split(seq_len(nrow(node_df)), node_df$col), seq_along
  ))

  COL_SP <- 3.2; ROW_SP <- 1.2
  node_df$x <- node_df$col * COL_SP
  col_n      <- table(node_df$col)
  node_df$y  <- vapply(seq_len(nrow(node_df)), function(i) {
    n <- as.integer(col_n[[as.character(node_df$col[i])]])
    (node_df$row_in_col[i] - (n + 1) / 2) * ROW_SP
  }, numeric(1))

  gx <- setNames(node_df$x, node_df$gene)
  gy <- setNames(node_df$y, node_df$gene)

  edge_df$x0 <- gx[edge_df$from]; edge_df$y0 <- gy[edge_df$from]
  edge_df$x1 <- gx[edge_df$to];   edge_df$y1 <- gy[edge_df$to]

  # ── Shorten edge endpoints to node-box boundary ───────────────────────────
  BOX_HW <- 0.58; BOX_HH <- 0.21

  shorten <- function(x0, y0, x1, y1) {
    dx <- x1 - x0; dy <- y1 - y0
    d  <- sqrt(dx^2 + dy^2)
    if (d < 1e-9) return(c(x0, y0, x1, y1))
    ux <- dx / d; uy <- dy / d
    off <- min(if (abs(ux) > 1e-9) BOX_HW / abs(ux) else Inf,
               if (abs(uy) > 1e-9) BOX_HH / abs(uy) else Inf)
    c(x0 + ux * off, y0 + uy * off, x1 - ux * off, y1 - uy * off)
  }
  seg <- t(mapply(shorten, edge_df$x0, edge_df$y0, edge_df$x1, edge_df$y1))
  edge_df$xs0 <- seg[, 1L]; edge_df$ys0 <- seg[, 2L]
  edge_df$xs1 <- seg[, 3L]; edge_df$ys1 <- seg[, 4L]
  edge_df$xm  <- (edge_df$xs0 + edge_df$xs1) / 2
  edge_df$ym  <- (edge_df$ys0 + edge_df$ys1) / 2
  edge_df$blabel <- sprintf("%+.3f", edge_df$beta)

  # ── Edge aesthetics ────────────────────────────────────────────────────────
  edge_df$ecol <- ifelse(!is.na(edge_df$beta) & edge_df$beta >= 0,
                          "#3B6BA5", "#C84B31")
  edge_df$elty <- ifelse(edge_df$via_t2, "dashed", "solid")
  babs <- abs(edge_df$beta); bmax <- max(babs, na.rm = TRUE)
  edge_df$elwd <- if (bmax > 0) 0.55 + babs / bmax * 1.65 else rep(1.0, nrow(edge_df))

  pos_e <- edge_df[!is.na(edge_df$beta) & edge_df$beta >= 0, ]
  neg_e <- edge_df[!is.na(edge_df$beta) & edge_df$beta <  0, ]

  TBAR <- 0.10
  if (nrow(neg_e) > 0L) {
    dx  <- neg_e$xs1 - neg_e$xs0; dy <- neg_e$ys1 - neg_e$ys0
    d   <- sqrt(dx^2 + dy^2); d[d < 1e-9] <- 1
    ux  <- dx / d; uy <- dy / d
    neg_e$tb_x0 <- neg_e$xs1 - uy * TBAR
    neg_e$tb_y0 <- neg_e$ys1 + ux * TBAR
    neg_e$tb_x1 <- neg_e$xs1 + uy * TBAR
    neg_e$tb_y1 <- neg_e$ys1 - ux * TBAR
  }

  # ── Node aesthetics ────────────────────────────────────────────────────────
  t2_int <- setdiff(
    unique(c(edge_df$from[edge_df$via_t2], edge_df$to[edge_df$via_t2])),
    c(source_gene, core_genes)
  )
  node_df$fill <- ifelse(node_df$is_source, "#E87722",
                   ifelse(node_df$is_core,            "#2B4590",
                   ifelse(node_df$gene %in% t2_int,   "#6FAF89",
                                                       "#5B9BD5")))

  # ── Plot ───────────────────────────────────────────────────────────────────
  p <- ggplot2::ggplot() +

    ggplot2::geom_segment(
      data = pos_e,
      ggplot2::aes(x = xs0, y = ys0, xend = xs1, yend = ys1,
                   colour = ecol, linewidth = elwd, linetype = elty),
      arrow   = ggplot2::arrow(length = ggplot2::unit(0.18, "cm"), type = "closed"),
      lineend = "round"
    ) +

    {if (nrow(neg_e) > 0L)
      ggplot2::geom_segment(
        data = neg_e,
        ggplot2::aes(x = xs0, y = ys0, xend = xs1, yend = ys1,
                     colour = ecol, linewidth = elwd, linetype = elty),
        lineend = "round"
      )
    } +

    {if (nrow(neg_e) > 0L)
      ggplot2::geom_segment(
        data = neg_e,
        ggplot2::aes(x = tb_x0, y = tb_y0, xend = tb_x1, yend = tb_y1,
                     colour = ecol, linewidth = elwd),
        lineend = "round"
      )
    } +

    ggplot2::geom_text(
      data     = edge_df,
      ggplot2::aes(x = xm, y = ym, label = blabel),
      size     = 2.4,
      colour   = "gray25",
      nudge_y  = 0.13,
      fontface = "italic"
    ) +

    ggplot2::geom_label(
      data          = node_df,
      ggplot2::aes(x = x, y = y, label = gene, fill = fill),
      colour        = "white",
      fontface      = "bold",
      size          = 3.2,
      label.r       = ggplot2::unit(0.28, "lines"),
      label.padding = ggplot2::unit(0.42, "lines"),
      label.size    = 0.3
    ) +

    ggplot2::scale_colour_identity() +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_linewidth_identity() +
    ggplot2::scale_linetype_identity() +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_void() +

    ggplot2::labs(
      title = sprintf("%s  →  {%s}",
                      source_gene, paste(core_genes, collapse = ", ")),
      subtitle = sprintf(
        "All shortest paths — T1 causal network (max %d hops)", max_hops
      ),
      caption = paste0(
        "Blue → : positive  |  Orange ⊣ : inhibition  |  Thicker = larger |β|",
        if (any(edge_df$via_t2, na.rm = TRUE))
          "  |  Dashed / green nodes: dense network (T2 fallback)" else ""
      )
    ) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face   = "bold", size = 13,
                                            hjust  = 0.5,
                                            margin = ggplot2::margin(b = 4)),
      plot.subtitle = ggplot2::element_text(size   = 9, hjust = 0.5,
                                            colour = "gray45",
                                            margin = ggplot2::margin(b = 14)),
      plot.caption  = ggplot2::element_text(size   = 7.5, colour = "gray50",
                                            hjust  = 0),
      plot.margin   = ggplot2::margin(18, 40, 14, 40)
    )

  # ── Save ───────────────────────────────────────────────────────────────────
  if (!is.null(output_file)) {
    n_cols <- max(node_df$col, na.rm = TRUE) + 1L
    width  <- max(7.0, n_cols * 3.4)
    height <- max(3.5, max(as.integer(col_n)) * ROW_SP * 1.6 + 2.5)
    out_dir <- dirname(output_file)
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    ggplot2::ggsave(output_file, p, width = width, height = height,
                    device = "pdf")
    message("Saved: ", output_file)
  }

  invisible(p)
}
