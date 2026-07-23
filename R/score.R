#' Score genes by network propagation to a candidate set
#'
#' Given a set of genes associated with a trait, scores all proteins in the
#' network by how strongly their activity propagates to those candidates.
#' Ranks by absolute score (rank 1 = strongest influence).
#'
#' @param candidates Character vector of gene symbols associated with the trait
#'   (e.g. GWAS hits, eQTL targets). Genes absent from the network are ignored
#'   with a warning.
#' @param lof_betas Named numeric vector of LoF effect sizes, keyed by gene
#'   symbol (e.g. `c(PCSK9 = -0.985, LDLR = 0.284)`). The weight applied to
#'   candidate `c` is `-lof_betas[c]`, matching the UDON formula
#'   `w(c) = -beta_LoF[c]`. Candidates missing from this vector fall back to
#'   unit weight (1.0). If `NULL` (default), all candidates receive unit weight
#'   (UDON_approx mode).
#'
#' @return A data.frame with one row per network gene, sorted by rank:
#'   \describe{
#'     \item{gene}{Gene symbol.}
#'     \item{score}{Signed UDON score. Positive = activity correlates with
#'       candidate set in the network's propagation direction.}
#'     \item{abs_score}{Absolute score used for ranking.}
#'     \item{rank}{Rank by abs_score descending; rank 1 = strongest influence.}
#'     \item{is_candidate}{Logical; TRUE if the gene is in \code{candidates}.}
#'   }
#'
#' @examples
#' # UDON_approx: unit weights
#' res <- score(c("PCSK9", "LDLR", "ANGPTL3"))
#' head(res, 10)
#'
#' # UDON: provide LoF betas where available; missing ones fall back to 1.0
#' res <- score(
#'   c("PCSK9", "LDLR", "ANGPTL3"),
#'   lof_betas = c(PCSK9 = -0.985, LDLR = 0.284, ANGPTL3 = -0.314)
#' )
#' head(res, 10)
#'
#' @export
score <- function(candidates, lof_betas = NULL) {
  stopifnot(is.character(candidates), length(candidates) > 0)
  if (!is.null(lof_betas)) {
    stopifnot(is.numeric(lof_betas), !is.null(names(lof_betas)))
  }

  .load_network()

  T1     <- .udon_env$T1
  T2     <- .udon_env$T2
  genes1 <- .udon_env$genes1
  genes2 <- .udon_env$genes2

  cands <- unique(candidates)

  # Warn about candidates absent from both networks
  missing <- setdiff(cands, union(genes1, genes2))
  if (length(missing) > 0) {
    warning(sprintf(
      "%d candidate(s) not found in either network and will not contribute ",
      "to scores: %s", length(missing), paste(sort(missing), collapse = ", ")
    ), call. = FALSE)
  }

  # Weight function: -lof_betas[g] if provided, else 1.0
  weight <- function(g) {
    if (!is.null(lof_betas) && g %in% names(lof_betas)) -lof_betas[[g]]
    else 1.0
  }

  # Weight vectors aligned to each network's gene order
  w1 <- ifelse(genes1 %in% cands, vapply(genes1, weight, numeric(1)), 0.0)
  w2 <- ifelse(genes2 %in% cands, vapply(genes2, weight, numeric(1)), 0.0)

  # Score = T %*% w (matrix-vector multiply)
  s1 <- as.numeric(T1 %*% w1)
  s2 <- as.numeric(T2 %*% w2)

  # Merge: sum contributions from both networks
  df1 <- data.frame(gene = genes1, s1 = s1, stringsAsFactors = FALSE)
  df2 <- data.frame(gene = genes2, s2 = s2, stringsAsFactors = FALSE)
  df  <- merge(df1, df2, by = "gene", all = TRUE)
  df[is.na(df$s1), "s1"] <- 0.0
  df[is.na(df$s2), "s2"] <- 0.0
  df$score     <- df$s1 + df$s2
  df$abs_score <- abs(df$score)
  df$rank      <- rank(-df$abs_score, ties.method = "min")
  df$is_candidate <- df$gene %in% cands
  df <- df[order(df$rank), c("gene","score","abs_score","rank","is_candidate")]
  rownames(df) <- NULL
  df
}

#' List all gene symbols in the network
#'
#' @return Character vector of gene symbols (union of both T matrices).
#' @export
network_genes <- function() {
  .load_network()
  sort(union(.udon_env$genes1, .udon_env$genes2))
}
