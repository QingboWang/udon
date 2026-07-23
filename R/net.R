# Internal environment: T matrices cached here after first load.
.udon_env <- new.env(parent = emptyenv())

# Loads T1, T2, genes1, genes2 from inst/extdata/ once per session.
.load_network <- function() {
  if (!is.null(.udon_env$T1)) return(invisible(NULL))

  message("Loading UDON network (runs once per session) ...")
  extdata <- system.file("extdata", package = "udon")

  genes1 <- readLines(file.path(extdata, "genes1.txt"))
  genes2 <- readLines(file.path(extdata, "genes2.txt"))

  n1 <- length(genes1)
  n2 <- length(genes2)

  T1_raw <- readBin(file.path(extdata, "T1.bin"),
                    what = "double", size = 4L, n = n1 * n1,
                    endian = "little")
  T2_raw <- readBin(file.path(extdata, "T2.bin"),
                    what = "double", size = 4L, n = n2 * n2,
                    endian = "little")

  .udon_env$T1     <- matrix(T1_raw, nrow = n1, ncol = n1, byrow = TRUE)
  .udon_env$T2     <- matrix(T2_raw, nrow = n2, ncol = n2, byrow = TRUE)
  .udon_env$genes1 <- genes1
  .udon_env$genes2 <- genes2

  message(sprintf("  Network ready: %d genes (causal), %d genes (dense).",
                  n1, n2))
  invisible(NULL)
}
