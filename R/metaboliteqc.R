#' Copy NAs from one matrix to another
#'
#' @param from Matrix with NAs
#' @param to Matrix to which to copy NAs from \code{from}
#' @return A matrix like \code{to} but with NAs where \code{from} has
#'     NAs.
#' @export
copy_NAs <- function(from, to) {
    if (!identical(dim(from), dim(to)))
        stop("Arguments FROM and TO must have the same dimensions.")
    if (!identical(dimnames(from), dimnames(to)))
        stop("Arguments FROM and TO must have the same dimnames.")
    missing_values <- is.na(from)
    to[missing_values] <- NA
    to
}

#' Find number of missing values and missingness rate
#'
#' @param mat Matrix
#' @param by An integer.  Dimension for which to compute missingness
#'     statistics (1 = rows, 2 = columns).  Same as \code{MARGIN}
#'     argument to \code{\link[base]{apply}}.
#' @return A data frame with columns "id", "missings", and
#'     "missingness_rate" containing row (column) names, number of
#'     missings per row (column), and missingness rate per row
#'     (column) if \code{by} is 1 (2).
#' @export
summarize_NAs <- function(mat, by) {
    id <- dimnames(mat)[[by]]
    if (is.null(id))
        id <- as.character(seq_len(dim(mat)[by]))
    missings <- apply(mat, by, count_NAs)
    missingness_rate <- missings / dim(mat)[by]
    data.frame(id, missings, missingness_rate,
        stringsAsFactors = FALSE, row.names = NULL)
}

#' Count number of NAs
#'
#' @param x An object
#' @return Number of NAs in \code{x}.
#' @export
count_NAs <- function(x) {
    sum(is.na(x))
}

#' Histogram of number of NAs
#'
#' @param mat Matrix
#' @param by An integer.  Whether to count NAs by rows or by columns
#'     (1 = rows, 2 = columns).
#' @param \dots Arguments passed on to \code{\link[graphics]{hist}}
#' @param main Title for plot
#' @param xlab Label for x-axis
#' @return Returns \code{NULL} invisibly.
#' @export
hist_NAs <- function(mat, by, ..., main = NULL, xlab = NULL) {
    if (is.null(main))
        main <- "Histogram of number of NAs"
    if (is.null(xlab))
        xlab <- "Number of NAs"
    d <- summarize_NAs(mat, by)
    hist(d$missings, ..., main = main, xlab = xlab)
    rug(d$missings, side = 3, ticksize = -.02)
    box()
    invisible()
}

#' Apply a function columnwise to groups of rows
#'
#' @param f Function
#' @param mat Matrix
#' @param by Factor of length \code{nrow(mat)} used to group rows
#' @param \dots Arguments passed on to \code{f}
#' @return A matrix with \code{ncol(mat)} columns and
#'     \code{length(unique(by)) * n} rows, where \code{n} is the
#'     length of the result returned by \code{f}.
#' @export
columnwise <- function(f, mat, by, ...) {
    g <- function(x) {
        do.call(f, c(list(x), list(...)))
    }
    matrices <- split.data.frame(mat, by)
    do.call(rbind, lapply(matrices, apply, 2, g))
}
