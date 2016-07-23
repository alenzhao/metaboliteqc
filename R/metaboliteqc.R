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
