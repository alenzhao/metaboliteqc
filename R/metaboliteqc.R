#' Copy missing values from one matrix to another
#'
#' @param from Matrix with missing values
#' @param to Matrix to which to copy missing values from \code{from}
#' @return A matrix like \code{to} but with missing values where
#'     \code{from} has missing values.
#' @export
copy_missing_values <- function(from, to) {
    if (!identical(dim(from), dim(to)))
        stop("Arguments FROM and TO must have the same dimensions.")
    if (!identical(dimnames(from), dimnames(to)))
        stop("Arguments FROM and TO must have the same dimnames.")
    missing_values <- is.na(from)
    to[missing_values] <- NA
    to
}

#' Count number of NAs
#'
#' @param x An object
#' @return Number of NAs in \code{x}.
#' @export
count_NAs <- function(x) {
    sum(is.na(x))
}
