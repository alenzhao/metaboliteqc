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

#' Count number of groups with measurements below percentile per sample
#'
#' @param mat Matrix with samples as columns and metabolites as rows
#' @param by Factor of length \code{nrow(mat)} used to group rows
#' @param f Function to "condense" the measurements of a sample for a
#'     given group of rows to a single number
#'     (\code{\link[base]{mean}}, \code{\link[stats]{median}},
#'     \code{\link[base]{min}}, etc.)
#' @param percentile A number between 0 and 100.
#' @param \dots Passed on to \code{f}
#' @return A integer vector with \code{colnames(mat)} as names
#'     containing the number of categories, that is the levels of
#'     \code{by}, for which a sample has measurements (as "condensed"
#'     by \code{f}) below \code{percentile}.  The percentiles are
#'     computed separately for every metabolite over the "condensed"
#'     metabolite measurements of all samples.
#' @export
count_below_percentile <- function(mat, by, f, percentile, ...) {
    values <- columnwise(f, mat, by, ...)
    low <- t(apply(values, 1, function(x) {
        x <= quantile(x, percentile / 100, na.rm = TRUE)
    }))
    apply(low, 2, sum, na.rm = TRUE)
}

#' Show how many samples have low measurements in how many groups
#'
#' Plot the number of samples for which a given number of groups has
#' mean metabolite measurements below a given percentile.
#'
#' @param mat Matrix with samples as columns and metabolites as rows
#' @param by Factor of length \code{nrow(mat)} used to group rows
#' @param percentiles Integer vector with numbers between 0 and 100.
#' @param by_label Label used in plot annotations to refer to the
#'     groups defined by \code{by}
#' @param xlim Limits plotting range on x-axis
#' @param ylim Limits plotting range on y-axis
#' @details Algorithm: For every sample and every group of rows the
#'     mean metabolite measurement is computed.  This results in every
#'     sample having one value per group (the mean metabolite
#'     measurement).  We then count for every sample the number of
#'     groups in which the sample's mean metabolite measurement is
#'     less that or equal to the given percentile of the distribution
#'     of all mean metabolite measurements in that group.  Finally we
#'     plot the number of groups N against the number of samples with
#'     mean metabolite measurements less than or equal to the given
#'     percentile in N groups.
#' @return Returns \code{NULL} invisibly.
#' @export
plot_samples_with_low_metabolites <- function(mat, by, percentiles, by_label = NULL, xlim = NULL, ylim = NULL) {
    number_of_categories <- length(unique(by))
    xs <- lapply(percentiles, function(percentile) {
        x <- count_below_percentile(mat, by, mean, percentile, na.rm = TRUE)
        f <- factor(x, levels = seq(0, number_of_categories))
        tbl <- table(f)
        data.frame(categories = as.integer(names(tbl)), persons = as.integer(tbl))
    })
    if (is.null(by_label))
        by_label <- "categories"
    if (is.null(xlim))
        xlim <- c(0, number_of_categories)
    if (is.null(ylim))
        ylim <- c(0, do.call(max, lapply(xs, `[[`, "persons")))
    plot(1, type = "n", xlim = xlim, ylim = ylim,
        xlab = sprintf("Number of %s (N = %d) with low metabolites", by_label, number_of_categories),
        ylab = sprintf("Number of samples (N = %d)", ncol(mat)),
        main = sprintf("Counting samples by number of %s with low metabolites", by_label))
    pch <- rep_len(c(21, 24, 22, 25, 23), length(xs))
    transparent_gray <- rgb(0, 0, 0, .5)
    bg <- rep_len(c("white", transparent_gray), length(xs))
    for (i in seq_along(xs)) {
        x <- xs[[i]]
        lines(x$categories, x$persons, type = "b", pch = pch[i], bg = bg[i])
    }
    legend("topright", pch = pch[seq_along(xs)], pt.bg = bg[seq_along(xs)],
        legend = paste("mean <", percentiles, "% percentile"))
    invisible()
}

#' Replace missing values with column/row means
#'
#' @param mat Matrix
#' @return A matrix just like \code{mat} but with NAs in column/row k
#'     replaced by the mean over column/row k.
#' @export
replace_NAs_with_column_means <- function(mat) {
    apply(mat, 2, function(x) {
        missings <- is.na(x)
        if (any(missings))
            x[missings] <- mean(x, na.rm = TRUE)
        x
    })
}

#' @rdname replace_NAs_with_column_means
#' @export
replace_NAs_with_row_means <- function(mat) {
    t(replace_NAs_with_column_means(t(mat)))
}

#' Correct metabolite measurements for a factor
#'
#' Replace the measurements of every metabolite by the residuals from
#' regressing the metabolite on \code{x}, where \code{x} is converted
#' to a factor.
#'
#' @param mat Matrix with samples as columns and metabolites as rows
#' @param x Factor for which to correct the metabolite measurements in
#'     \code{mat}
#' @return A matrix of the same dimensions as \code{mat} where every
#'     metabolite's measurements were replaced by the residuals from
#'     regressing the metabolite on the \code{x} treated as a factor.
#'     If a metabolite consists only of NAs, the residuals for that
#'     metabolite will also be NAs.  The same happens when a
#'     metabolite has no valid variance (because there is only a
#'     single non-NA measurement) or when the metabolite has zero
#'     variance.
#' @export
correct_for_factor <- function(mat, x) {
    compute_residuals <- function(xs) {
        variance <- var(xs, na.rm = TRUE)
        if (all(is.na(xs)) || is.na(variance) || variance == 0)
            return(all_missing)
        residuals(lm(xs ~ factor(x), na.action = na.exclude))
    }
    all_missing <- rep_len(NA, ncol(mat))
    result <- t(apply(mat, 1, compute_residuals))
    colnames(result) <- colnames(mat)
    result
}
