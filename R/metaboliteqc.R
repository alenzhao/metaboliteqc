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
#' @param main Title
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
plot_samples_with_low_metabolites <- function(mat, by, percentiles, by_label = NULL, xlim = NULL, ylim = NULL, main = NULL) {
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
    if (is.null(main))
        main <- sprintf("Counting samples by number of %s with low metabolites", by_label)
    plot(1, type = "n", xlim = xlim, ylim = ylim,
        xlab = sprintf("Number of %s (N = %d) with low metabolites", by_label, number_of_categories),
        ylab = sprintf("Number of samples (N = %d)", ncol(mat)),
        main = main)
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

#' Select rows/columns of a matrix using a predicate
#'
#' @param predicate A function taking a row/column as first argument
#'     and returning \code{TRUE} or \code{FALSE}
#' @param mat Matrix
#' @param margin Whether to select rows or columns (1 = rows, 2 =
#'     columns).
#' @param \dots Passed on to \code{predicate}
#' @return If \code{margin} is 1, a matrix with the same number of
#'     columns as \code{mat} and including all rows of \code{mat} for
#'     which \code{predicate} returned \code{TRUE}.  If \code{margin}
#'     is 2, a matrix with the same number of rows as \code{mat} and
#'     including all columns of \code{mat} for which \code{predicate}
#'     returned \code{TRUE}.
#' @export
select_if <- function(predicate, mat, margin, ...) {
    selected <- apply(mat, margin, predicate, ...)
    if (margin == 1)
        mat[selected, , drop = FALSE]
    else
        mat[, selected, drop = FALSE]
}

#' Draw outliers for boxplot
#'
#' @param outliers y-coordinates of outliers
#' @param at Length-1 x-coordinate of outliers
#' @return Returns \code{NULL} invisibly.
draw_outliers <- function(outliers, at) {
    if (length(outliers) == 0)
        return()
    points(rep_len(at, length(outliers)), outliers, pch = ".")
    invisible()
}

#' Draw whiskers for boxplot
#'
#' @param stats Numeric vector of length 5 containing the extreme of
#'     the lower whisker, the lower `hinge', the media, the upper
#'     `hinge' and the extreme of the upper whisker (see also
#'     \code{\link[grDevices]{boxplot.stats}}).
#' @param width Width of horizontal endings of whiskers
#' @param at x-coordinate of whiskers
#' @return Returns \code{NULL} invisibly.
draw_whiskers <- function(stats, width, at) {
    segments(x0 = at, y0 = stats[1], y1 = stats[2])
    segments(x0 = at, y0 = stats[4], y1 = stats[5])
    segments(x0 = at - width / 2, x1 = at + width / 2, y0 = stats[1])
    segments(x0 = at - width / 2, x1 = at + width / 2, y0 = stats[5])
    invisible()
}

#' Draw box for boxplot
#'
#' @param stats Numeric vector of length 5 containing the extreme of
#'     the lower whisker, the lower `hinge', the media, the upper
#'     `hinge' and the extreme of the upper whisker (see also
#'     \code{\link[grDevices]{boxplot.stats}}).
#' @param width Width of the box
#' @param at x-coordinate of the (center of the) box
#' @return Returns \code{NULL} invisibly.
draw_box <- function(stats, width, at) {
    rect(xleft = at - width / 2,
        xright = at + width / 2,
        ybottom = stats[2],
        ytop = stats[4],
        border = "black")
    segments(x0 = at - width / 2,
        x1 = at + width / 2,
        y0 = stats[3],
        y1 = stats[3])
    invisible()
}

#' Draw boxplot
#'
#' @param stats Output of \code{\link[grDevices]{boxplot.stats}}
#' @param at x-coordinate of the boxplot
#' @return Returns \code{NULL} invisibly.
draw_boxplot <- function(stats, at) {
    draw_box(stats$stats, .6, at)
    draw_whiskers(stats$stats, .4, at)
    draw_outliers(stats$out, at)
    invisible()
}

#' Draw boxplots for subsets of columns
#'
#' @param mat Matrix
#' @param by Factor of length \code{ncol(mat)} used to group columns
#'     into subsets
#' @param ylab Label for y-axis
#' @param ylim Range of y-axis
#' @return Returns \code{NULL} invisibly.
#' @export
byboxplot <- function(mat, by, ylab, ylim = NULL) {
    xs <- lapply(split.data.frame(t(mat), by), as.vector)
    stats <- lapply(xs, boxplot.stats)
    if (is.null(ylim))
        ylim <- range(xs, na.rm = TRUE)
    xlim <- c(1, length(xs))
    par(mar = c(0, 3, 0, 3))
    plot(1, xlim = xlim, ylim = ylim, type = "n", axes = FALSE, ann = FALSE)
    for (i in seq_along(stats)) {
        draw_boxplot(stats[[i]], at = i)
    }
    mtext(ylab, side = 2, line = 1)
    axis(4)
    box()
    invisible()
}

#' Draw boxplots for rectangular subsets of matrix
#'
#' @param mat Matrix
#' @param by.x Factor to group columns of \code{mat} into subsets
#' @param by.y Factor to group rows of \code{mat} into subsets
#' @param main Title
#' @param ylim y-axis range to use for all subplots
#' @details The \code{by.x} argument has suffix "x" because the
#'     subsets of columns that \code{by.x} creates will be drawn along
#'     the x-axis.  Similarly, the \code{by.y} argument has suffix "y"
#'     because the subsets of rows that \code{by.y} creates will be
#'     drawn along the y-axis.
#' @return Returns \code{NULL} invisibly.
#' @export
bybyboxplot <- function(mat, by.x, by.y, main = NULL, ylim = NULL) {
    number_of_plots <- length(unique(by.y))
    ds <- split.data.frame(mat, by.y)
    if (is.null(main))
        main <- "Boxplots"
    par(mfrow = c(number_of_plots, 1),
        oma = c(3, 0, 6, 0))
    for (i in seq_along(ds)) {
        byboxplot(ds[[i]], by.x, names(ds)[i], ylim = ylim)
        if (i == 1L)
            axis(3, at = seq_along(levels(by.x)), labels = levels(by.x), cex.axis = 1)
    }
    axis(1, at = seq_along(levels(by.x)), labels = levels(by.x), cex.axis = 1)
    title(main, line = 4, outer = TRUE, cex.main = 1.5)
    invisible()
}

#' Find indices of cells that satisfy a predicate
#'
#' @param mat Matrix
#' @param margin Whether to apply \code{predicate} to rows or to
#'     columns (1 = rows, 2 = columns)
#' @param predicate A function that takes a vector (a row or a column
#'     of \code{mat}) and returns a logical vector of the same length
#' @param \dots Passed on to \code{predicate}
#' @details Note that \code{predicate} can use information on all
#'     elements in a row (column) to decide whether a particular
#'     element of the row (column) satisfies the predicate or not.
#' @return A two-column character matrix that can be used as index to
#'     \code{mat} to extract elements that satisfy \code{predicate}.
#' @export
find_cell_indices <- function(mat, margin, predicate, ...) {
    x <- apply(mat, margin, predicate, ...)
    if (margin == 1L)
        x <- t(x)
    x[is.na(x)] <- FALSE
    cbind(rownames(x)[row(x)[x]],
        colnames(x)[col(x)[x]])
}

#' Save dataset as tab-delimited plain text file
#'
#' @param x Dataset
#' @param file Filename
#' @param \dots Passed on to \code{\link[utils]{write.table}}
#' @details This function is a wrapper around
#'     \code{\link[utils]{write.table}}.  Refer to the latter's
#'     documentation for more details.
#' @return Returns \code{NULL} invisibly.
#' @export
write.delim <- function(x, file, ...) {
    write.table(x, file, quote = FALSE, sep = "\t", ...)
    invisible()
}

#' Mark outliers
#'
#' @param x A vector of numbers
#' @return A logical vector of the same length as \code{x} where
#'     elements are \code{TRUE} if the corresponding element of
#'     \code{x} is an outlier among the the values in \code{x}, and
#'     otherwise \code{FALSE}.
#' @export
mark_outliers <- function(x) {
    mean <- mean(x, na.rm = TRUE)
    sd <- sd(x, na.rm = TRUE)
    abs(x - mean) > 4 * sd
}

jpeg_width <- 800
jpeg_height <- 1100
jpeg_quality <- 100

plot_NAs_by_run_day <- function(main, filename, mat, by.x) {
    jpeg(filename, width = jpeg_width, height = 300, quality = jpeg_quality)
    x <- apply(mat, 2, count_NAs)
    bybyboxplot(t(x), main = main, by.x = by.x, by.y = "Number of NAs")
    dev.off()
}

plot_NAs_by_run_day_and_pathway <- function(main, filename, mat, run_days, pathways, ylim) {
    jpeg(filename, width = jpeg_width, height = jpeg_height, quality = jpeg_quality)
    x <- t(columnwise(count_NAs, t(mat), run_days))
    bybyboxplot(x, main = main,
        by.x = ordered(as.integer(colnames(x))),
        by.y = pathways, ylim = ylim)
    dev.off()
}

plot_measurements_by_run_day <- function(main, filename, mat, run_days, ylab) {
    jpeg(filename, width = jpeg_width, height = 300, quality = jpeg_quality)
    bybyboxplot(mat, by.x = run_days, by.y = ylab, main = main)
    dev.off()
}

plot_measurements_by_run_day_and_pathway <- function(main, filename, mat, run_days, pathways) {
    jpeg(filename, width = jpeg_width, height = jpeg_height, quality = jpeg_quality)
    bybyboxplot(mat, ylim = c(-3, 3), by.x = run_days, by.y = pathways, main = main)
    dev.off()
}

plot_NAs_per_sample <- function(main, filename, mat) {
    jpeg(filename, height = 350)
    hist_NAs(mat, 2, main = main,
        xlab = sprintf("Number of missing metabolites (N = %d)", nrow(mat)),
        ylab = sprintf("Number of samples (N = %d)", ncol(mat)))
    dev.off()
}

is_in_lower_tail <- function(x, percentage) {
    x < quantile(x, percentage, na.rm = TRUE)
}

find_percent_samples_with_low_metabolites <- function(percentage, mat) {
    is_low_in_metabolite <- t(apply(mat, 1, is_in_lower_tail, percentage))
    low_metabolites <- colSums(is_low_in_metabolite, na.rm = TRUE)
    non_NA_metabolites <- colSums(!is.na(mat))
    percent_low_metabolites <- 100 * low_metabolites / non_NA_metabolites
    # Percentage of samples with at least 0, 1, ... percent low
    # metabolites.
    vapply(0:100, function(percent) {
        sum(percent_low_metabolites >= percent, na.rm = TRUE)
    }, integer(1)) / ncol(mat) * 100
}

plot_percent_samples_with_low_metabolites <- function(main, filename, mat, percentages, xlim = NULL, ylim = NULL) {
    number_of_samples <- ncol(mat)
    xs <- lapply(percentages, find_percent_samples_with_low_metabolites, mat)
    if (is.null(xlim))
        xlim <- c(0, 100)
    if (is.null(ylim))
        ylim <- c(0, 100)
    cols <- rep_len(c("magenta", "green", "blue", "orange"), length(xs))
    jpeg(filename, width = jpeg_width, height = 350)
    par(mar = c(5, 4, 4, 4) + .5)
    plot(1, type = "n", axes = FALSE, ann = FALSE, xaxs = "i", xlim = xlim, ylim = ylim)
    abline(h = 0, lty = "dotted")
    Map(function(x, col) lines(seq(0, 100), x, col = col), xs, cols)
    title(main)
    mtext("minimum percentage of low metabolites", side = 1, line = 3)
    mtext("number of samples", side = 4, line = 3)
    mtext("percentage of samples", side = 2, line = 3)
    axis(1)
    axis(2)
    at <- pretty(ylim / 100 * number_of_samples)
    axis(4, at = at / number_of_samples * 100, at)
    legend("topright", fill = cols, bty = "n", inset = .01,
        title = "where low means less than the",
        legend = sprintf("%d %% percentile", 100 * percentages))
    box()
    dev.off()
}

plot_samples_with_low_pathways <- function(main, filename, mat, pathways) {
    jpeg(filename, height = 350)
    plot_samples_with_low_metabolites(mat, by = pathways,
        percentiles = c(1, 3, 5, 10), by_label = "pathways",
        xlim = c(3, 9), ylim = c(0, 100), main = main)
    dev.off()
}
