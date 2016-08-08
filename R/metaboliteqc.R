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
plot_samples_with_low_median_metabolites <- function(mat, by, percentiles, by_label = NULL, xlim = NULL, ylim = NULL, main = NULL) {
    number_of_categories <- length(unique(by))
    xs <- lapply(percentiles, function(percentile) {
        x <- count_below_percentile(mat, by, median, percentile, na.rm = TRUE)
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
        title = "where low means", title.adj = 0,
        legend = paste("mean <", percentiles, "% percentile"), bty = "n")
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
#' @param width Width of box
#' @return Returns \code{NULL} invisibly.
draw_boxplot <- function(stats, at, width) {
    draw_box(stats$stats, width, at)
    draw_whiskers(stats$stats, width * 2/3, at)
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
    xlim <- c(1, length(xs)) + c(-.5, .5)
    par(mar = c(0, 3, 0, 3))
    plot(1, xlim = xlim, ylim = ylim, type = "n", axes = FALSE, ann = FALSE)
    for (i in seq_along(stats)) {
        draw_boxplot(stats[[i]], at = i, width = .6)
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
#' @param xlab Label for x-axis
#' @details The \code{by.x} argument has suffix "x" because the
#'     subsets of columns that \code{by.x} creates will be drawn along
#'     the x-axis.  Similarly, the \code{by.y} argument has suffix "y"
#'     because the subsets of rows that \code{by.y} creates will be
#'     drawn along the y-axis.
#' @return Returns \code{NULL} invisibly.
#' @export
bybyboxplot <- function(mat, by.x, by.y, main = NULL, ylim = NULL, xlab = NULL) {
    number_of_plots <- length(unique(by.y))
    ds <- split.data.frame(mat, by.y)
    if (is.null(main))
        main <- "Boxplots"
    par(mfrow = c(number_of_plots, 1), oma = c(4, 0, 4, 0), cex = 1)
    for (i in seq_along(ds)) {
        byboxplot(ds[[i]], by.x, names(ds)[i], ylim = ylim)
    }
    axis(1, at = seq_along(levels(by.x)), labels = levels(by.x), cex.axis = 1)
    title(main, outer = TRUE)
    if (!is.null(xlab))
        title(xlab = xlab, line = 3, outer = TRUE)
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

#' @export
plot_NAs_by_run_day <- function(main, filename, mat, by.x, width) {
    jpeg(filename, width = width, height = 300, quality = jpeg_quality)
    x <- apply(mat, 2, count_NAs)
    bybyboxplot(t(x), main = main, by.x = by.x, by.y = "Number of NAs", xlab = "Run day")
    dev.off()
}

#' @export
plot_NAs_by_run_day_and_pathway <- function(main, filename, mat, run_days, pathways, ylim, width) {
    jpeg(filename, width = width, height = jpeg_height, quality = jpeg_quality)
    x <- t(columnwise(count_NAs, t(mat), run_days))
    bybyboxplot(x, main = main,
        by.x = ordered(as.integer(colnames(x))),
        by.y = pathways, ylim = ylim, xlab = "Run day")
    dev.off()
}

#' @export
plot_measurements_by_run_day <- function(main, filename, mat, run_days, ylab, width) {
    jpeg(filename, width = width, height = 300, quality = jpeg_quality)
    bybyboxplot(mat, by.x = run_days, by.y = ylab, main = main, xlab = "Run day")
    dev.off()
}

#' @export
plot_measurements_by_run_day_and_pathway <- function(main, filename, mat, run_days, pathways, width, ylim = NULL) {
    jpeg(filename, width = width, height = jpeg_height, quality = jpeg_quality)
    bybyboxplot(mat, ylim = ylim, by.x = run_days, by.y = pathways, main = main, xlab = "Run day")
    dev.off()
}

#' @export
plot_NAs_per_sample <- function(main, filename, mat, width) {
    jpeg(filename, width = width, height = 350)
    hist_NAs(mat, 2, main = main,
        xlab = sprintf("Number of missing metabolites (N = %d)", nrow(mat)),
        ylab = sprintf("Number of samples (N = %d)", ncol(mat)))
    dev.off()
}

#' @export
is_in_lower_tail <- function(x, percentage) {
    x < quantile(x, percentage, na.rm = TRUE)
}

#' @export
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

#' @export
plot_percent_samples_with_low_metabolites <- function(main, filename, mat, xlim = NULL, ylim = NULL, width) {
    number_of_samples <- ncol(mat)
    percentages <- c(.1, .05, .03, .01)
    xs <- lapply(percentages, find_percent_samples_with_low_metabolites, mat)
    if (is.null(xlim))
        xlim <- c(0, 100)
    if (is.null(ylim))
        ylim <- c(0, 100)
    cols <- rep_len(c("magenta", "green", "blue", "orange"), length(xs))
    jpeg(filename, width = width, height = 350)
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
        title = "where low means", title.adj = 0,
        legend = sprintf("< %d %% percentile", 100 * percentages))
    box()
    dev.off()
}

#' @export
plot_samples_with_low_pathways <- function(main, filename, mat, pathways, width, xlim = NULL, ylim = NULL) {
    if (is.null(xlim))
        xlim <- c(0, length(unique(pathways)))
    if (is.null(ylim))
        ylim <- c(0, ncol(mat))
    jpeg(filename, width = width, height = 350)
    plot_samples_with_low_median_metabolites(mat, by = pathways,
        percentiles = c(10, 5, 3, 1), by_label = "pathways",
        xlim = xlim, ylim = ylim, main = main)
    dev.off()
}

#' @export
save_NAs_per_sample <- function(filename, mat) {
    d <- summarize_NAs(mat, 2)
    names(d)[1] <- "SAMPLE_ID"
    write.delim(d, filename, row.names = FALSE)
}

#' @export
find_outliers <- function(mat) {
    outliers <- find_cell_indices(mat, 1, mark_outliers)
    data.frame(
        COMP_ID = outliers[, 1, drop = FALSE],
        SAMPLE_ID = outliers[, 2, drop = FALSE],
        stringsAsFactors = FALSE)
}

#' @export
count_outliers_by_pathway <- function(outliers, pathways, labels) {
    named_pathways <- lapply(pathways, function(pathways) {
        with(pathways, setNames(SUPER_PATHWAY, COMP_ID))
    })
    ds <- Map(function(outliers, pathways, label) {
        as.data.frame(table(pathways[outliers$COMP_ID],
            dnn = "Outliers"), responseName = label)
    }, outliers, named_pathways, labels)
    Reduce(function(x, y) merge(x, y, by = "Outliers"), ds)
}

#' @export
count_outliers_by_run_day <- function(outliers, run_days, labels) {
    named_run_days <- lapply(run_days, function(run_days) {
        with(run_days, setNames(RUN_DAY, SAMPLE_ID))
    })
    ds <- Map(function(outliers, run_days, label) {
        d <- as.data.frame(table(run_days[outliers$SAMPLE_ID]))
        names(d) <- c(sprintf("Run Day (%s)", label), "Outliers")
        d[order(d$Outliers, decreasing = TRUE), , drop = FALSE]
    }, outliers, named_run_days, labels)
    do.call(cbind, ds)
}

#' @export
plot_NAs_per_metabolite <- function(main, filename, mat, width) {
    jpeg(filename, width = width, height = 350)
    hist_NAs(mat, 1, main = main,
        xlab = sprintf("Number of missing samples (N = %d)", ncol(mat)),
        ylab = sprintf("Number of metabolites (N = %d)", nrow(mat)))
    dev.off()
}

#' @export
plot_NAs_per_metabolite_by_pathway <- function(main, filename, mat, pathways, width) {
    jpeg(filename, width = width, height = 1.2 * jpeg_height, quality = jpeg_quality)
    ds <- split.data.frame(mat, pathways)
    par(mfrow = c(length(ds), 1), mar = c(3, 5, 0, 1), oma = c(2, 0, 4, 0), cex = 1)
    for (i in seq_along(ds)) {
        hist_NAs(ds[[i]], 1, breaks = 10, xlim = c(0, ncol(ds[[i]])),
            ylab = names(ds)[i], xlab = "", main = "", xpd = NA)
    }
    title(main, outer = TRUE)
    title(xlab = sprintf("Number of missing samples (N = %d)", ncol(mat)), line = 0, outer = TRUE)
    dev.off()
}

#' @export
show_distribution_of_NAs_by_pathway <- function(mat, pathways) {
    x <- apply(mat, 1, count_NAs)
    y <- do.call(rbind, by(x, pathways, summary))
    z <- apply(100 * y / ncol(mat), 2, function(x) {
        sprintf("%.1f", x)
    })
    rownames(z) <- rownames(y)
    z
}

#' @export
save_NAs_per_metabolite <- function(filename, mat) {
    d <- summarize_NAs(mat, 1)
    names(d)[1] <- "COMP_ID"
    write.delim(d, filename, row.names = FALSE)
}

#' @export
cv <- function(x) {
    sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)
}

#' @export
compare_cv_strategies <- function(mat, run_days, outliers, percentage) {
    set_high_cvs_to_NA <- function(x, percentage) {
        too_big <- x > quantile(x, percentage, na.rm = TRUE)
        x[too_big] <- NA
        x
    }
    mat2 <- mat
    mat2[outliers] <- NA
    cvs1 <- t(columnwise(cv, t(mat), run_days))
    cvs2 <- t(columnwise(cv, t(mat2), run_days))
    cvs3 <- set_high_cvs_to_NA(cvs1, percentage)
    cvs4 <- set_high_cvs_to_NA(cvs2, percentage)
    cv_threshold <- .25
    any2 <- function(x, ...) {
        if (any(x > cv_threshold, na.rm = TRUE)) 1 else 0
    }
    f <- too_high_cv <- function(cvs, f) {
        too_high <- apply(cvs, 1, function(x) {
            f(x, na.rm = TRUE) > cv_threshold
        })
        sum(too_high, na.rm = TRUE)
    }
    data.frame(
        type = c("any", "median", "mean"),
        "1" = c(f(cvs1, any2), f(cvs1, median), f(cvs1, mean)),
        "2" = c(f(cvs2, any2), f(cvs2, median), f(cvs2, mean)),
        "3" = c(f(cvs3, any2), f(cvs3, median), f(cvs3, mean)),
        "4" = c(f(cvs4, any2), f(cvs4, median), f(cvs4, mean)),
        check.names = FALSE)
}

#' @export
read_data <- function(filename, sample_id = "SAMPLE_ID") {
    lines_before_data <- find_lines_to_skip(filename)
    d <- read.delim(filename, na.strings = "",
        skip = lines_before_data, stringsAsFactors = FALSE,
        check.names = FALSE)
    ids <- find_header_names(filename, sample_id)
    d <- use_ids_as_column_names(d, ids)
    d
}

find_lines_to_skip <- function(filename) {
    min(grep("^\t+", readLines(filename), invert = TRUE)) - 1L
}

#' @export
find_header_names <- function(filename, keyword, sep = "\t") {
    contents <- readLines(filename)
    header <- grep(keyword, contents, value = TRUE)
    garbage <- sprintf("^%s+%s%s+", sep, keyword, sep)
    header_names <- unlist(strsplit(sub(garbage, "", header), sep), use.names = FALSE)
    header_names[header_names == ""] <- NA
    header_names
}

use_ids_as_column_names <- function(data, ids) {
    cols <- colnames(data)
    id_columns <- seq(length(cols) - length(ids) + 1L, length(cols))
    cols[id_columns] <- ids
    colnames(data) <- cols
    data
}

#' @export
extract_data_matrix <- function(data, starting_at_column, metabolite_id) {
    m <- as.matrix(data[, seq(starting_at_column, ncol(data))])
    rownames(m) <- data[[metabolite_id]]
    m
}

#' @export
scale_metabolites_to_median_one <- function(mat) {
    metabolite_medians <- apply(mat, 1, median, na.rm = TRUE)
    mat / metabolite_medians
}

#' @export
extract_and_order_metabolites <- function(mat, metabolites) {
    mat[metabolites, , drop = FALSE]
}

#' @export
extract_and_order_samples <- function(mat, samples) {
    mat[, samples, drop = FALSE]
}

#' Check whether objects are identical
#'
#' @param \dots Any number of objects
#' @param objects List of objects to be compared
#' @param key Function to apply to objects before comparison (see
#'     Details)
#' @details If the objects to be compared are referenced by separate
#'     variables, use the \dots argument.  If you have a list of
#'     objects, use the \code{objects} argument.  You are free to use
#'     both the \dots and \code{objects} arguments at the same time.
#'
#'     If \code{key} is \code{NULL}, objects are compared directly.
#'     If \code{key} is a function, it will be applied to the objects
#'     and the return values will be compared.  This makes it possible
#'     to compare objects in terms of arbitrary attributes.
#' @return Returns \code{TRUE} if all objects are identical, otherwise
#'     returns \code{FALSE}.
#' @examples
#' \dontrun{
#' # Do data frames d1, d2, and d3 have the same row names?
#' same(d1, d2, d3, key = rownames)
#' # Equivalent to:
#' same(objects = list(d1, d2, d3), key = rownames)
#' # Or:
#' same(d1, objects = list(d2, d3), key = rownames)
#' }
#' @export
same <- function(..., objects = NULL, key = NULL) {
    xs <- c(list(...), objects)
    if (is.null(key))
        key <- identity
    reference <- key(xs[[1]])
    for (x in xs)
        if (!identical(key(x), reference))
            return(FALSE)
    TRUE
}

#' @export
set_negative_and_zero_to_NA <- function(x) {
    x[x <= 0] <- NA
    x
}

#' @export
read_as_matrix <- function(filename, colClasses = NA) {
    as.matrix(read.delim(filename, check.names = FALSE, colClasses = colClasses))
}

#' @export
count_metabolites_with_high_global_cv <- function(mat) {
    sum(apply(mat, 1, cv) > .25, na.rm = TRUE)
}

#' @export
find_sample_type <- function(filename, sample_id, sample_type) {
    data.frame(
        SAMPLE_ID = find_header_names(filename, sample_id),
        SAMPLE_TYPE = find_header_names(filename, sample_type),
        stringsAsFactors = FALSE)
}

#' @export
count_samples_per_run_day <- function(study_run_days, reference_run_days) {
    f <- function(x, label) {
        d <- as.data.frame(table(x))
        names(d) <- c("Run day", label)
        d
    }
    d1 <- f(study_run_days, "Study")
    d2 <- f(reference_run_days, "Reference")
    cbind(d1, d2)
}

#' @export
count_samples_and_metabolites <- function(matrices) {
    counts <- do.call(cbind, lapply(matrices, function(m) c(ncol(m), nrow(m))))
    dimnames(counts) <- list(c("samples", "metabolites"), names(matrices))
    counts
}

#' @export
count_negative_and_zero_measurements <- function(matrices) {
    counts <- do.call(cbind, lapply(matrices, function(m) {
        x <- as.vector(m)
        c(sum(x < 0, na.rm = TRUE), sum(x == 0, na.rm = TRUE))
    }))
    dimnames(counts) <- list(c("negative", "zero"), names(matrices))
    counts
}

#' @export
`%c%` <- function(f, g) {
    function(x) f(g(x))
}

#' @export
count_measurements_per_run_day <- function(matrices, run_days) {
    stopifnot(same(objects = run_days, key = levels))
    do.call(cbind, unname(Map(function(mat, run_days, label) {
        d <- data.frame(run_days = levels(run_days),
            measurements = tapply(run_days, run_days, length) * nrow(mat))
        names(d) <- c("Run day", label)
        d
    }, matrices, run_days, names(matrices))))
}

#' Apply function to every sample
#'
#' @param mat Matrix with samples as columns and metabolites as rows
#' @param f Function to apply to samples
#' @param \dots Further arguments passed to \code{f}
#' @return If \code{f} returns a vector of length 1,
#'     \code{each_sample} returns a vector with one element per
#'     sample.  Otherwise if \code{f} returns a vector of the same
#'     length as its argument, \code{each_sample} returns a matrix
#'     with one column per sample, that is, the result will have the
#'     same dimensions as \code{mat}.
#' @export
each_sample <- function(mat, f, ...) {
    apply(mat, 2, f, ...)
}

#' Apply function to every metabolite
#'
#' @param mat Matrix with samples as columns and metabolites as rows
#' @param f Function to apply to metabolites
#' @param \dots Further arguments passed to \code{f}
#' @return If \code{f} returns a vector of length 1,
#'     \code{each_metabolite} returns a vector with one element per
#'     metabolite.  Otherwise if \code{f} returns a vector of the same
#'     length as its argument, \code{each_metabolite} returns a matrix
#'     with one row per metabolite, that is, the result will have the
#'     same dimensions as \code{mat}.
#' @export
each_metabolite <- function(mat, f, ...) {
    x <- apply(mat, 1, f, ...)
    if (is.vector(x))
        x
    else
        t(x)
}

#' Apply function to every run day
#'
#' @param mat Matrix with samples as columns and metabolites as rows
#' @param run_days Factor of length \code{ncol(mat)} used to split
#'     samples by run day
#' @param f Function to apply to run day subsets of \code{mat}
#' @param \dots Further arguments passed to \code{f}
#' @return A list with one element per run day.  Elements are named
#'     using the corresponding level of \code{run_days}.
#' @export
each_run_day <- function(mat, run_days, f, ...) {
    if (is.vector(mat))
        mat <- matrix(mat, nrow = 1)
    lapply(lapply(split.data.frame(t(mat), run_days), t), f, ...)
}

#' Apply function to every pathway
#'
#' @param mat Matrix with samples as columns and metabolites as rows
#' @param pathways Factor of length \code{nrow(mat)} used to split
#'     metabolites by pathway
#' @param f Function to apply to pathway subsets of \code{mat}
#' @param \dots Further arguments passed to \code{f}
#' @return A list with one element per pathway.  Elements are named
#'     using the corresponding level of \code{pathways}.
#' @export
each_pathway <- function(mat, pathways, f, ...) {
    if (is.vector(mat))
        mat <- matrix(mat, ncol = 1)
    lapply(split.data.frame(mat, pathways), f, ...)
}
