context("find_cell_indices")

mat <- as.matrix(read.table(text = "
   1  2  3
a -1  4  5
b  2  3 NA
", check.names = FALSE))

is_maximum_or_special <- function(x, special = NULL) {
    x == max(x, na.rm = TRUE) | x %in% special
}

expect_pairs <- function(mat, ...) {
    pairs <- list(...)
    for (pair in pairs) {
        found <- FALSE
        for (row in seq_len(nrow(mat))) {
            if (identical(mat[row, ], pair)) {
                found <- TRUE
                break
            }
        }
        if (!found) {
            fail(sprintf("Pair not found: (\"%s\", \"%s\")", pair[1], pair[2]))
        }
    }
    expect_true(TRUE)
}

test_that("by row", {
    actual <- find_cell_indices(mat, 1, is_maximum_or_special)
    expect_pairs(actual, c("a", "3"), c("b", "2"))
})

test_that("by column", {
    actual <- find_cell_indices(mat, 2, is_maximum_or_special)
    expect_pairs(actual, c("b", "1"), c("a", "2"), c("a", "3"))
})

test_that("using '...'", {
    actual <- find_cell_indices(mat, 1, is_maximum_or_special, special = "-1")
    expect_pairs(actual, c("a", "1"), c("a", "3"), c("b", "2"))
})
