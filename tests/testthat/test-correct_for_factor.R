context("correct_for_factor")

set.seed(12345)
n <- 5
mat <- 100 * rbind(rnorm(n), rnorm(n), rnorm(n))
rownames(mat) <- c("a", "b", "c")
colnames(mat) <- seq_len(ncol(mat))
variable <- c(1, 1, 1, 2, 2)
actual <- correct_for_factor(mat, variable)

change_first_row_to <- function(mat, row) {
    mat[1, ] <- row
    mat
}

expect_nothing_but_NAs_in_first_row <- function(actual, number_of_columns) {
    expect_equal(unname(actual[1, ]), rep(NA_real_, number_of_columns))
}

test_that("dimnames", {
    expect_equal(dimnames(actual), dimnames(mat))
})

test_that("the residuals for a given metabolite average to 0", {
    expect_lt(max(abs(apply(actual, 1, mean))), 1e-10)
})

test_that("the residuals when averaged over samples are different from 0", {
    expect_gt(min(abs(apply(actual, 2, mean))), 10)
})

test_that("an NA-only metabolite gets NAs as residuals", {
    mat2 <- change_first_row_to(mat, rep(NA, ncol(mat)))
    actual <- correct_for_factor(mat2, variable)
    expect_nothing_but_NAs_in_first_row(actual, ncol(mat2))
})

test_that("a metabolite with NA variance gets NAs as residuals", {
    mat2 <- change_first_row_to(mat, c(1, rep(NA, ncol(mat) - 1)))
    actual <- correct_for_factor(mat2, variable)
    expect_nothing_but_NAs_in_first_row(actual, ncol(mat2))
})

test_that("a metabolite with zero variance gets NAs as residuals", {
    mat2 <- change_first_row_to(mat, rep(1, ncol(mat)))
    actual <- correct_for_factor(mat2, variable)
    expect_nothing_but_NAs_in_first_row(actual, ncol(mat2))
})
