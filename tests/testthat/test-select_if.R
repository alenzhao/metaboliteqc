context("select_if")

mat <- rbind(1:3, 4:6)
rownames(mat) <- c("a", "b")
colnames(mat) <- 1:3

test_that("select rows", {
    mean_greater_than_3 <- function(x) mean(x) > 3
    actual <- select_if(mean_greater_than_3, mat, 1)
    expected <- mat[2, , drop = FALSE]
    expect_equal(actual, expected)
})

test_that("select columns", {
    divisible_by_3 <- function(x) all(x %% 3L == 0L)
    actual <- select_if(divisible_by_3, mat, 2)
    expected <- mat[, 3, drop = FALSE]
    expect_equal(actual, expected)
})

test_that("we get a matrix even if no row is selected", {
    under_no_circumstances <- function(x) FALSE
    actual <- select_if(under_no_circumstances, mat, 1)
    expected <- matrix(integer(0), nrow = 0, ncol = ncol(mat),
        dimnames = list(NULL, colnames(mat)))
    expect_equal(actual, expected)
})

test_that("we get a matrix even if no column is selected", {
    under_no_circumstances <- function(x) FALSE
    actual <- select_if(under_no_circumstances, mat, 2)
    expected <- matrix(integer(0), ncol = 0, nrow = nrow(mat),
        dimnames = list(rownames(mat), NULL))
    expect_equal(actual, expected)
})

test_that("using '...'", {
    all_is_precious <- function(x, precious) {
        all(x %in% precious)
    }
    actual <- select_if(all_is_precious, mat, 1, 1:4)
})
