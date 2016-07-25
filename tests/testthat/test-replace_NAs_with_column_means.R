context("replace_NAs_with_column_means")

mat <- cbind(c(1, 2, NA), c(NA, 2, 3))
rownames(mat) <- c("a", "b", "c")
colnames(mat) <- c("1", "2")

test_that("replace_NAs_with_column_means", {
    actual <- replace_NAs_with_column_means(mat)
    expected <- cbind(c(1, 2, 1.5), c(2.5, 2, 3))
    dimnames(expected) <- dimnames(mat)
    expect_equal(actual, expected)
})

test_that("replace_NAs_with_row_means", {
    actual <- replace_NAs_with_row_means(mat)
    expected <- cbind(c(1, 2, 3), c(1, 2, 3))
    dimnames(expected) <- dimnames(mat)
    expect_equal(actual, expected)
})

test_that("all NA stays NA", {
    mat[1, ] <- NA
    actual <- replace_NAs_with_row_means(mat)
    expected <- cbind(c(NA, 2, 3), c(NA, 2, 3))
    dimnames(expected) <- dimnames(mat)
    expect_equal(actual, expected)
})
