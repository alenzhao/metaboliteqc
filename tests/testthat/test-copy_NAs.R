context("copy_NAs")

matrix1 <- matrix(1:12, nrow = 3)
colnames(matrix1) <- c("w", "x", "y", "z")
matrix2 <- matrix1
matrix2[1, 1] <- NA

test_that("happy path", {
    m <- copy_NAs(from = matrix2, to = matrix1)
    expect_true(is.na(m[1, 1]))
    expect_true(all(m == matrix1, na.rm = TRUE))
})

test_that("matrices must have the same dimensions", {
    expect_error(copy_NAs(from = matrix2, to = matrix2[, 1]),
        "Arguments FROM and TO must have the same dimensions.")
})

test_that("matrices must have the same dimnames", {
    expect_error(copy_NAs(from = matrix2, to = unname(matrix1)),
        "Arguments FROM and TO must have the same dimnames.")
})
