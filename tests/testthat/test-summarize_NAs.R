context("summarize_NAs")

mat <- matrix(1:4, nrow = 2)
colnames(mat) <- c("x", "y")
mat[1, 2] <- NA

test_that("by column", {
    actual <- summarize_NAs(mat, by = 2)
    expected <- data.frame(
        id = colnames(mat),
        missings = c(0L, 1L),
        missingness_rate = c(0, .5),
        stringsAsFactors = FALSE)
    expect_equal(actual, expected)
})

test_that("by row", {
    actual <- summarize_NAs(mat, by = 1)
    expected <- data.frame(
        id = as.character(seq_len(nrow(mat))),
        missings = c(1L, 0L),
        missingness_rate = c(.5, 0),
        stringsAsFactors = FALSE)
    expect_equal(actual, expected)
})
