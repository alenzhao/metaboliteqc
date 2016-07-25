context("count_below_percentile")

mat <- as.matrix(read.table(text = "
            id-1 id-2
pathway-1.a    1    2
pathway-1.b    1   NA
pathway-2.a   10   20
pathway-2.b   10   NA
", check.names = FALSE))

pathway <- sub("..$", "", rownames(mat))

test_that("all samples are <= 100 % percentile for all categories", {
    percentile <- 100
    actual <- count_below_percentile(mat, by = pathway, mean, percentile, na.rm = TRUE)
    expected <- setNames(c(2, 2), colnames(mat))
    expect_equal(actual, expected)
})

test_that("one sample is <= 50 % percentile for all categories, the other for none", {
    percentile <- 50
    actual <- count_below_percentile(mat, by = pathway, mean, percentile, na.rm = TRUE)
    expected <- setNames(c(2, 0), colnames(mat))
    expect_equal(actual, expected)
})
