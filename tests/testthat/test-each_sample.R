context("each_sample and each_metabolite")

mat <- rbind(
    c(-1, +1, +1, +1),
    c(-1, -1, +1, +1),
    c(-1, -1, -1, +1))
rownames(mat) <- paste0("m", seq_len(nrow(mat)))
colnames(mat) <- paste0("p", seq_len(ncol(mat)))

test_that("each_sample with function returning non-scalar", {
    actual <- each_sample(mat, abs)
    expected <- abs(mat)
    expect_equal(actual, expected)
})

test_that("each_sample with function returning scalar", {
    actual <- each_sample(mat, sum)
    expected <- setNames(c(-3, -1, 1, 3), colnames(mat))
    expect_equal(actual, expected)
})

test_that("each_metabolite with function returning non-scalar", {
    actual <- each_metabolite(mat, abs)
    expected <- each_sample(mat, abs)
    expect_equal(actual, expected)
})

test_that("each_sample with function returning scalar", {
    actual <- each_metabolite(mat, sum)
    expected <- setNames(c(2, 0, -2), rownames(mat))
    expect_equal(actual, expected)
})

test_that("each_sample and each_metabolite", {
    actual <- each_sample(each_metabolite(mat, `<`, 0), sum)
    expected <- setNames(c(3, 2, 1, 0), colnames(mat))
    expect_equal(actual, expected)
})
