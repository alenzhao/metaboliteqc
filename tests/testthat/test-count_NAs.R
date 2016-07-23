context("count_NAs")

f <- count_NAs

test_that("vectors", {
    expect_equal(f(c(1, NA)), 1)
})

test_that("list", {
    expect_equal(f(list(1, NA)), 1)
})

test_that("matrix", {
    expect_equal(f(rbind(1, NA)), 1)
})

test_that("data frame", {
    expect_equal(f(data.frame(1, NA)), 1)
})
