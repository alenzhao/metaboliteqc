context("same")

test_that("1 == 1", {
    expect_true(same(1, 1))
})

test_that("1 != 2", {
    expect_false(same(1, 2))
})

test_that("with key argument", {
    expect_true(same(c(1, 2), c(3, 4), key = length))
})

test_that("NULL == NULL", {
    expect_true(same(NULL, NULL))
})

test_that("NA != NULL", {
    expect_false(same(NA, NULL))
})

test_that("integer(0) != NULL", {
    expect_false(same(integer(0), NULL))
})

test_that("same(x, y) <==> same(y, x)", {
    x <- 1:2
    y <- rev(x)
    expect_false(same(x, y))
    expect_false(same(y, x))
})

test_that("with objects argument", {
    expect_true(same(-1, objects = list(1, -1), key = abs))
})
