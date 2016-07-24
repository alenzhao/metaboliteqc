context("columnwise")

mat <- as.matrix(read.table(text = "
            id-1 id-2 id-3
pathway-1.a    1    2    3
pathway-1.b    1    2   NA
pathway-2.a   10   20   30
pathway-2.b   10   20   NA
", check.names = FALSE))

pathway <- sub("..$", "", rownames(mat))

test_that("function with scalar result", {
    actual <- columnwise(mean, mat, by = pathway, na.rm = TRUE)
    expected <- matrix(c(c(1, 2, 3), c(10, 20, 30)), nrow = 2, byrow = TRUE,
        dimnames = list(c("pathway-1", "pathway-2"), c("id-1", "id-2", "id-3")))
    expect_equal(actual, expected)
})

test_that("function with non-scalar result", {
    expect_equal(columnwise(identity, mat, by = pathway), mat)
})
