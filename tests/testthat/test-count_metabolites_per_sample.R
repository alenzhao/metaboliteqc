context("count_metabolites_per_sample")

mat <- rbind(
    c(1, 2, 3, 4),  # metabolite 1
    c(1, 1, 2, 3),  # metabolite 2
    c(1, 1, 1, 2))  # metabolite 3
colnames(mat) <- paste0("p", seq_len(ncol(mat)))

test_that("count_metabolites_per_sample", {
    actual <- count_metabolites_per_sample(mat, `==`, 1)
    expected <- setNames(c(3, 2, 1, 0), colnames(mat))
    expect_equal(actual, expected)
})
