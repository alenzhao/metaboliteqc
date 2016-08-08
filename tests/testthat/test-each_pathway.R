context("each_pathway")

mat <- rbind(
    c(1, 2, 3, 4),
    c(1, 2, 3, 4),
    c(1, 2, 3, 4))
colnames(mat) <- paste0("p", seq_len(ncol(mat)))
rownames(mat) <- paste0("m", seq_len(nrow(mat)))

run_days <- c("rd1", "rd1", "rd2", "rd2")
pathways <- c("pw1", "pw1", "pw2")

test_that("happy path", {
    actual <- each_pathway(mat, pathways, length %c% as.vector)
    expected <- list(pw1 = 8, pw2 = 4)
    expect_equal(actual, expected)
})

test_that("using ...", {
    actual <- each_pathway(mat, pathways, `<`, 2.5)
    pw1 <- matrix(c(TRUE, TRUE, FALSE, FALSE), nrow = 2, ncol = ncol(mat),
        byrow = TRUE, dimnames = list(c("m1", "m2"), colnames(mat)))
    pw2 <- matrix(c(TRUE, TRUE, FALSE, FALSE), nrow = 1, ncol = ncol(mat),
        byrow = TRUE, dimnames = list("m3", colnames(mat)))
    expected <- list(pw1 = pw1, pw2 = pw2)
    expect_equal(actual, expected)
})

test_that("a vector is interpreted as a matrix with one sample and multiple metabolites", {
    x <- mat[, 1, drop = FALSE]
    actual <- each_pathway(as.vector(x), pathways, `<`, 2.5)
    expected <- lapply(each_pathway(x, pathways, `<`, 2.5), unname)
    expect_equal(actual, expected)
})
