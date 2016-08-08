context("each_run_day and each_pathway")

mat <- rbind(
    c(1, 2, 3, 4),
    c(1, 2, 3, 4),
    c(1, 2, 3, 4))
colnames(mat) <- paste0("p", seq_len(ncol(mat)))
rownames(mat) <- paste0("m", seq_len(nrow(mat)))

run_days <- c("rd1", "rd1", "rd2", "rd2")
pathways <- c("pw1", "pw1", "pw2")

test_that("each_run_day", {
    actual <- each_run_day(mat, run_days, mean %c% as.vector)
    expected <- list(rd1 = 1.5, rd2 = 3.5)
    expect_equal(actual, expected)
})

test_that("each_pathway", {
    actual <- each_pathway(mat, pathways, length %c% as.vector)
    expected <- list(pw1 = 8, pw2 = 4)
    expect_equal(actual, expected)
})

test_that("each_run_day using ...", {
    actual <- each_run_day(mat, run_days, `<`, 2.5)
    rd1 <- matrix(TRUE, nrow = nrow(mat), ncol = 2,
        dimnames = list(rownames(mat), c("p1", "p2")))
    rd2 <- matrix(FALSE, nrow = nrow(mat), ncol = 2,
        dimnames = list(rownames(mat), c("p3", "p4")))
    expected <- list(rd1 = rd1, rd2 = rd2)
    expect_equal(actual, expected)
})

test_that("each_pathway using ...", {
    actual <- each_pathway(mat, pathways, `<`, 2.5)
    pw1 <- matrix(c(TRUE, TRUE, FALSE, FALSE), nrow = 2, ncol = ncol(mat),
        byrow = TRUE, dimnames = list(c("m1", "m2"), colnames(mat)))
    pw2 <- matrix(c(TRUE, TRUE, FALSE, FALSE), nrow = 1, ncol = ncol(mat),
        byrow = TRUE, dimnames = list("m3", colnames(mat)))
    expected <- list(pw1 = pw1, pw2 = pw2)
    expect_equal(actual, expected)
})
