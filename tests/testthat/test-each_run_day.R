context("each_run_day")

mat <- rbind(
    c(1, 2, 3, 4),
    c(1, 2, 3, 4),
    c(1, 2, 3, 4))
colnames(mat) <- paste0("p", seq_len(ncol(mat)))
rownames(mat) <- paste0("m", seq_len(nrow(mat)))

run_days <- c("rd1", "rd1", "rd2", "rd2")
pathways <- c("pw1", "pw1", "pw2")

test_that("happy path", {
    actual <- each_run_day(mat, run_days, mean %c% as.vector)
    expected <- list(rd1 = 1.5, rd2 = 3.5)
    expect_equal(actual, expected)
})

test_that("using ...", {
    actual <- each_run_day(mat, run_days, `<`, 2.5)
    rd1 <- matrix(TRUE, nrow = nrow(mat), ncol = 2,
        dimnames = list(rownames(mat), c("p1", "p2")))
    rd2 <- matrix(FALSE, nrow = nrow(mat), ncol = 2,
        dimnames = list(rownames(mat), c("p3", "p4")))
    expected <- list(rd1 = rd1, rd2 = rd2)
    expect_equal(actual, expected)
})

test_that("a vector is interpreted as a matrix with one metabolite and multiple samples", {
    x <- mat[1, , drop = FALSE]
    actual <- each_run_day(as.vector(x), run_days, `<`, 2.5)
    expected <- lapply(each_run_day(x, run_days, `<`, 2.5), unname)
    expect_equal(actual, expected)
})
