library(testthat)
library(MTFM)
set.seed(42)
M <- matrix(sample(c(0,1), 100, replace=TRUE, prob=c(0.8,0.2)),10)
o11 <- modoverlap(c(1,0,0,1),c(1,0,0,0))
o12 <- modoverlap(c(1,0,0,1),c(1,0,0,2))
o23 <- modoverlap(c(0,0,1,1),c(1,1,0,0))
test_that("modoverlap produces correct results", {
    expect_equal(o11,1)
    expect_equal(o12,1)
    expect_equal(o23,0)
})

