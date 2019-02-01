library(testthat)
library(MTFM)

context("modoverlap")
o11 <- MFM:::modoverlap(c(1,0,0,1),c(1,0,0,0))
o12 <- MFM:::modoverlap(c(1,0,0,1),c(1,0,0,2))
o23 <- MFM:::modoverlap(c(0,0,1,1),c(1,1,0,0))
test_that("modoverlap produces correct results", {
    expect_equal(o11,1)
    expect_equal(o12,1)
    expect_equal(o23,0)
})

context("stroverlap")
o11 <- MFM:::stroverlap(c(1,2,3),c(4,5,6,7))
o12 <- MFM:::stroverlap(c(1,2,3),c(5,6,7,1))
o23 <- MFM:::stroverlap(numeric(0),c(1,0))
test_that("modoverlap produces correct results", {
    expect_equal(o11,0)
    expect_equal(o12,1)
    expect_equal(o23,0)
})
