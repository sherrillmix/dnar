context("Rainbow functions")

test_that("Test rainbow.lab",{
	expect_equal(length(rainbow.lab(100)),100)
	expect_equal(sum(grepl('#[0-9A-F]{8}',rainbow.lab(100))),100)
})

