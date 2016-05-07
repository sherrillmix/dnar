#avoid confusing subprocesses in R CMD check
Sys.setenv("R_TESTS" = "")
library(testthat)
library(dnar)
test_check("dnar")

