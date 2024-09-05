test_that("infoset function works as expected", {
  skip_on_cran()  # Skip this test on CRAN

  # Your test code here
  gross.ret<-as.data.frame(lapply(sample.data, g_ret))
  infoset(gross.ret$ETF_1)

})
