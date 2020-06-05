context("Obtain data")
library(GwasDataImport)

test_that("data obtain pipeline", {

	ebi_id <- "GCST005522"
	x <- EbiDataset$new(
		ebi_id = ebi_id, 
		ftp_path=get_ftp_path(ebi_id),
		igd_id = "ebi-a-EBITEST",
		wd = paste0("uploads/", ebi_id)
	)
	x$pipeline()

	expect_true(file.exists(x$metadata_file))
	expect_true(file.exists(x$datainfo_file))
	expect_true(x$nsnp == 93613)
	expect_true(x$metadata$ncase == 1886)
	expect_true(x$metadata$ncontrol == 10421)
	expect_true(x$metadata$population == "European")
	expect_true(x$metadata$trait == "Narcolepsy")
	expect_true(x$metadata$unit == "logOR")
	x$api_gwasdata_delete()

})




