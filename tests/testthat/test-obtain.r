context("EBI upload")
library(GwasDataImport)

test_that("EBI upload", {
	skip_on_travis()
	ebi_id <- "GCST005522"
	x <- EbiDataset$new(
		ebi_id = ebi_id, 
		ftp_path = get_ftp_path(ebi_id),
		igd_id = "ebi-a-EBITEST",
		wd = paste0("uploads/", ebi_id)
	)
	# x$pipeline()

	# expect_true(file.exists(x$gwas_out))
	# expect_true(x$nsnp > x$nsnp_read * 0.8)
	# expect_true(x$metadata$ncase == 1886)
	# expect_true(x$metadata$ncontrol == 10421)
	# expect_true(x$metadata$population == "European")
	# expect_true(x$metadata$trait == "Narcolepsy")
	# expect_true(x$metadata$unit == "logOR")
	# x$api_gwasdata_delete()
})







