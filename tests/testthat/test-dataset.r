context("Standard upload")
library(GwasDataImport)

test_that("Standard upload", {

	filename <- system.file(package="GwasDataImport", "extdata/pos_0002.txt.gz")
	id <- paste0("ieu-b-testing_", as.numeric(Sys.time())) %>% gsub("\\.", "_", .)
	x <- Dataset$new(filename=filename, igd_id=id)
	x$determine_columns(list(chr_col=1, snp_col=2, pos_col=3, oa_col=4, ea_col=5, eaf_col=6, beta_col=7, se_col=8, pval_col=9))
	x$format_dataset()
	x$collect_metadata(list(trait="hello", build="HG19/GRCh37", category="NA", subcategory="NA", group_name="public", population="European", sex="Males", note='asdasd', ontology="EFO_1234;EFO_2345"))
	x$api_metadata_upload()
	x$api_gwasdata_upload()
	# x$api_report()
	# x$api_gwas_release()
	x$api_gwasdata_delete()


})
