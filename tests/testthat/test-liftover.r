context("Liftover")
library(GwasDataImport)

test_that("liftover", {
	filename <- system.file(package="GwasDataImport", "extdata/pos_0002.txt.gz")
	dat <- data.table::fread(filename)
	a <- liftover_gwas(dat, to=38, chr_col="CHROM", pos_col="POS", snp_col="RSID", ea_col="ALT", oa_col="REF")
	expect_true(nrow(a) > nrow(dat)*0.8)
	a <- liftover_gwas(dat, to=36, chr_col="CHROM", pos_col="POS", snp_col="RSID", ea_col="ALT", oa_col="REF")
	expect_true(nrow(a) > nrow(dat)*0.8)
	a <- determine_build_biomart(dat$RSID, dat$CHROM, dat$POS)
})

