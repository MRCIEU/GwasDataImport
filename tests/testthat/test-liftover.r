context("Liftover")
library(GwasDataImport)

test_that("liftover", {
  skip("Tests for GH to run interactively")
	filename <- system.file(package="GwasDataImport", "extdata/pos_0002.txt.gz")
	dat <- data.table::fread(filename)
	a <- liftover_gwas(dat, to=38, chr_col="CHROM", pos_col="POS", snp_col="RSID", ea_col="ALT", oa_col="REF")
#	b <- liftover_gwas2(dat, to=38, chr_col="CHROM", pos_col="POS", snp_col="RSID", ea_col="ALT", oa_col="REF")
	expect_true(nrow(a) > nrow(dat)*0.8)
	a <- liftover_gwas(dat, to=36, chr_col="CHROM", pos_col="POS", snp_col="RSID", ea_col="ALT", oa_col="REF")
	expect_true(nrow(a) > nrow(dat)*0.8)
	a <- determine_build_biomart(dat$RSID, dat$CHROM, dat$POS)

	b <- liftover_gwas(a, to=37, chr_col="CHROM", pos_col="POS", snp_col="RSID", ea_col="ALT", oa_col="REF")
})

# microbenchmark::microbenchmark(liftover_gwas(dat, to=38, chr_col="CHROM", pos_col="POS", snp_col="RSID", ea_col="ALT", oa_col="REF"), liftover_gwas2(dat, to=38, chr_col="CHROM", pos_col="POS", snp_col="RSID", ea_col="ALT", oa_col="REF"), times=10)


# dg <- GRanges(seqnames=paste0("chr",dat$CHROM), ranges=IRanges(start=dat$POS, end=dat$POS), mcols=DataFrame(ind=1:nrow(dat))
# dg


# dat2 <- data.table::fread("~/Downloads/Stroke_Bothsex_eur_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz")
# microbenchmark::microbenchmark(a <- liftover_gwas(dat2, to=37, chr_col="#CHR", pos_col="POS", snp_col="rsid", ea_col="ALT", oa_col="REF"), times=1)
# microbenchmark::microbenchmark(b <- liftover_gwas2(dat2, to=37, chr_col="#CHR", pos_col="POS", snp_col="rsid", ea_col="ALT", oa_col="REF"), times=1)

