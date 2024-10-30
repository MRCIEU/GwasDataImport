af_conflicts_function<-function(out=NULL)
{
	af_conflicts<-CheckSumStats::flag_af_conflicts(target_dat=out)
	af_conflicts_test<-"test not possible"
	if(af_conflicts$number_of_snps>10 & af_conflicts$proportion_conflicts>=0.10){
		af_conflicts_test<-"fail"
	}		
	if(af_conflicts$number_of_snps>10 & af_conflicts$proportion_conflicts<0.10){
		af_conflicts_test<-"pass"
	}		
	# if(af_conflicts_test=="fail") warning("allele frequency conflicts identified. It looks like you have incorrectly specified the effect allele frequency column.")
	af_conflicts<-c(af_conflicts,"test"=af_conflicts_test)
	Plot1<-CheckSumStats::make_plot_maf(ref_1000G=c("AFR","AMR","EAS","EUR","SAS","ALL"),target_dat=out)
	# af_conflicts<-c(af_conflicts,"plot"=Plot1)
	return(list(af_conflicts,Plot1))
}

# stop("allele frequency conflicts identified. It looks like you have incorrectly specified the effect allele frequency column.")
# , efo_id=NULL,trait=NULL,ignore_conflict=FALSE,map_association_to_study=FALSE
gc_conflicts_function<-function(out=NULL,metadata=NULL)
{

	gwas_catalog<-CheckSumStats::gwas_catalog_hits(trait=metadata$trait)
	if(length(gwas_catalog) == 1)
	{
		if(gwas_catalog == "no results found")
		{
			gwas_catalog<-CheckSumStats::gwas_catalog_hits(efo_id=metadata$ontology)
		}
	}
	gc_dat<-CheckSumStats::compare_effect_to_gwascatalog2(dat=out,gwas_catalog=gwas_catalog,trait=metadata$trait)
	# if(gc_dat == "no results found"){
	# 	gc_dat<-CheckSumStats::compare_effect_to_gwascatalog2(dat=out,efo_id=metadata$ontology)
	# }

	gc_conflicts<-CheckSumStats::flag_gc_conflicts2(gc_dat=gc_dat)	
	
	N_snps<-gc_conflicts$effect_size_conflicts$n_snps
	
	n_high<-unlist(gc_conflicts$effect_size_conflicts["high conflict"])
	n_moderate<-unlist(gc_conflicts$effect_size_conflicts["moderate conflict"])
	N_conflicts<-n_high + n_moderate
	proportion_conflicts<-N_conflicts/N_snps

	gc_conflicts_test<-"Number of SNPs may be too low for reliable test. Interpret with caution"
	if(N_snps>10){
		if(proportion_conflicts >=0.8){ 
			gc_conflicts_test<-"strong_fail"
			# warning("Very strong effect size conflicts with GWAS catalog identified. The effect allele column looks wrong")
		}
		if(proportion_conflicts<0.8 & proportion_conflicts>0.3){
			gc_conflicts_test<-"moderate_fail"
			# warning("Effect size conflicts with GWAS catalog identified. It looks like you may have incorrectly specified the effect allele column")
		}
		if(proportion_conflicts<0.3){
			gc_conflicts_test<-"pass"
		}
	}
	gc_conflicts<-c(gc_conflicts,"test"=gc_conflicts_test)
	Plot2<-CheckSumStats::make_plot_gwas_catalog(dat=out,gwas_catalog=gwas_catalog,gc_dat=gc_dat)
	return(list(gc_conflicts,Plot2))
}

reported_gwas_hits_in_gwascatalog <- function(out=NULL,clump=TRUE,eaf_test=NULL,metadata=NULL)
{	
	Pop<-infer_ancestry(target_dat=out)
	l2<-unlist(Pop)
	Pop<-names(l2[which(l2==max(l2))])
	dat<-out[which(out$p<5e-8),]
	# l<-unlist(eaf_test)
	# if(l[names(l)=="test"] =="fail") warning("eaf test failed, which might undermine the infer_ancestry test for ld_clumping")
	if(nrow(dat)==0) return("no GWAS hits at 5e-8 threshold")

	Clump<-ieugwasr::ld_clump(dat=dplyr::tibble(rsid=dat$rsid, pval=dat$p),clump_r2 = 0.001,clump_p=5e-8,pop=Pop) #it seems like there is a bug such that when none of the variants are present in the reference panel this function returns a Server code:503.  
	gwas_hits<-Clump$rsid
	if(length(Clump$rsid)<5)
	{
		gwas_hits<-dat$rsid
	}

	gc_list<-find_hits_in_gwas_catalog(gwas_hits=gwas_hits,trait=metadata$trait,distance_threshold=50000) 
	#if no results are found on reported trait, then use the user provided ontology to find associations
	if(!is(gc_list, "list")){
		if(gc_list == "no results found")
		{
			metadata$ontology<-gsub(":","_",metadata$ontology) #currently there is a bug in the gwasrapidd code that does not allow EFOs with semicolons. The semicolon must be replaced with an underscore. I have writtent to the developers to see if this can be fixed.  
			gc_list<-find_hits_in_gwas_catalog(gwas_hits=gwas_hits,efo_id=metadata$ontology,distance_threshold=50000) 
		}
	}

	if(!is(gc_list, "list")){
		if(gc_list == "no results found") return(gc_list)
	}

	# metadata$ontology

	N_in_gc<-length(gc_list$in_gc)
	N_not_in_gc<-length(gc_list$not_in_gc )
	N<-N_not_in_gc+N_in_gc
	test_hits_gc<-"test not possible"
	if(N_not_in_gc/N>0.3)
	{
		test_hits_gc<-fail
	}
	if(N_not_in_gc/N<0.3)
	{
		test_hits_gc<-"pass"
	}
	false_positive_hits<-c(gc_list,"test"=test_hits_gc)

	return(false_positive_hits)
}
