BiocManager::install("sevenbridges")
BiocManager::install("maftools")
library(tidyr)
library(maftools)
library(sevenbridges)

args = commandArgs(trailingOnly = TRUE)
input_txt <- args[1]
input_genes <- args[2]
threshold <- args[3]
output_txt <- sub(".snpEff.vcf.gz.hg38_multianno.txt",".snpEff.hg38_multianno.filt.txt", basename(input_txt))
print(output_txt)

#Convert Annovar fields to VCF
testmaf <- as.data.frame(annovarToMaf(input_txt, refBuild = "hg38"))
colnames(testmaf) <- make.unique(names(testmaf))

#Calculate Allele Fraction
testmaf <-separate(testmaf, V190 ,sep=":", c("genotype","allelic_depth","read_depth2","genotype_quality","phred_likelihoods"))
testmaf <- separate(testmaf, allelic_depth, sep= "," , c("refnum", "altnum"), remove = F)
testmaf$allele_fraction <- (as.numeric(testmaf$altnum)/(as.numeric(testmaf$refnum)+as.numeric(testmaf$altnum)))
testmaf <- cbind(testmaf[7], testmaf[2:6], testmaf[191], testmaf[195], testmaf[202:203], testmaf[205], testmaf[8:14], testmaf[17:47], testmaf[51:185], testmaf[196], testmaf[198:201])
colnames(testmaf) <- c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele2", "rsID", "FILTER", "read_depth", "quality", "allele_fraction", "Variant_Classification", "tx", "exon", "txChange", "aaChange", "Variant_Type", "Func.refGene", "ExonicFunc.refGene", "AAChange.refgene", "InterVar_automated", "PVS1", "PS1", "PS2", "PS3", "PS4", "PM1", "PM2", "PM3", "PM4", "PM5", "PM6", "PP1", "PP2", "PP3", "PP4", "PP5", "BA1", "BS1", "BS2", "BS3", "BS4", "BP1", "BP2", "BP3", "BP4", "BP5", "BP6", "BP7", "nci60", "Kaviar_AF", "Kaviar_AC", "Kaviar_AN", "HRC_AF", "HRC_AC", "HRC_AN", "HRC_non1000G_AF", "HRC_non1000G_AC", "HRC_non1000G_AN", "gnomad_exome_AF", "gnomad_exome_AF_popmax", "gnomad_exome_AF_male", "gnomad_exome_AF_female", "gnomad_exome_AF_raw", "gnomad_exome_AF_afr", "gnomad_exome_AF_sas", "gnomad_exome_AF_amr", "gnomad_exome_AF_eas", "gnomad_exome_AF_nfe", "gnomad_exome_AF_fin", "gnomad_exome_AF_asj", "gnomad_exome_AF_oth", "gnomad_exome_non_topmed_AF_popmax", "gnomad_exome_non_neuro_AF_popmax", "gnomad_exome_non_cancer_AF_popmax", "gnomad_exome_controls_AF_popmax", "gnomad_genome_AF", "gnomad_genome_AF_popmax", "gnomad_genome_AF_male", "gnomad_genome_AF_female", "gnomad_genome_AF_raw", "gnomad_genome_AF_afr", "gnomad_genome_AF_sas", "gnomad_genome_AF_amr", "gnomad_genome_AF_eas", "gnomad_genome_AF_nfe", "gnomad_genome_AF_fin", "gnomad_genome_AF_asj", "gnomad_genome_AF_oth", "gnomad_genome_non_topmed_AF_popmax", "gnomad_genome_non_neuro_AF_popmax", "gnomad_genome_non_cancer_AF_popmax", "gnomad_genome_controls_AF_popmax", "SIFT_score", "SIFT_converted_rankscore", "SIFT_pred", "Polyphen2_HDIV_score", "Polyphen2_HDIV_rankscore", "Polyphen2_HDIV_pred", "Polyphen2_HVAR_score", "Polyphen2_HVAR_rankscore", "Polyphen2_HVAR_pred", "LRT_score", "LRT_converted_rankscore", "LRT_pred", "MutationTaster_score", "MutationTaster_converted_rankscore", "MutationTaster_pred", "MutationAssessor_score", "MutationAssessor_score_rankscore", "MutationAssessor_pred", "FATHMM_score", "FATHMM_converted_rankscore", "FATHMM_pred", "PROVEAN_score", "PROVEAN_converted_rankscore", "PROVEAN_pred", "VEST3_score", "VEST3_rankscore", "MetaSVM_score", "MetaSVM_rankscore", "MetaSVM_pred", "MetaLR_score", "MetaLR_rankscore", "MetaLR_pred", "M-CAP_score", "M-CAP_rankscore", "M-CAP_pred", "REVEL_score", "REVEL_rankscore", "MutPred_score", "MutPred_rankscore", "CADD_raw", "CADD_raw_rankscore", "CADD_phred", "DANN_score", "DANN_rankscore", "fathmm-MKL_coding_score", "fathmm-MKL_coding_rankscore", "fathmm-MKL_coding_pred", "Eigen_coding_or_noncoding", "Eigen-raw", "Eigen-PC-raw", "GenoCanyon_score", "GenoCanyon_score_rankscore", "integrated_fitCons_score", "integrated_fitCons_score_rankscore", "integrated_confidence_value", "GERP++_RS", "GERP++_RS_rankscore", "phyloP100way_vertebrate", "phyloP100way_vertebrate_rankscore", "phyloP20way_mammalian", "phyloP20way_mammalian_rankscore", "phastCons100way_vertebrate", "phastCons100way_vertebrate_rankscore", "phastCons20way_mammalian", "phastCons20way_mammalian_rankscore", "SiPhy_29way_logOdds", "SiPhy_29way_logOdds_rankscore", "Interpro_domain", "GTEx_V6p_gene", "GTEx_V6p_tissue", "avsnp150", "CLNALLELEID", "CLNDN", "CLNDISDB", "CLNREVSTAT", "CLNSIG", "dbscSNV_ADA_SCORE", "dbscSNV_RF_SCORE", "GME_AF", "GME_NWA", "GME_NEA", "GME_AP", "GME_Israel", "GME_SD", "GME_TP", "GME_CA", "regsnp_fpr", "regsnp_disease", "regsnp_splicing_site", "cosmic89_coding", "cosmic89_noncoding", "info", "genotype", "allelic_depth","ref_depth", "alt_depth")

#Add bsID
a <- Auth(platform = "cavatica", token = "d8554e151de14d9b998af569994f4bfb")
p <- a$project(id = "gaonkark/sv-test")
file <- p$file(name = basename(input_txt), exact = TRUE)
testmaf$Tumor_Sample_Barcode <- file$meta()$`Kids First Biospecimen ID`

#Filtering
testmaf <- testmaf[as.numeric(testmaf$read_depth)>4 & testmaf$allele_fraction>0.15 & as.character(testmaf$FILTER)=="PASS",]
genes <- read.delim(input_genes, header = F, sep = "", stringsAsFactors = F)
testmaf <- testmaf[is.element(testmaf$Hugo_Symbol,genes$V1), ]
rm(genes)
testmaf <- testmaf[grep("Frame_Shift_Del|Frame_Shift_Ins|In_Frame_Del|In_Frame_Ins|Missense_Mutation|Nonsense_Mutation|Nonstop_Mutation|Splice_Region|Splice_Site|Targeted_Region|Translation_Start_Site",
                           testmaf$Variant_Classification, ignore.case = T),]

#Predictions
siftfilter <- testmaf[grep("D", testmaf$SIFT_pred, ignore.case = TRUE),]
lrtfilter <- testmaf[grep("D", testmaf$LRT_pred, ignore.case = TRUE),]
muttfilter <- testmaf[grep("A|D", testmaf$MutationTaster_pred, ignore.case = TRUE),]
mutafilter <- testmaf[grep("H|M", testmaf$MutationAssessor_pred, ignore.case = TRUE),]
pphdivfilter <- testmaf[grep("D|P", testmaf$Polyphen2_HDIV_pred, ignore.case = TRUE),]
pphvarfilter <- testmaf[grep("D|P", testmaf$Polyphen2_HVAR_pred, ignore.case = TRUE),]
metasvmfilter <- testmaf[grep("D", testmaf$MetaSVM_pred, ignore.case = TRUE),]
metalrfilter <- testmaf[grep("D", testmaf$MetaLR_pred, ignore.case = TRUE),]
mcapfilter <- testmaf[grep("D", testmaf$M.CAP_pred, ignore.case = TRUE),]
fathmmfilter <- testmaf[grep("D", testmaf$FATHMM_pred, ignore.case = TRUE),]
proveanfilter <- testmaf[grep("D", testmaf$PROVEAN_pred, ignore.case = TRUE),]
intervar <- testmaf[grep("pathogenic|uncertain",testmaf$InterVar_automated, ignore.case = T),]
clinvar <- testmaf[grep("pathogenic|uncertain|risk", testmaf$CLNSIG, ignore.case = T),]
merged <- rbind(clinvar, intervar, siftfilter,lrtfilter,muttfilter,pphdivfilter,pphvarfilter,metasvmfilter,metalrfilter,mcapfilter,fathmmfilter,proveanfilter,mutafilter)
rm("siftfilter", "lrtfilter", "muttfilter", "mutafilter", "pphdivfilter", "pphvarfilter", "metasvmfilter", "metalrfilter","mcapfilter", "fathmmfilter", "proveanfilter", "intervar", "clinvar")
merged <- unique(merged)

write.table(merged, output_txt, sep = "\t", quote = F, row.names = F)



