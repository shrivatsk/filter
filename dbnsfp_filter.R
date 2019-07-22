BiocManager::install("sevenbridges")
library(sevenbridges)
library(magrittr)
library(tidyr)
library(bedr)
args = commandArgs(trailingOnly = TRUE)
input_txt <- "testoutput_flags.txt"
input_genes <- "predisposition-genes.txt"
input_bed <- "pten_exon.bed"
threshold <- 0.001
output_txt <- sub(".snpEffdbnsfp_anno.txt",".snpEff.dbnsfp_anno.filt.all.txt", basename(input_txt))
print(output_txt)
#read input
input <- read.delim(input_txt, sep = "\t", comment.char = "#", stringsAsFactors = F)
print("input txt read")
df <- as.data.frame(input)
print("input txt converted to df")
rm(input)
#pass filter
passfilter <- df[grep("pass", df$FILTER, ignore.case = T),]
rm(df)
genes <- read.delim(input_genes, header = F, sep = "", stringsAsFactors = F)
#check against genes
genefilter <- passfilter[is.element(passfilter$ANN.0..GENE.,genes$V1), ]
rm(passfilter)
qfilter <- genefilter[grep("high|moderate",
                           genefilter$ANN.0..IMPACT., ignore.case = T),]
print("gene and pass filtered")
#check against threshold
qfilter[71:102] <- sapply(qfilter[71:102],as.numeric)
thresholdfilter <- qfilter[apply((qfilter[,71:102] <= threshold) | (is.na(qfilter[,71:102])), 1, all),]
rm(genefilter)
rm(qfilter)
rm(genes)
thresholdfilter <- thresholdfilter %>% separate_rows(ANN.0..EFFECT, sep = "&")
# siftfilter <- thresholdfilter[grep("D", thresholdfilter$dbNSFP_SIFT_pred, ignore.case = TRUE),]
# sift4Gfilter <- thresholdfilter[grep("D", thresholdfilter$dbNSFP_SIFT4G_pred, ignore.case = TRUE),]
# lrtfilter <- thresholdfilter[grep("D", thresholdfilter$dbNSFP_LRT_pred, ignore.case = TRUE),]
# muttfilter <- thresholdfilter[grep("A|D", thresholdfilter$dbNSFP_MutationTaster_pred, ignore.case = TRUE),]
# pphdivfilter <- thresholdfilter[grep("D|P", thresholdfilter$dbNSFP_Polyphen2_HDIV_pred, ignore.case = TRUE),]
# pphvarfilter <- thresholdfilter[grep("D|P", thresholdfilter$dbNSFP_Polyphen2_HVAR_pred, ignore.case = TRUE),]
# metasvmfilter <- thresholdfilter[grep("D", thresholdfilter$dbNSFP_MetaSVM_pred, ignore.case = TRUE),]
# metalrfilter <- thresholdfilter[grep("D", thresholdfilter$dbNSFP_MetaLR_pred, ignore.case = TRUE),]
# mcapfilter <- thresholdfilter[grep("D", thresholdfilter$dbNSFP_M_CAP_pred, ignore.case = TRUE),]
# primateaifilter <- thresholdfilter[grep("D", thresholdfilter$dbNSFP_PrimateAI_pred, ignore.case = TRUE),]
# fathmmfilter <- thresholdfilter[grep("D", thresholdfilter$dbNSFP_FATHMM_pred, ignore.case = TRUE),]
# proveanfilter <- thresholdfilter[grep("D", thresholdfilter$dbNSFP_PROVEAN_pred, ignore.case = TRUE),]
# mutafilter <- thresholdfilter[grep("H|M", thresholdfilter$dbNSFP_MutationAssessor_pred, ignore.case = TRUE),]
# aloftfilter <- thresholdfilter[grep("D", thresholdfilter$dbNSFP_Aloft_pred, ignore.case = TRUE),]
# deogenfilter <- thresholdfilter[grep("D", thresholdfilter$dbNSFP_DEOGEN2_pred, ignore.case = TRUE),]
# clinvarfilter <- thresholdfilter[grep("pathogen", thresholdfilter$dbNSFP_clinvar_clnsig, ignore.case = T),]
# merged <- rbind(clinvarfilter, siftfilter, sift4Gfilter,lrtfilter,muttfilter,pphdivfilter,pphvarfilter,metasvmfilter,metalrfilter,mcapfilter,primateaifilter,fathmmfilter,proveanfilter,mutafilter,aloftfilter,deogenfilter)
# merged <- unique(merged)
# nonmissensefilter <- thresholdfilter[grep("missense", thresholdfilter$ANN.0..EFFECT, ignore.case = T, invert = T),]
# merged <- rbind(merged, nonmissensefilter)
# merged <- unique(merged)
names(thresholdfilter)[names(thresholdfilter)== "ANN.0..GENE."] <- "Hugo_Symbol"
names(thresholdfilter)[names(thresholdfilter)== "CHROM"] <- "Chromosome"
names(thresholdfilter)[names(thresholdfilter)== "POS"] <- "Start_Position"
for (i in 1:nrow(thresholdfilter)) {
  thresholdfilter$End_Position[i] <- thresholdfilter$Start_Position[i] + max(unlist(lapply(as.list(unlist(strsplit(thresholdfilter$ALT[i], ","))), nchar))) - nchar(thresholdfilter$REF[i])
  startpos <- thresholdfilter$Start_Position[i]
  endpos <- thresholdfilter$End_Position[i]
  if(endpos-startpos > 0) {
    thresholdfilter$Variant_Type[i] <- "INS"
  } else if (endpos-startpos < 0) {
    thresholdfilter$Variant_Type[i] <- "DEL"
  } else {
    thresholdfilter$Variant_Type[i] <- "SNP"
  }
}
a <- Auth(platform = "cavatica", token = "d8554e151de14d9b998af569994f4bfb")
p <- a$project(id = "gaonkark/sv-test")

file <- p$file(name = basename(input_txt), exact = TRUE)
thresholdfilter$Tumor_Sample_Barcode <- file$meta()$`Kids First Biospecimen ID`
thresholdfilter <- cbind( thresholdfilter[,(1:3)], thresholdfilter[,(103:105)], thresholdfilter[,4:102] )

names(thresholdfilter)[names(thresholdfilter)== "ANN.0..EFFECT"] <- "Variant_Classification"
names(thresholdfilter)[names(thresholdfilter)== "REF"] <- "Reference_Allele"
names(thresholdfilter)[names(thresholdfilter)== "ALT"] <- "Tumor_Seq_Allele2"



# merged$Variant_Classification[merged$Variant_Classification == "missense_variant"] <- "Missense_Mutation"
# merged$Variant_Classification[merged$Variant_Classification == "protein_protein_contact"] <- "Missense_Mutation"
# merged$Variant_Classification[merged$Variant_Classification == "stop_gained"] <- "Nonsense_Mutation"
# merged$Variant_Classification[merged$Variant_Classification == "start_lost"] <- "Translational_Start_Site"
# merged$Variant_Classification[merged$Variant_Classification == "splice_donor_variant"] <- "Splice_Site"
# merged$Variant_Classification[merged$Variant_Classification == "splice_acceptor_variant"] <- "Splice_Site"
# merged$Variant_Classification[(merged$Variant_Classification == "conservative_inframe_deletion") | (merged$Variant_Classification == "disruptive_inframe_deletion")] <- "In_Frame_Del"
# merged$Variant_Classification[(merged$Variant_Classification == "conservative_inframe_insertion") | (merged$Variant_Classification == "disruptive_inframe_insertion")] <- "In_Frame_Ins"
# merged$Variant_Classification[(merged$Variant_Classification == "frameshift_variant") & (merged$Variant_Type == "INS")] <- "Frame_Shift_Ins"
# merged$Variant_Classification[(merged$Variant_Classification == "frameshift_variant") & (merged$Variant_Type == "DEL")] <- "Frame_Shift_Del"
write.table(thresholdfilter, output_txt, sep = "\t", quote = F, row.names = F)

# "chromosome_number_variation","exon_loss_variant","frameshift_variant","rare_amino_acid_variant","splice_acceptor_variant","splice_donor_variant","start_lost","stop_gained","stop_lost","transcript_ablation","coding_sequence_variant","conservative_inframe_deletion","conservative_inframe_insertion","disruptive_inframe_deletion","disruptive_inframe_insertion","missense_variant","regulatory_region_ablation","splice_region_variant","TFBS_ablation"




