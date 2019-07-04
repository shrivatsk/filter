install.packages("sevenbridges")
install.packages("tidyverse")
library(tidyverse)
library(sevenbridges)
args = commandArgs(trailingOnly = TRUE)
input_txt <- args[1]
input_genes <- args[2]
threshold <- args[3]
output_txt <- sub("snpeff.dbnsfp_anno","snpEff.dbnsfp_anno.filt", basename(input_txt))
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
#keep SIFT predictions
siftfilter <- thresholdfilter[grep("D", thresholdfilter$dbNSFP_SIFT_pred, ignore.case = TRUE),]
sift4Gfilter <- thresholdfilter[grep("D", thresholdfilter$dbNSFP_SIFT4G_pred, ignore.case = TRUE),]
lrtfilter <- thresholdfilter[grep("D", thresholdfilter$dbNSFP_LRT_pred, ignore.case = TRUE),]
muttfilter <- thresholdfilter[grep("A|D", thresholdfilter$dbNSFP_MutationTaster_pred, ignore.case = TRUE),]
pphdivfilter <- thresholdfilter[grep("D|P", thresholdfilter$dbNSFP_Polyphen2_HDIV_pred, ignore.case = TRUE),]
pphvarfilter <- thresholdfilter[grep("D|P", thresholdfilter$dbNSFP_Polyphen2_HVAR_pred, ignore.case = TRUE),]
metasvmfilter <- thresholdfilter[grep("D", thresholdfilter$dbNSFP_MetaSVM_pred, ignore.case = TRUE),]
metalrfilter <- thresholdfilter[grep("D", thresholdfilter$dbNSFP_MetaLR_pred, ignore.case = TRUE),]
mcapfilter <- thresholdfilter[grep("D", thresholdfilter$dbNSFP_M_CAP_pred, ignore.case = TRUE),]
primateaifilter <- thresholdfilter[grep("D", thresholdfilter$dbNSFP_PrimateAI_pred, ignore.case = TRUE),]
fathmmfilter <- thresholdfilter[grep("D", thresholdfilter$dbNSFP_FATHMM_pred, ignore.case = TRUE),]
proveanfilter <- thresholdfilter[grep("D", thresholdfilter$dbNSFP_PROVEAN_pred, ignore.case = TRUE),]
mutafilter <- thresholdfilter[grep("H|M", thresholdfilter$dbNSFP_MutationAssessor_pred, ignore.case = TRUE),]
aloftfilter <- thresholdfilter[grep("D", thresholdfilter$dbNSFP_Aloft_pred, ignore.case = TRUE),]
deogenfilter <- thresholdfilter[grep("D", thresholdfilter$dbNSFP_DEOGEN2_pred, ignore.case = TRUE),]
clinvarfilter <- thresholdfilter[grep("pathogen", thresholdfilter$dbNSFP_clinvar_clnsig, ignore.case = T),]
merged <- rbind(clinvarfilter, siftfilter, sift4Gfilter,lrtfilter,muttfilter,pphdivfilter,pphvarfilter,metasvmfilter,metalrfilter,mcapfilter,primateaifilter,fathmmfilter,proveanfilter,mutafilter,aloftfilter,deogenfilter)
merged <- unique(merged)

names(merged)[names(merged)== "ANN.0..GENE."] <- "Hugo_Symbol"
names(merged)[names(merged)== "CHROM"] <- "Chromosome"
names(merged)[names(merged)== "POS"] <- "Start_Position"
for (i in 1:nrow(merged)) {
  merged$End_Position[i] <- merged$Start_Position[i] + max(unlist(lapply(as.list(unlist(strsplit(merged$ALT[i], ","))), nchar))) - nchar(merged$REF[i])
  startpos <- merged$Start_Position[i]
  endpos <- merged$End_Position[i]
  if(endpos-startpos > 0) {
    merged$Variant_Type[i] <- "INS"
  } else if (endpos-startpos < 0) {
    merged$Variant_Type[i] <- "DEL"
  } else {
    merged$Variant_Type[i] <- "SNP"
  }
}
a <- Auth(platform = "cavatica", token = "d8554e151de14d9b998af569994f4bfb")
p <- a$project(id = "gaonkark/sv-test")

file <- p$file(name = basename(input_txt), exact = TRUE)
merged$Tumor_Sample_Barcode <- file$meta()$`Kids First Biospecimen ID`
merged <- cbind( merged[,(1:3)], merged[,(103:105)], merged[,4:102] )

names(merged)[names(merged)== "ANN.0..EFFECT"] <- "Variant_Classification"
names(merged)[names(merged)== "REF"] <- "Reference_Allele"
names(merged)[names(merged)== "ALT"] <- "Tumor_Seq_Allele2"


merged <- merged %>% separate_rows(Variant_Classification, sep = "&")
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
write.table(merged, output_txt, sep = "\t", quote = F, row.names = F)

# "chromosome_number_variation","exon_loss_variant","frameshift_variant","rare_amino_acid_variant","splice_acceptor_variant","splice_donor_variant","start_lost","stop_gained","stop_lost","transcript_ablation","coding_sequence_variant","conservative_inframe_deletion","conservative_inframe_insertion","disruptive_inframe_deletion","disruptive_inframe_insertion","missense_variant","regulatory_region_ablation","splice_region_variant","TFBS_ablation"




