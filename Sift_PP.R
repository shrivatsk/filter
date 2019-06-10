args = commandArgs(trailingOnly = TRUE)
input_maf <- args[1]
input_genes <- args[2]
threshold <- args[3]
output_maf <- sub(".maf",".filt.maf", input_maf)
#read input
input <- read.delim(input_maf, comment.char = "#", stringsAsFactors = F)
df <- as.data.frame(input)
#filter only keep pass
passfilter <- df[grep("pass", df$FILTER, ignore.case = T),]
#check against threshold
thresholdfilter <- passfilter[apply((passfilter[,124:132]<= 0.001) | (is.na(passfilter[124:132])), 1, all),]
genes <- read.delim(input_genes, header = F, sep = "", stringsAsFactors = F)
#check against genes
genefilter <- thresholdfilter[is.element(thresholdfilter$Hugo_Symbol,genes$V1), ]
#keep nonsynonymous 
mfilter <- genefilter[grep("Frame_Shift_Del|Frame_Shift_Ins|In_Frame_Del|In_Frame_Ins|Missense_Mutation|Nonsense_Mutation|Nonstop_Mutation|Splice_Region|Splice_Site|Targeted_Region|Translation_Start_Site",
                         genefilter$Variant_Classification, ignore.case = T),]
#keep SIFT predictions
siftfilter <- mfilter[grep( "tolerated", mfilter$SIFT, ignore.case = TRUE, invert = TRUE),]
#keep PolyPhen predictions
ppfilter <- mfilter[grep("benign", mfilter$PolyPhen, ignore.case = TRUE, invert = TRUE), ]
#Union of SIFT and Polyphen
dupl <- rownames(siftfilter) %in% rownames(ppfilter)
filtered <- rbind(siftfilter, ppfilter[!dupl,])
write.table(filtered, output_maf)
