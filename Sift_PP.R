args = commandArgs(trailingOnly = TRUE)
input_maf <- args[1]
input_genes <- args[2]
threshold <- args[3]
output_maf <- sub(".maf",".filt.maf", basename(input_maf))
print(output_maf)
#read input
input <- read.delim(input_maf, comment.char = "#", stringsAsFactors = F)
print("input maf read")
df <- as.data.frame(input)
print("input maf converted to df")
#filter only keep pass
passfilter <- df[grep("pass", df$FILTER, ignore.case = T),]
print("pass filter done")
#check against threshold
thresholdfilter <- passfilter[apply((passfilter[,124:132]<= threshold) | (is.na(passfilter[124:132])), 1, all),]
print("threshold filter done")
genes <- read.delim(input_genes, header = F, sep = "", stringsAsFactors = F)
print("genes read from txt")
#check against genes
genefilter <- thresholdfilter[is.element(thresholdfilter$Hugo_Symbol,genes$V1), ]
print("gene filter done")
#keep nonsynonymous 
mfilter <- genefilter[grep("Frame_Shift_Del|Frame_Shift_Ins|In_Frame_Del|In_Frame_Ins|Missense_Mutation|Nonsense_Mutation|Nonstop_Mutation|Splice_Region|Splice_Site|Targeted_Region|Translation_Start_Site",
                         genefilter$Variant_Classification, ignore.case = T),]
print("nonsynonymous filter done")
#keep SIFT predictions
siftfilter <- mfilter[grep("tolerated", mfilter$SIFT, ignore.case = TRUE, invert = TRUE),]
print("sift filter done")
#keep PolyPhen predictions
ppfilter <- mfilter[grep("benign", mfilter$PolyPhen, ignore.case = TRUE, invert = TRUE), ]
print("polyphen filter done")
#Union of SIFT and Polyphen
dupl <- rownames(siftfilter) %in% rownames(ppfilter)
filtered <- rbind(siftfilter, ppfilter[!dupl,])
print("sift and polyphen combined")
write.table(filtered, output_maf, sep = "\t", quote = F, row.names = F)
print("written to output file")
