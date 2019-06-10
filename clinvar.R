args = commandArgs(trailingOnly = TRUE)
input_maf <- args[1]
output_maf <- sub(".maf",".filt.maf", input_maf)

#read input
input <- read.delim(input_maf, comment.char = "#", stringsAsFactors = F)
df <- as.data.frame(input)
#filter clin_sig
clinefilter <- mfilter[grep("benign", mfilter$CLIN_SIG, ignore.case = T, invert = T),]
#output
write.table(filtered, output_maf)