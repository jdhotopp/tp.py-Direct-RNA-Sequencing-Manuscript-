#!/usr/bin/env Rscript

library(GenomicAlignments)

args = commandArgs(trailingOnly=TRUE)

readbam = readGAlignments(args[1], use.names=TRUE)
bam_df <- as.data.frame(readbam)
read_ends = data.frame("end" = bam_df$end,"read_id" = rownames(bam_df))

write.table(read_ends, file="read_ends.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
~
