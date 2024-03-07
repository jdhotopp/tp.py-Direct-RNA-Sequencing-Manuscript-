#!/usr/bin/env Rscript

# use bam2bed to get bed files
# arguments required: bedfile, sample name, region start, region end, strand, svg/jpeg
# no spaces in sample names

options(warn=-1)

library(ggplot2)
library(ggh4x)
library(gggenes)
library(patchwork)
library(svglite)

args = commandArgs(trailingOnly=TRUE)
if (length(args)<6) {
    stop("Please supply all required arguments:\n
    bedfile1, sample name, region start,
    region end, strand, svg/jpeg\n
    Please avoid using spaces in sample names.\n", call.=FALSE)
}


cat("\nReading in bed file...\n")
bed1 <- read.delim(args[1], header=FALSE)
colnames(bed1) <- c("Chrom", "Start", "End", "ReadID", "Score", "Strand")

region = c(as.numeric(args[3]), as.numeric(args[4]))

label1 = paste0(args[2], "_Start")
strand = args[5]

cat("Formatting the data...\n")
# sort by read starts
start_sort1 <- bed1[order(bed1$Start),]

    # separate positive strands
    start_sort_strand1 <- start_sort1[(start_sort1$Strand == strand),]
    #start_sort_strand2 <- start_sort2[(start_sort2$Strand == strand2),]

    # filter by region of interest
    start_region1 <- start_sort_strand1[(start_sort_strand1$Start >= region[1]) &
        (start_sort_strand1$End <= region[2]),]
    start_region1$pos <- c(nrow(start_region1):1)
    start_region1$label <- as.factor(label1)
    start_region1$plot_label <- as.factor(args[2])

# sort by read ends
# to use read end sorting, just change "start_region1" to "end_region1" in the plot below
end_sort1 <- bed1[order(bed1$End),]

    # separate positive strands
    end_sort_strand1 <- end_sort1[(end_sort1$Strand == strand),]

    # filter by region of interest
    end_region1 <- end_sort_strand1[(end_sort_strand1$Start >= region[1]) &
        (end_sort_strand1$End <= region[2]),]
    end_region1$pos <- c(1:nrow(end_region1))
    end_region1$label <- as.factor(label1)
    end_region1$plot_label <- as.factor(args[2])

cat("Plotting...\n")

p1 <- ggplot(start_region1) +
        geom_segment(aes(x=Start, y=pos, xend=End, yend=pos), size=1) +
        theme_classic() +
        scale_color_manual(values=c("black")) +
        xlim(as.numeric(args[3]), as.numeric(args[4])) +
        theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
            axis.title=element_blank(),
            panel.border=element_rect(fill=NA),
            legend.position = "none",
            strip.text = element_text(size = 12),
            axis.text.x=element_text(size=15)) #+

jpeg_filename = paste0(args[2],".jpeg")
svg_filename = paste0(args[2],".svg")
if(args[6] == "svg"){
    ggsave(svg_filename, p1, width = 10, height = 10)
}
if(args[6] == "jpeg"){
    ggsave(jpeg_filename, p1, width = 10, height = 10)
    }


