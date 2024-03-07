#!/usr/bin/env Rscript

# use bam2bed to get bed files
# arguments required: bedfile1, bedfile2, sample1 name, sample2 name, region start, region end, strand
# no spaces in sample names

options(warn=-1)

library(ggplot2)
library(ggh4x)
library(gggenes)
library(patchwork)
library(svglite)

args = commandArgs(trailingOnly=TRUE)
if (length(args)<11) {
    stop("Please supply all required arguments:\n
    bedfile1, bedfile2, sample1 name, sample2 name, region start,
    region end, strand1, strand2, predictions bedfile, gff file, seqid, svg/jpeg\n
    Please avoid using spaces in sample names. Seqid argument must be identical to seqid in gff file.\n", call.=FALSE)
}


cat("\nReading in bed file 1...\n")
bed1 <- read.delim(args[1], header=FALSE)
colnames(bed1) <- c("Chrom", "Start", "End", "ReadID", "Score", "Strand")
bed1 <- bed1[(bed1$Chrom == args[11]),]

cat("Reading in bed file 2...\n")
bed2 <- read.delim(args[2], header=FALSE)
colnames(bed2) <- c("Chrom", "Start", "End", "ReadID", "Score", "Strand")
bed2 <- bed2[(bed2$Chrom == args[11]),]

cat("Reading in prediction file...\n")
pred <- read.delim(args[9], header=FALSE)
colnames(pred) <- c("Chrom", "Start", "End", "ReadID", "Score", "Strand")
pred <- pred[(pred$Chrom == args[11]),]

cat("Reading in gff file...\n")
gff_unfilt <- read.delim(args[10], header=FALSE, comment.char="#")
colnames(gff_unfilt) <- c("Chrom", "Source", "Type", "Start", "End", "Score",
    "Strand", "Phase", "Attributes")
gff <- gff_unfilt[(gff_unfilt$Chrom == args[11]),]

region = c(as.numeric(args[5]), as.numeric(args[6]))

label1 = paste0(args[3], "_Start")
label2 = paste0(args[4], "_Start")
label3 = paste0(args[3], "_End")
label4 = paste0(args[4], "_End")
label5 = "Prediction_fwd"
label6 = "Prediction_rev"
data_labels = c(label1, label2, label3, label4, label5, label6)

strand1 = args[7]
strand2 = args[8]

cat("Formatting the data...\n")
# sort by read starts
start_sort1 <- bed1[order(bed1$Start),]
start_sort2 <- bed2[order(bed2$Start),]

    # separate positive strands
    # just change sign if needing antisense strand
    start_sort_strand1 <- start_sort1[(start_sort1$Strand == strand1),]
    start_sort_strand2 <- start_sort2[(start_sort2$Strand == strand2),]

    # filter by region of interest
    start_region1 <- start_sort_strand1[(start_sort_strand1$Start >= region[1]) &
        (start_sort_strand1$End <= region[2]),]
    start_region1$pos <- c(1:nrow(start_region1))
    start_region1$label <- as.factor(data_labels[1])
    start_region1$plot_label <- as.factor(args[3])
    start_region2 <- start_sort_strand2[(start_sort_strand2$Start >= region[1]) &
        (start_sort_strand2$End <= region[2]),]
    start_region2$pos <- c(1:nrow(start_region2))
    start_region2$label <- as.factor(data_labels[2])
    start_region2$plot_label <- as.factor(args[4])


# sort by read ends
end_sort1 <- bed1[order(bed1$End),]
end_sort2 <- bed2[order(bed2$End),]

    # separate positive strands
    end_sort_strand1 <- end_sort1[(end_sort1$Strand == strand1),]
    end_sort_strand2 <- end_sort2[(end_sort2$Strand == strand2),]

    # filter by region of interest
    end_region1 <- end_sort_strand1[(end_sort_strand1$Start >= region[1]) &
        (end_sort_strand1$End <= region[2]),]
    end_region1$pos <- c(1:nrow(end_region1))
    end_region1$label <- as.factor(data_labels[3])
    end_region1$plot_label <- as.factor(args[3])
    end_region2 <- end_sort_strand2[(end_sort_strand2$Start >= region[1]) &
        (end_sort_strand2$End <= region[2]),]
    end_region2$pos <- c(1:nrow(end_region2))
    end_region2$label <- as.factor(data_labels[4])
    end_region2$plot_label <- as.factor(args[4])

# format predictions
   
    # filter by region of interest
    pred_region <- pred[(pred$Start >= region[1]) &
        (pred$End <= region[2]),]
    
    # switch rev coordinates for plotting
    pred_fwd <- pred_region[(pred_region$Strand == "+"),]
    pred_rev <- pred_region[(pred_region$Strand == "-"),]
    colnames(pred_rev) <- c("Chrom", "End", "Start", "ReadID", "Score", "Strand")
    pred_plot_df <- rbind(pred_fwd, pred_rev)

# format gff

    # use only rows with "gene" under the "Type" column
    gff_gene <- gff[(gff$Type == "gene"),]

    # filter by region of interest
    gff_region <- gff_gene[(gff_gene$Start >= region[1]) &
        (gff_gene$End <= region[2]),]

    # get "gene" attribute from last column
    getAttributeField <- function (x, field, attrsep = ";") {
     s = strsplit(x, split = attrsep, fixed = TRUE)
     sapply(s, function(atts) {
         a = strsplit(atts, split = "=", fixed = TRUE)
         m = match(field, sapply(a, "[", 1))
         if (!is.na(m)) {
             rv = a[[m]][2]
         }
         else {
             rv = as.character(NA)
         }
         return(rv)
     })
}

    gff_region$Gene <- getAttributeField(as.character(gff_region$Attributes), "gene")

    # create a coordinate in the middle of start and end for the gene label position
    gff_region$Mid <- gff_region$Start + ((gff_region$End - gff_region$Start)/2)

    # switch rev coordinates for plotting
    gff_fwd <- gff_region[(gff_region$Strand == "+"),]
    gff_rev <- gff_region[(gff_region$Strand == "-"),]
    colnames(gff_rev) <- c("Chrom", "Source", "Type", "End", "Start", "Score",
    "Strand", "Phase", "Attributes", "Gene", "Mid")
    gff_plot_df <- rbind(gff_fwd, gff_rev)

cat("Plotting...\n")
    
# create plot for predictions

pred_arrow_size <- ifelse(nrow(pred_plot_df) > 20, 1.5, 3)

p1 <- ggplot(pred_plot_df, aes(xmin = Start, xmax = End, y = ReadID, fill=Strand, forward = TRUE)) +
         geom_gene_arrow(arrowhead_height = unit(pred_arrow_size, "mm"), 
            arrow_body_height = unit(pred_arrow_size, "mm"),
            arrowhead_width = unit(1, "mm")) +
            theme_classic() +
            xlim(as.numeric(args[5]), as.numeric(args[6])) +
            theme(axis.text=element_blank(), 
            axis.ticks=element_blank(),
            axis.title=element_blank(),
            panel.border=element_rect(fill=NA),
            legend.position = "none") +
            scale_fill_manual(values=c("olivedrab4", "olivedrab2"))


# create plot for gff input

nudge_amt <- ifelse(nrow(gff_plot_df) < 10, 1, 1.6) 
arrow_size <- ifelse(nrow(gff_plot_df) > 20, 2, 3)

p2 <- ggplot(gff_plot_df, aes(xmin = Start, xmax = End, y = Gene, fill=Strand, forward = TRUE)) +
         geom_gene_arrow(arrowhead_height = unit(arrow_size, "mm"), 
            arrowhead_width = unit(1, "mm"),
            arrow_body_height = unit(arrow_size, "mm")) +
            geom_text(aes(x = Mid, label = Gene), 
               nudge_y = nudge_amt, check_overlap = T) +
            theme_classic() +
            xlim(as.numeric(args[5]), as.numeric(args[6])) +
            theme(axis.text.y=element_blank(), 
                axis.ticks.y=element_blank(),
                axis.title=element_blank(),
                axis.text.x=element_text(size = 15),
                panel.border=element_rect(fill=NA),
                legend.position = "none") +
            scale_fill_manual(values=c("grey45", "grey70")) +
            scale_y_discrete(expand = c(0.3, 0.3))
 

# create plot for alignments
plot_df <- rbind(start_region1, start_region2, end_region1, end_region2)

size_diff = round(nrow(start_region1) / nrow(start_region2), digits=0)
panel_diff = ifelse(size_diff <= 4, size_diff, 4)
panel_size = c(panel_diff, 1, panel_diff, 1)

p3 <- ggplot(plot_df) +
        geom_segment(aes(x=Start, y=pos, xend=End, yend=pos, color=label), size=1) +
        theme_classic() + 
        scale_color_manual(values=c("tomato", "red3", "steelblue2", "blue3")) +
        xlim(as.numeric(args[5]), as.numeric(args[6])) +
        theme(axis.text=element_blank(), axis.ticks=element_blank(),
            axis.title=element_blank(),
            panel.border=element_rect(fill=NA),
            legend.position = "none",
            strip.text = element_text(size = 12)) +
        facet_wrap(~ label, ncol = 1, 
            scales = "free_y", 
            strip.position = "right") +
        force_panelsizes(rows = panel_size)

pred_height <- ifelse(nrow(pred_plot_df) > 30, 2, 1)
p_svg <- p3 / p1 / p2 + plot_layout(heights = c(5, pred_height, 1))

jpeg_filename = paste0(args[3],"_",args[4],".jpeg")
svg_filename = paste0(args[3],"_",args[4],".svg")
pdf_filename = paste0(args[3],"_",args[4],".pdf")
if(args[12] == "svg"){
    ggsave(svg_filename, p_svg, width = 14, height = 18)
}
if(args[12] == "jpeg"){
    ggsave(jpeg_filename, p_svg, width = 14, height = 18)
    }

if(args[12] == "pdf"){
    ggsave(pdf_filename, p_svg, width = 14, height = 18)
    }
        
