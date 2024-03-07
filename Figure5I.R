
## Read in Libraries
library(tidyverse)
devtools::source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")


## Read in File
df <- readr::read_delim("021524_E2348_69_6salmonIlluminashort_2salmonaONTlong_6FADUIlluminashort.txt",
                        skip = 1, delim = "\t", col_names = F)

df <- df %>% unite("X1", X1:X3, sep = "_", remove = TRUE)
df <- df %>% column_to_rownames("X1")
colnames(df) <- str_to_upper(letters[1:14])

## remove columns C, K
## Swap H and G
df <- df %>% select(-'C', -'K')
df <- df %>% select("A", "B", "D", "E", "F", "H", "G", "I", "J", "L", "M", "N")


heatmap_matrix <- as.matrix(log2(df+1))

groups <- data.frame((colnames(df)))
colnames(groups) <- "Sample"


## Plot Heatmap
col_fun = circlize::colorRamp2(c(0, 7, 18), c("white", "firebrick1", "firebrick3"))
breaks <- seq(min(heatmap_matrix),max(heatmap_matrix),by=.5)
hmcol.log <- colorRampPalette(c("white","firebrick1","firebrick3"))(length(breaks)-1)
colcol <- as.matrix(as.character(groups$Colour))


pdf("Figure5I.pdf",
    width = 10,
    height = 10)
heatmap.3(heatmap_matrix,
          col=hmcol.log,
          trace="none",
          labRow = NA,
          labCol = colnames(heatmap_matrix),
          Rowv = T,
          Colv = F,
          ColSideColors=colcol,
          lhei = c(2,8,1),
          breaks = breaks,
          density.info = 'histogram',
          dendrogram = "row")

invisible(dev.off())