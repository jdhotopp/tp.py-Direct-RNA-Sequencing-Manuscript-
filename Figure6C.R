# ---
# Title: "Script to generate Figure 6C"
# ---

# Dependancies

library("ggplot2")

# Function to determine mode

# Input files 

Depth23S<-"merged_all.23S.sorted.depth"
Depth16S<-"merged_all.16S.sorted.depth"
DepthSINIVT<-"SINV_IVT_20220512_allpos.txt"

# Loading data into data frame

Depth_23S <- as.data.frame(read.delim(Depth23S, sep="\t", header=FALSE))
Depth_16S <- as.data.frame(read.delim(Depth16S, sep="\t", header=FALSE))
SINVDepth <- as.data.frame(read.delim(DepthSINIVT,sep="\t", header=FALSE))

# Labeling the data frames

Depth_23S$V1<-"23S"
Depth_16S$V1<-"16S"
SINVDepth$V1<-"SINV"

# Normalizing the depth

Depth_23S$V3<-Depth_23S$V3/max(Depth_23S$V3)
Depth_16S$V3<-Depth_16S$V3/max(Depth_16S$V3)
SINVDepth$V3<-SINVDepth$V3/max(SINVDepth$V3)

# Setting the index for the x-axis to "1"

Depth_23S$V2<-rev(c(1:length(Depth_23S$V2)))
Depth_16S$V2<-rev(c(1:length(Depth_16S$V2)))
SINVDepth$V2<-rev(c(1:length(SINVDepth$V2)))

# File Processing for Ribosomal depth comparison

FullDepth<-rbind(Depth_23S,Depth_16S,SINVDepth)
colnames(FullDepth)<-c("Transcript","RelativePosition","NormalizedDepth")

# Ribosomal Depth comparison plot as PDF

pdf("Transcript.Depth.Comparisons.allpos.pdf",useDingbats=FALSE)
ggplot(FullDepth[FullDepth$Transcript == "23S" | FullDepth$Transcript == "SINV" | FullDepth$Transcript == "16S",], aes(x=RelativePosition, y=NormalizedDepth, color=Transcript)) + geom_smooth() + theme_bw(base_size=24)+xlim(0,3000)
dev.off()


