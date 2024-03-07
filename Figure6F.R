library("ggplot2")
library(scales)

data <- as.data.frame(read.delim("Table_for_Figure6F_header.txt", sep="\t", header=TRUE))

pdf("Figure6F.pdf",useDingbats=FALSE)

ggplot(data, aes(x=Length, y=log2Ratio)) + stat_density_2d(aes(fill = ..level..), geom = "polygon") + scale_fill_gradient(low=muted("blue"), high="yellow") + geom_smooth(method = "lm", formula = y ~ x) + theme_bw(base_size=24) + xlab("Transcript Length (kbp)") + ylab("log2(ONT TPM / Illumina TPM)") + xlim(0,5) + ylim(-7,7)

ggplot(data, aes(x=Length, y=log2Ratio)) + geom_point(size=.1) + geom_smooth(method = "lm", formula = y ~ x) + theme_bw(base_size=20) + xlab("Transcript Length (kbp)") + ylab("log2(ONT TPM / Illumina TPM)") + xlim(0,15) + ylim(-12,12) 

dev.off()

lm(data$log2Ratio ~ data$Length)