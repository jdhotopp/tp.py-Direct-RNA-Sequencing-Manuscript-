# ---
# title: "Script to plot a compariosn of TPMs for Salmon/ONT, Salmon/Illumina, FADU/Illumina
# ---

# Dependancies

library("ggplot2")
library(scales)

# Input file

de <- as.data.frame(read.delim("021524_E2348_69_6salmonIlluminashort_2salmonaONTlong_6FADUIlluminashort_with_header.txt", sep="\t", header=TRUE))


de$logSI1 <- log2(de$SalIll1)
de$logFI1 <- log2(de$FADUIll1)

de$logSI2 <- log2(de$SalIll2)
de$logFI2 <- log2(de$FADUIll2)

de$logSI3 <- log2(de$SalIll3)
de$logFI3 <- log2(de$FADUIll3)

de$logSI4 <- log2(de$SalIll4)
de$logFI4 <- log2(de$FADUIll4)

de$logSI5 <- log2(de$SalIll5)
de$logFI5 <- log2(de$FADUIll5)

de$logSI6 <- log2(de$SalIll6)
de$logFI6 <- log2(de$FADUIll6)

de$logSO1 <- log2(de$SalONT1)
de$logSO2 <- log2(de$SalONT2)


# Plot

pdf("SalmonIllumina1_FADUIllumina1.pdf",useDingbats=FALSE)

ggplot(de, aes(x=logSI1, y=logFI1))  + stat_density_2d(aes(fill = ..level..), geom = "polygon") + scale_fill_gradient(low=muted("blue"), high="yellow") + theme_bw(base_size=24) + xlab("Salmon Illumina ERR3932853 TPM") + ylab("FADU Illumina ERR3932853 TPM") 

ggplot(de, aes(x=logSI1, y=logFI1))  + geom_point(size=.1) + theme_bw(base_size=20) + xlab("log2(Salmon Illumina ERR3932853 TPM)") + ylab("log2(FADU Illumina ERR3932853 TPM)")

dev.off()

pdf("SalmonIllumina2_FADUIllumina2.pdf",useDingbats=FALSE)

ggplot(de, aes(x=logSI2, y=logFI2))  + stat_density_2d(aes(fill = ..level..), geom = "polygon") + scale_fill_gradient(low=muted("blue"), high="yellow") + theme_bw(base_size=24) + xlab("Salmon Illumina ERR3932854 TPM 2") + ylab("FADU Illumina ERR3932854 TPM 2") 

ggplot(de, aes(x=logSI2, y=logFI2))  + geom_point(size=.1) + theme_bw(base_size=20) + xlab("log2(Salmon Illumina ERR3932854 TPM 2)") + ylab("log2(FADU Illumina ERR3932854 TPM 2)")

dev.off()

pdf("SalmonIllumina3_FADUIllumina3.pdf",useDingbats=FALSE)

ggplot(de, aes(x=logSI3, y=logFI3))  + stat_density_2d(aes(fill = ..level..), geom = "polygon") + scale_fill_gradient(low=muted("blue"), high="yellow") + theme_bw(base_size=24) + xlab("Salmon Illumina ERR3932855 TPM 3") + ylab("FADU Illumina ERR3932855 TPM 3") 

ggplot(de, aes(x=logSI3, y=logFI3))  + geom_point(size=.1) + theme_bw(base_size=20) + xlab("log2(Salmon Illumina ERR3932855 TPM 3)") + ylab("log2(FADU Illumina ERR3932855 TPM 3)")

dev.off()

pdf("SalmonIllumina4_FADUIllumina4.pdf",useDingbats=FALSE)

ggplot(de, aes(x=logSI4, y=logFI4))  + stat_density_2d(aes(fill = ..level..), geom = "polygon") + scale_fill_gradient(low=muted("blue"), high="yellow") + theme_bw(base_size=24) + xlab("Salmon Illumina ERR3932847 TPM 4") + ylab("FADU Illumina ERR3932847 TPM 4") 

ggplot(de, aes(x=logSI4, y=logFI4))  + geom_point(size=.1) + theme_bw(base_size=20) + xlab("log2(Salmon Illumina ERR3932847 TPM 4)") + ylab("log2(FADU Illumina ERR3932847 TPM 4)")

dev.off()

pdf("SalmonIllumina5_FADUIllumina5.pdf",useDingbats=FALSE)

ggplot(de, aes(x=logSI5, y=logFI5))  + stat_density_2d(aes(fill = ..level..), geom = "polygon") + scale_fill_gradient(low=muted("blue"), high="yellow") + theme_bw(base_size=24) + xlab("Salmon Illumina ERR3932848 TPM 5") + ylab("FADU Illumina ERR3932848 TPM 5") 

ggplot(de, aes(x=logSI5, y=logFI5))  + geom_point(size=.1) + theme_bw(base_size=20) + xlab("log2(Salmon Illumina ERR3932848 TPM 5)") + ylab("log2(FADU Illumina ERR3932848 TPM 5)")

dev.off()

pdf("SalmonIllumina6_FADUIllumina6.pdf",useDingbats=FALSE)

ggplot(de, aes(x=logSI6, y=logFI6))  + stat_density_2d(aes(fill = ..level..), geom = "polygon") + scale_fill_gradient(low=muted("blue"), high="yellow") + theme_bw(base_size=24) + xlab("Salmon Illumina TPM 6") + ylab("FADU Illumina TPM 6") 

ggplot(de, aes(x=logSI6, y=logFI6))  + geom_point(size=.1) + theme_bw(base_size=20) + xlab("log2(Salmon Illumina ERR3932849 TPM 6)") + ylab("log2(FADU Illumina ERR3932849 TPM 6)")

dev.off()

pdf("SalmonIllumina1_SalmonONT1.pdf",useDingbats=FALSE)

ggplot(de, aes(x=logSI1, y=logSO1))  + stat_density_2d(aes(fill = ..level..), geom = "polygon") + scale_fill_gradient(low=muted("blue"), high="yellow") + theme_bw(base_size=24) + xlab("Salmon Illumina ERR3932853 TPM 1") + ylab("Salmon ONT SRR18061003 TPM 1") 

ggplot(de, aes(x=logSI1, y=logSO1))  + geom_point(size=.1) + theme_bw(base_size=20) + xlab("log2(Salmon Illumina ERR3932853 TPM 1)") + ylab("log2(Salmon ONT SRR18061003 TPM 1)")

dev.off()

pdf("SalmonIllumina2_SalmonONT2.pdf",useDingbats=FALSE)

ggplot(de, aes(x=logSI2, y=logSO2))  + stat_density_2d(aes(fill = ..level..), geom = "polygon") + scale_fill_gradient(low=muted("blue"), high="yellow") + theme_bw(base_size=24) + xlab("Salmon Illumina ERR3932854 TPM 2") + ylab("Salmon ONT SRR18061004 2") 

ggplot(de, aes(x=logSI2, y=logSO2))  + geom_point(size=.1) + theme_bw(base_size=20) + xlab("log2(Salmon Illumina ERR3932854 TPM 2)") + ylab("log2(Salmon ONT SRR18061004 TPM 2)")

dev.off()
