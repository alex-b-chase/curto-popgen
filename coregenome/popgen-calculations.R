setwd('/Volumes/Bennett_BACKUP/Research/curtobacterium/chapter4_popgen/roary_coregenome/coregenes90/')

######## ######## ######## ######## ######## ######## ######## 
####### calculate dNdS ratios
######## ######## ######## ######## ######## ######## ######## 

# run screen genes before to eliminate genes with "bad" reading frames (i.e. not multiple of 3bp)

library(seqinr)
library(Biostrings)
library(ggplot2)

genes <- list.files("goodalignments/")
X <- read.alignment(file = paste("goodalignments/", genes[1], sep=""), format = "fasta")
Y <- matrix(nrow = 0, ncol = (X$nb - 1))
dN <- rep(0, length(genes))
dS <- rep(0, length(genes))
vdN <- rep(0, length(genes))
vdS <- rep(0, length(genes))
z <- rep(0, length(genes))
for(i in 1:length(genes)){
  X <- read.alignment(file = paste("goodalignments/", genes[i], sep=""), format = "fasta")
  #Y <- rbind(Y, kaks(X)$ka[1:(X$nb - 1)])
  dN[i] <- kaks(X)$ka[1]
  dS[i] <- kaks(X)$ks[1]
  vdN[i] <- kaks(X)$vka[1]
  vdS[i] <- kaks(X)$vks[1]
  z[i] <- (dN[i] - dS[i])/(sqrt(vdN[i] + vdS[i]))
}
p <- 2*(1-pnorm(abs(z)))
padj <- p.adjust(p, method = "BH")
lp <- -log10(padj)
plot(1:length(p), -log10(padj))

# remove the infinities
dS[dS < 0] <- 0
dNdS <- dN/dS
dNdS[!is.finite(dNdS)] <- 0

# Make table
genelist <- sapply(strsplit((list.files("goodalignments/")), split="\\."), function(x) x[[1]])

Dtable <- data.frame(genelist,  dN, dS, dNdS, z, p, padj)
colnames(Dtable) <- c("gene", "dN", "dS", "dNdS", "Z", "P", "Adjusted.P")

desc <- read.csv(file = "../coregeneIDs90.txt", sep = '\t', header = T)

z <- merge(Dtable, desc, by = "gene")

geneloc <- read.table('../coregenes90/genelocation.txt', header = T, sep = '\t')
final <- merge(z, geneloc, by = 'gene')

library(cowplot)

pdf("dnds-frequencies.pdf", height = 10, width = 10)

q1 <- ggplot(data=final, aes(log2(final$dNdS))) + 
  geom_histogram(breaks=seq(-8, 1, by = 0.5), 
                 col="black", 
                 fill="lightgray") + 
  labs(title="Directional Selection") +
  labs(x="dN/dS", y="Count") +
  theme_bw()

q2 <- ggplot(final, aes(y = -log10(Adjusted.P), x = ref.location)) +
  geom_point(aes(colour = log2(dNdS) >= 0)) +
  scale_colour_manual(name = 'dNdS > 0', values = setNames(c('red','lightgray'),c(T, F))) +
  theme_bw() +
  ylim(c(0, 15)) +
  ylab(label = "-log(P)") +
  xlab(label = "Reference Genome Position")

plot_grid(q1, q2, labels = "AUTO", ncol = 1, align = 'v')

dev.off()

write.csv(final, file = "dnds.allgenes.csv", qmethod='escape', quote=TRUE, row.names = F)


