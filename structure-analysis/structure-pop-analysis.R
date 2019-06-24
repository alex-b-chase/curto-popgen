setwd('/Volumes/Bennett_BACKUP/Research/curtobacterium/chapter4_popgen/mauve_new/')

install.packages(c("fields","RColorBrewer","mapplots"), repos = "http://cran.us.r-project.org")
source("http://bioconductor.org/biocLite.R")
biocLite("LEA")

source("http://membres-timc.imag.fr/Olivier.Francois/Conversion.R")
source("http://membres-timc.imag.fr/Olivier.Francois/POPSutilities.R")

library(LEA)

# struct2geno(file = 'core_alignment1500.snp', TESS = FALSE, diploid = FALSE, FORMAT = 1,
#             extra.row = 0, extra.col = 0, output = "genotype.geno")

obj.at = snmf("./genotype.geno", K = 1:20, ploidy = 1, entropy = T,
              CPU = 6, project = "new")

pdf('entropy.pdf', height = 10, width = 8)

plot(obj.at, col = "blue4", cex = 1.4, pch = 19)

dev.off()

# polz lab clustering identified 7 populations, so use that since it is vague from entropy output...
obj.snmf = snmf("./genotype.geno", K = 7, alpha = 100, project = "new")
qmatrix = Q(obj.snmf, K = 7)

pdf('structure-analysis_K7.pdf', height = 6, width = 12)

barplot(t(qmatrix), col = c("orange","violet","yellow", "black", "red", "lightgreen", "blue"), border = NA, space = 0,
        xlab = "Individuals", ylab = "Admixture coefficients")

dev.off()

obj.snmf = snmf("./genotype.geno", K = 4, alpha = 1000, project = "new", ploidy = 1)
qmatrix = Q(obj.snmf, K = 4)

write.table(qmatrix, file = "K4-lea-admixture.txt", sep = '\t', quote = F)

pdf('structure-analysis_K4.pdf', height = 6, width = 12)

barplot(t(qmatrix), col = c("dodgerblue","coral","yellow", "green"), border = NA, space = 0,
        xlab = "Individuals", ylab = "Admixture coefficients")

 dev.off()




dat.str  <- matrix(sample(c(101:105,-9),
                          200, prob = c(rep(1,5), 0.1),
                          replace = TRUE),
                   nrow = 10, ncol = 20)
write.table(dat.str,
            file = "dat.str",
            col.names = FALSE,
            row.names = FALSE,
            quote = FALSE)

### Conversion
struct2geno("dat.str", diploid = 2, FORMAT = 1)
