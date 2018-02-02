library(reshape2) 
source("libbootstrap.r")

OneOverLatCop <- load.male.quality.measure.1("1.MaleQualityMeasures_20160526.txt")
Data <- load.male.quality.measure.2("2.MaleQualityMeasures_Averages_20160331.txt")
melted <- load.huang.data("DGRP_Huang_expression.txt")

AllGenes <- merge(Data, melted, by.x='RALline', by.y='Line')
AllGenes$index <- 1:length(AllGenes$gene)

LatCopTransformed <- merge(OneOverLatCop, melted, by.x='RALline', by.y='Line')
LatCopTransformed$index <- 1:length(LatCopTransformed$gene)

GeneExpression <- melted[,c("gene", "av_males", "SigSexBias", "Log2FC")]
GeneExpression <- subset(GeneExpression, !duplicated(gene))

######

minValue <- -10
maxValue <- 10
step <- 0.1

######
# Real data

rho.Scale_Prom <- get_corrcoef_per_gene(AllGenes, "gene", "value", "Scale_Prom")
Scale_Prom.df <- data.frame(gene = names(rho.Scale_Prom), rho = rho.Scale_Prom) 
rownames(Scale_Prom.df) <- 1:length(Scale_Prom.df$gene)
SlopeDataProm <- merge(Scale_Prom.df, GeneExpression, by.x = 'gene', by.y = 'gene')
alldfProm <- binit(SlopeDataProm$rho, SlopeDataProm$Log2FC, minValue, maxValue, step)
lmfitProm <- lm(ResponseMean ~ midpoints , data=alldfProm, weights = freq, na.action = na.exclude)
real.pval <- anova(lmfitProm)[[5]][[1]]
real.slope <- lmfitProm$coefficients[[2]][1]

#####
# Bootstrap

niters <- 1000
boot.slope.prom <- numeric(niters)
boot.pvals.prom <- numeric(niters)

run.start <- date()
#boot.df <- iterate.rho(AllGenes, GeneExpression, boot.slope.prom, boot.pvals.prom, niters, minValue, maxValue, step)
print(run.start)
print(date())

#write.table(boot.df, "boot.df")
