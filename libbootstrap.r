binit <- function(values, factors, minValue, maxValue, step) {
  breaks <- seq(minValue, maxValue - step, step)
  histcat <- as.factor(sapply(factors, function(x) { sum(x > breaks) }))
  
  newdf <- data.frame(binID = seq(1, length(breaks)), 
                      lower = breaks, 
                      midpoints = breaks + (step / 2))
  
  Responsemean <- tapply(values, histcat, mean)
  Responselength <- tapply(values, histcat, length)
  Response.df <- data.frame(binID = rownames(Responsemean), 
                            ResponseMean = Responsemean, 
                            N = Responselength, 
                            freq = Responselength / sum(Responselength))
  
  finaldf <- merge(newdf, Response.df, by = "binID", all.x = TRUE)
  finaldf
}

get_corrcoef_per_gene <- function(data, gene.colname, xname, yname) {
    tapply(data$index, data[, gene.colname], function(y) {
               cor.test(data[y, xname], data[y, yname], method = "spearman")[4][[1]]
  })
}

#get_rho_per_gene

get_lm_per_gene <- function(data, gene.colname, xname, yname) {
  tapply(data$index, data[, gene.colname], function(y) {
    lm(data[y, yname] ~ data[y, xname])
  })
}

get_pval_per_gene <- function(mod) {
  sapply(mod, function(x) {
    anova(x)[[5]][[1]]
  })
}

get_slope_per_gene <- function(mod) {
  sapply(mod, function(x) {
    summary(x)$coefficients[[2]]
  })
}

load.male.quality.measure.1 <- function(fname)
{
    data <- read.table(fname, header = TRUE)
    data$Rep <- as.factor(data$Rep)
    data$RALline <- as.factor(data$RALline)
    data$LabTempWinston <- as.numeric(data$LabTempWinston)
    data$totOffspring <- as.numeric(data$totOffspring)
    data$oneOverLatCop <- 1/(data$LatCop)

    #################################################################
    #For each measure, take the median per RALline - transform lat cop to 1/latCop
    #################################################################
    OneOverLatCop <- aggregate(data$oneOverLatCop, by=list(data$RALline), median)
    colnames(OneOverLatCop) <- c("RALline", "MedOneOverLatCop")
    
    #################################################################
    #Revision 2 - scale each measure - mean = 0, SD = 1
    #################################################################
    #get means of the two phenotypes
    MeanLatCop <- mean(OneOverLatCop$MedOneOverLatCop)
    #scale to make means of zero
    OneOverLatCop$MedLatCop_minus_MeanLatCop <- OneOverLatCop$MedOneOverLatCop - MeanLatCop
    #get sd of two phenotypes scaled for means of zero
    Stdev_newLatCop <- sd(OneOverLatCop$MedLatCop_minus_MeanLatCop)
    #scale to make sd of 1 (of the means scaled phenotypes)
    OneOverLatCop$Scale_LatCop <- OneOverLatCop$MedLatCop_minus_MeanLatCop / Stdev_newLatCop
    
    OneOverLatCop
}


load.male.quality.measure.2 <- function(fname)
{
    data <- read.table(fname, header = TRUE)
    #get means of the two phenotypes
    MeanProm <- mean(data$MedianPromiscuity)
    #scale to make means of zero
    data$MedProm_minus_MeanProm <- data$MedianPromiscuity - MeanProm
    #get sd of two phenotypes scaled for means of zero
    Stdev_newProm <- sd(data$MedProm_minus_MeanProm)
    #scale to make sd of 1 (of the means scaled phenotypes)
    data$Scale_Prom <- data$MedProm_minus_MeanProm / Stdev_newProm

    data
}

load.huang.data <- function(fname)
{
    b <- read.table(fname)
    ColumnsIwant <- c("gene", "X208.M", "X301.M", "X303.M", "X304.M", "X307.M", "X313.M", "X315.M", "X324.M", "X335.M", "X357.M", "X358.M", "X360.M", "X362.M", "X365.M", "X375.M", "X379.M", "X380.M", "X391.M", "X399.M", "X427.M", "X437.M", "X486.M", "X517.M", "X555.M", "X639.M", "X705.M", "X707.M", "X712.M", "X714.M", "X730.M", "X765.M", "X774.M", "X799.M", "X820.M", "X852.M","av_females", "av_males", "Log2FC", "sexBias_pval", "sexBias_pvalfdr", "SexBias", "SigSexBias")
    
    #loop through each gene and calculate correlation between GE and male trait 
    #just do 4 genes at first in the dataset Genes and then do for all
    b2 <- b[, ColumnsIwant]
    dim(subset(b2, SigSexBias == 'M'))
    dim(subset(b2, SigSexBias == 'F'))
    melted <- melt(b2, id.vars=c("gene","av_females", "av_males", "Log2FC", "sexBias_pval", "sexBias_pvalfdr", "SexBias", "SigSexBias"))
    melted$Line <- substr(melted$variable, 2,4)

    melted
}


my.seq <- function(x, y) {
  seq(x, x + (y - 1))
}


iterate.slope <- function(data, expression, boot.slope, boot.pvals,
                          niters, minValue, maxValue, step)
{
    min.index <- tapply(data$index, data$RALline, min)
    ngenes <- length(unique(data$gene))
    
    for (i in 1:niters) {
        rlines <- sample(min.index, length(min.index), replace = TRUE)
        rindices <- lapply(rlines, my.seq, ngenes)
        sub.genes <- data[unlist(rindices), ]
        sub.genes$index <- seq(1:length(sub.genes$RALline))
        
        lms.Scale_Prom <- get_lm_per_gene(sub.genes, "gene", "value", "Scale_Prom")
        Scale_Prom.slope <- get_slope_per_gene(lms.Scale_Prom)
        
        Scale_Prom.df <- data.frame(gene = names(Scale_Prom.slope), 
                                    Slope = Scale_Prom.slope) 
        rownames(Scale_Prom.df) <- 1:length(Scale_Prom.df$gene)
        
        SlopeDataProm <- merge(Scale_Prom.df, expression, by.x = 'gene', by.y = 'gene')
        
        
        alldfProm <- binit(SlopeDataProm$Slope, SlopeDataProm$Log2FC, minValue, maxValue, step)
        lmfitProm <- lm(ResponseMean ~ midpoints , data=alldfProm, weights = freq, na.action = na.exclude)
        boot.pvals[i] <- anova(lmfitProm)[[5]][[1]]
        boot.slope[i] <- lmfitProm$coefficients[[2]][1]
        print(paste(i, boot.slope[i], boot.pvals[i], sep = " "))
    }
}


iterate.rho <- function(data, expression, boot.slope, boot.pvals,
                        niters, minValue, maxValue, step)
{
    min.index <- tapply(data$index, data$RALline, min)
    ngenes <- length(unique(data$gene))
    
    for (i in 1:niters) {
        rlines <- sample(min.index, length(min.index), replace = TRUE)
        rindices <- lapply(rlines, my.seq, ngenes)
        sub.genes <- data[unlist(rindices), ]
        sub.genes$index <- seq(1:length(sub.genes$RALline))
        
        Scale_Prom.rho <- get_corrcoef_per_gene(sub.genes, "gene", "value", "Scale_Prom")
        
        Scale_Prom.df <- data.frame(gene = names(Scale_Prom.rho), 
                                    rho = Scale_Prom.rho) 
        rownames(Scale_Prom.df) <- 1:length(Scale_Prom.df$gene)
        
        SlopeDataProm <- merge(Scale_Prom.df, expression, by.x = 'gene', by.y = 'gene')
        
        
        alldfProm <- binit(SlopeDataProm$rho, SlopeDataProm$Log2FC, minValue, maxValue, step)
        lmfitProm <- lm(ResponseMean ~ midpoints , data=alldfProm, weights = freq, na.action = na.exclude)
        boot.pvals[i] <- anova(lmfitProm)[[5]][[1]]
        boot.slope[i] <- lmfitProm$coefficients[[2]][1]
        print(paste(i, boot.slope[i], boot.pvals[i], sep = " "))
    }
    boot.df <- data.frame(pvals = boot.pvals, slope = boot.slope)
    boot.df
}
