library(mixOmics)
library(ggpubr)
library(gridExtra)
library(ggplot2)
library(PLSDAbatch)
library(vegan)

dfp.clr <- read.csv("dfp_batcheffects_100%_remove_686&695.csv",header = TRUE,row.names = NULL)


## Convert imported intensity file to correct format
dfp.clr <- as.matrix(dfp.clr)
dfp.clr <- t(dfp.clr)
colnames(dfp.clr) <- NULL
dfp.clr <- dfp.clr[-1,]


## Add in batch info as factor
dfp.batch.treat <- read.csv("dfp_batcheffects_batch&treat_remove686_695.csv",header = TRUE)

dfp.batch <- factor(dfp.batch.treat$Batch)
dfp.trt <- factor(dfp.batch.treat$Treatment)

## Create scatterplot of PCA
dfp.pca.before <- pca(dfp.clr, scale = TRUE)

#dfp.pca.before <- pca(dfp.corrected, scale = TRUE)

Scatter_Density(object = dfp.pca.before, batch = dfp.batch, trt = dfp.trt)


## Create boxplot for batch comparison
dfp.df <- data.frame(value = dfp.clr[,1], batch = dfp.batch)

box_plot(df = dfp.df, title = 'A1A5B6', x.angle = 30)

## Create a density plot for each batch
density_plot(df = dfp.df, title = 'A1A5B6')


dfp.clr.scale <- scale(dfp.corrected, center = T, scale = T)

dfp.factors.df <- data.frame(trt = dfp.trt, batch = dfp.batch)

rda.res <- varpart(dfp.clr.scale, ~ trt, ~ batch, data = dfp.factors.df)

dfp.prop.df <- data.frame(Treatment = NA, Intersection = NA, Batch = NA, Residuals = NA)
dfp.prop.df[1,] <- rda.res$part$indfract$Adj.R.squared

dfp.prop.df[dfp.prop.df < 0] = 0
dfp.prop.df <- as.data.frame(t(apply(dfp.prop.df, 1, function(x){x/sum(x)})))


partVar_plot(prop.df = dfp.prop.df)

#### Managing batch effects

# Linear Regression

install.packages("performance")

library(performance)

install.packages("lmerTest")

library(lmerTest)

dfp.lm <- linear_regres(data = dfp.clr, trt = dfp.trt, batch = dfp.batch, type = 'linear model')


n <- ncol(dfp.clr)

dfp.corrected <- matrix(, nrow = 21, ncol = n)

for(column in 1:n){
  dfp.corrected[,column] <- dfp.lm[["model"]][[column]][["fitted.values"]]
}

rownames <- row.names(dfp.clr)
row.names(dfp.corrected) <- rownames

#### PLSDA for batch effects.

dfp_plsda_batch <- PLSDA_batch(dfp.clr, dfp.trt, dfp.batch, ncomp.trt = 4, ncomp.bat = 1)

dfp_X.corrected <- dfp_plsda_batch$X.nobatch

write.csv(dfp_X.corrected,file="dfp_batcheffects_PLSDA_100%_remove_686&695.csv")