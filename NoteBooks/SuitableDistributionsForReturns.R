source("VaRESFunction.R")
source('ShapeTriangleFunction.R')
source('VaRBackTestFunction.R')
#########----------------Import Data-------------#########
path = 'D:/R_project/MarketRiskModeling/Raw/AdjPrice/'
dfList = list.files(path, pattern = '*.csv')

#########-----------Fitting Stock returns to the GHD --------#############
Vniki <- read.csv(paste(path,dfList[1], sep = ''))
VaRES(Vniki, DensityImageSave = TRUE, QQImageSave = TRUE, VaRImageSave = TRUE, ESImageSave = TRUE )

#########--------------Shape triangle for HYP distribution------------#########
for (df in dfList[1:1]){
  ShapeTriangle(path=path, df = df, ShapeSave=TRUE)
}
ShapeTriangle(path=path, df = dfList[1])

#########----------------Fitting QCOM returns to the GLD-------------#########
VaRBackTest(path = path, df = dfList[1], ShapeSave = TRUE, maxDataPoints=300)

########----- FTSE 100 stocks: shape triangle of standardized GLD-----########
#### This part is for some stocks not one stock
library(FRAPO)
library(fBasics)

## Loading data
data(INDTRACK3)
P <- INDTRACK3[,-1]
R <- returnseries(P, method = "discret", trim = TRUE)

## Fitting and calculating beta and lambda
Fit <- apply(R, 2, gldFit, method ="rob", doplot = FALSE,
             trace = FALSE)
DelaBetaParam <- matrix(unlist(lapply(Fit, function(x){
  l <- x@fit$estimate[c(3,4)]
  res <- c(l[2] - l[1], l[1]+l[2])
  res
})), ncol = 2, byrow = TRUE)

## Shape triangle
plot(DelaBetaParam, xlim = c(-2,2), ylim = c(-2,0),
     xlab = expression(delta == lambda[4]-lambda[3]),
     ylab = expression(beta == lambda[3]+lambda[4]),
     pch = 19, cex = 0.5)
segments(x0 = -2, y0 = -2, x1 = 0, y1 = 0,
         col = 'grey', lwd = 0.8, lty = 2)
segments(x0 = 2, y0 = -2, x1 = 0, y1 = 0,
         col = "grey", lwd = 0.8, lty = 2)
segments(x0 = 0, y0 = -2, x1 = 0, y1 = 0, col = "blue",
         lwd = 0.8, lty = 2)
segments(x0 = -0.5, y0 = -0.5, x1 = 0.5, y1 = -0.5,
         col = 'red', lwd = 0.8, lty = 2)
segments(x0 = -1, y0 = -1, x1 = 1, y1 = -1,
         col = 'red', lwd = 0.8, lty = 2)

