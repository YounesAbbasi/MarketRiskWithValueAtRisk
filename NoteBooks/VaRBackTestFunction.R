VaRBackTest <- function(path, df, ShapeSave = FALSE, windowlen = 104, maxDataPoints = 1000){
options(warn=-1)
library(lmomco)
library(FRAPO)
Vniki <- read.csv(paste(path,df, sep = ''))
## Return calculation
if (dim(Vniki)[1]>maxDataPoints){
  endpoint=dim(Vniki)[1]
  startpoint=endpoint-maxDataPoints
} else{
  startpoint=1
  endpoint=dim(Vniki)[1]
}
y <- timeSeries(Vniki[startpoint:endpoint,"Adj.Close"],
                charvec = as.character(Vniki[startpoint:endpoint,2]))
yret <- na.omit(diff(log(y))*100)

rd <- c(1, 5, 10, 20, 40)
yrets <- na.omit(matrix(unlist(lapply(rd, function(x) diff(log(y),
                                                           lag = x))), ncol = 5))
Idx <- yrets[,2]
L <- -100 * Idx

## Computing VaR ( Normal & GLD ) 99%, moving window
ep <- windowlen:length(L)
sp <- 1:length(ep)
level <- 0.99

VaR <- matrix(NA, ncol = 2, nrow = length(ep))
for(i in 1:length(sp)){
  x <- L[sp[i]:ep[i]]
  lmom <- lmom.ub(x)
  fit <- pargld(lmom)
  VaRGld <- quagld(level, fit)
  VaRNor <- qnorm(level, mean(x), sd(x))
  VaR[i,] <- c(VaRGld, VaRNor)
}
## Summarising results
Res <-  cbind(L[windowlen+1:length(L)], VaR[-nrow(VaR), ])
colnames(Res) <- c("Loss", "VaRGld", "VaRNor")
  
  ## Plot of backtest results
if (ShapeSave == TRUE){  jpeg(paste(Vniki[1,5],"_VaRBackTest.jpeg", sep = ''), )
plot(Res[,"Loss"], type = "p", xlab = paste("Time Index_", df, sep = ''),
     ylim = c(-15, max(Res)))
abline(h = 0, col = "grey")
lines(Res[, "VaRGld"], col = "blue", lwd = 2)
lines(Res[, "VaRNor"], col = "red", lwd = 2)

legend("topleft", legend = c("Losses", "VaR GLD", "VaR Normal"),
       col = c("black", "blue", "red"),
       lty = c(NA, 1, 1), pch = c(19, NA, NA), bty = "n")

dev.off()
} else{
plot(Res[,"Loss"], type = "p", xlab = paste("Time Index_", df, sep = ''),
     ylim = c(-15, max(Res)))
abline(h = 0, col = "grey")
lines(Res[, "VaRGld"], col = "blue", lwd = 2)
lines(Res[, "VaRNor"], col = "red", lwd = 2)

legend("topleft", legend = c("Losses", "VaR GLD", "VaR Normal"),
       col = c("black", "blue", "red"),
       lty = c(NA, 1, 1), pch = c(19, NA, NA), bty = "n")

}
  }
