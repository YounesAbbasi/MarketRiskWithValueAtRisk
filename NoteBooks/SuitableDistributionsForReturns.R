library(NormalLaplace)
library(numDeriv)
library(MASS)
library(timeDate)
library( ghyp )
library (timeSeries )
library ( fBasics )
library(e1071)
#########----------------Fitting HPW returns to the GHD-------------#########

 Vniki <- read.csv('D:/R_project/MarketRiskModeling/Raw/Vniki.csv')
## Return calculation
startpoint=2100
endpoint=4806
y <- timeSeries(Vniki[startpoint:endpoint,"Adj.Close"],
                charvec = as.character(Vniki[startpoint:endpoint,2]))
yret <- na.omit(diff(log(y))*100)

## Fitting
ef <- density(yret)
ghdfit <- fit.ghypuv(yret, symmetric = FALSE, control = list(maxit = 1000))
hypfit <- fit.hypuv(yret, symmetric = FALSE, control = list(maxit = 1000))
nigfit <- fit.NIGuv(yret, symmetric = FALSE, control = list(maxit = 1000))

## Densities
ghddens <- dghyp(ef$x, ghdfit)
hypdens <- dghyp(ef$x, hypfit)
nigdens <- dghyp(ef$x, nigfit)
nordens <- dnorm(ef$x, mean = mean(yret), sd=sd(c(yret[,1])))
col.def <- c('black', 'red', 'blue', 'green', 'orange')
plot(ef, xlab = '', ylab = expression(f(x)), ylim = c(0, 0.45))
lines(ef$x, ghddens, col = 'red')
lines(ef$x, hypdens, col = 'blue')
lines(ef$x, nigdens, col = 'green')
lines(ef$x, nordens, col = 'orange')
legend('topleft', legend = c("empirical" , "GHD" , "HYP" , "NIG" , "NORM"),
       col = col.def,  lty = 1)

## QQ-Plots
qqghyp(ghdfit, line = TRUE, ghyp.col = 'red',
       plot.legend = FALSE, gaussian = FALSE,
       main = '', cex = 0.8)
qqghyp(hypfit, add = TRUE, ghyp.pch = 2, ghyp.col = 'blue',
       gaussian = FALSE, line = FALSE, cex = 0.8)
qqghyp(nigfit, add = TRUE, ghyp.pch = 3, ghyp.col = 'green',
       gaussian = FALSE, line = FALSE, cex = 0.8)
legend('topleft', legend = c("GHD" , "HYP" , "NIG"),
       col = col.def[-c(1,5)], pch = 1:3)

## Diagnostics
AIC <- stepAIC.ghyp(yret, dist = c('ghyp', 'hyp', 'NIG'), symmetric = FALSE, 
                    control = list(maxit=1000))
LRghdnig <- lik.ratio.test(ghdfit , nigfit)
LRghdyp <- lik.ratio.test(ghdfit , hypfit)

#########----------VaR and ES derived from the GHD, HYP, and NIG------#########


## Probabilities
p <- seq(0.001 , 0.05 , 0.001)

## VaR
ghd.VaR <- abs(qghyp(p, ghdfit))
hyp.VaR <- abs(qghyp(p, hypfit))
nig.VaR <- abs(qghyp(p, nigfit))
nor.VaR <- abs(qnorm(p, mean = mean(yret), 
                     sd = sd(c(yret[,1]))))
emp.VaR <- abs(quantile(x = yret, probs = p))

## Plot of VaR
plot(emp.VaR, type  = 'l', xlab = '', ylab = 'VaR',
     axes = FALSE, ylim=range(c(hyp.VaR, nig.VaR, ghd.VaR, 
                                nor.VaR, emp.VaR)))
box()
axis(1, at = seq(along=p), labels = names(emp.VaR),
     tick = FALSE)
axis(2, at = pretty(range(emp.VaR, ghd.VaR, hyp.VaR,
                          nig.VaR, nor.VaR)))
lines(seq(along = p), ghd.VaR, col = 'red')
lines(seq(along = p), hyp.VaR, col = 'blue')
lines(seq(along = p), nig.VaR, col = 'green')
lines(seq(along = p), nor.VaR, col = 'orange')
legend('topright',
       legend = c('Empirical', 'GHD', 'HYP', 'NIG', 'Normal'),
       col = col.def, lty = 1)

## ES
ghd.ES <- abs(ESghyp(p, ghdfit))
hyp.ES <- abs(ESghyp(p, hypfit))
nig.ES <- abs(ESghyp(p, nigfit))
nor.ES <- abs(mean(yret)-sd(c(yret[,1]))*
                dnorm(qnorm(1-p))/p)
obs.p <- ceiling(p * length(yret))
emp.ES <- sapply(obs.p, function(x) abs(mean(sort(c(yret))
                                             [1:x])))

## Plot of ES
plot(emp.ES, type = 'l', xlab = '', ylab = 'ES', axes = FALSE,
     ylim = range(c(hyp.ES, nig.ES, ghd.ES, nor.ES, emp.ES)))
box()
axis(1, at = 1:length(p), labels = names(emp.VaR),
     tick = FALSE)
axis(2, at = pretty(range(emp.ES, ghd.ES, hyp.ES, nig.ES,
                          nor.ES)))
lines(1:length(p), ghd.ES, col = 'red')
lines(1:length(p), hyp.ES, col = 'blue')
lines(1:length(p), nig.ES, col = 'green')
lines(1:length(p), nor.ES, col = 'orange')
legend('topright',
       legend = c('Empirical', 'GHD', 'HYP', 'NIG', 'Normal'),
       col = col.def, lty = 1)

#########--------------Shape triangle for HYP distribution------------#########
rd <- c(1, 5, 10, 20, 40)
yrets <- na.omit(matrix(unlist(lapply(rd, function(x) diff(log(y),
                                      lag = x))), ncol = 5))

## Function for xi/chi coefficients
xichi <-  function(x){
  param <-  coef(x, type = 'alpha.delta')
  rho <- param[['beta']]/param[['alpha']]
  zeta <- param[['delta']] * sqrt(param[['alpha']]^2 -
                                    param[['beta']]^2)
  xi <- 1/sqrt(1+zeta)
  chi <- xi * rho
  result <- c(chi, xi)
  names(result) <- c("chi","xi")
  return(result)
}

## HYP fitting
hypfits <- apply(yrets, 2, fit.hypuv, symmetric = FALSE)
points <- matrix(unlist(lapply(hypfits, xichi)), ncol = 2, byrow = TRUE)

## Shape Triangle
col.def <- c("black", "blue", "red", "green", "orange")
leg.def <- paste(rd, rep("day return", 5))
plot(points, ylim = c(-0.5, 1.2), xlim = c(-1.5, 1.5), col = col.def, pch = 16,
     ylab = expression(xi), xlab = expression(chi))

lines(x = c(0, -1), y = c(0, 1))
lines(x = c(0, 1), y = c(0, 1))
lines(x = c(-1, 1), y = c(1,1))

legend("bottomright", legend = leg.def, col = col.def, pch = 16)
text(x = 0.0, y = 1.05, label = "Laplace", srt = 0)
text(x = -1.0, y = 1.05, labels = "Exponential", srt = 0)
text(x = 1.0, y = 1.05, labels = "Exponential", srt = 0)
text(x = 0.0, y = -0.1, labels = "Normal", srt = 0)
text(x = -0.9, y = 0.5, labels = "Hyperbolic, left skewed", srt = 302)
text(x = 0.9, y = 0.5, labels = "Hyperbolic, right skewed", srt = 57)


#########----------------Fitting QCOM returns to the GLD-------------#########

## Loading of packages 
library(lmomco )
library(FRAPO )

## Data loading

Idx <- yrets[,2]
L <- -100 * Idx

## Computing VaR ( Normal & GLD ) 99%, moving window
windowlen = 104
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
plot(Res[,"Loss"], type = "p", xlab = "Time Index",
     ylim = c(-15, max(Res)))
abline(h = 0, col = "grey")
lines(Res[, "VaRGld"], col = "blue", lwd = 2)
lines(Res[, "VaRNor"], col = "red", lwd = 2)

legend("bottomright", legend = c("Losses", "VaR GLD", "VaR Normal"),
       col = c("black", "blue", "red"),
       lty = c(NA, 1, 1), pch = c(19, NA, NA), bty = "n")

########----- FTSE 100 stocks: shape triangle of standardized GLD-----########

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
