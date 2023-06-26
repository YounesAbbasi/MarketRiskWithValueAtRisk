#########----------------Fitting Stock returns to the GHD-------------#########
VaRES <- function(Vniki, DensityImageSave = FALSE, QQImageSave = FALSE, VaRImageSave = FALSE, ESImageSave = FALSE){
  library( ghyp )
  library (timeSeries )
  library ( fBasics )
  if (dim(Vniki)[1]>3000){
    endpoint=dim(Vniki)[1]
    startpoint=endpoint-3000
  } else{
    startpoint=1
    endpoint=dim(Vniki)[1]
  }
  
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
  if (DensityImageSave==TRUE){
    jpeg(paste(Vniki[1,5],"_density.jpeg", sep = ''), )
    plot(ef, xlab = '', ylab = expression(f(x)), ylim = c(0, 0.45))
    lines(ef$x, ghddens, col = 'red')
    lines(ef$x, hypdens, col = 'blue')
    lines(ef$x, nigdens, col = 'green')
    lines(ef$x, nordens, col = 'orange')
    legend('topleft', legend = c("empirical" , "GHD" , "HYP" , "NIG" , "NORM"),
           col = col.def,  lty = 1)
    name = paste("density_",Vniki[1,5], sep = '')
    dev.off()
  } else{plot(ef, xlab = '', ylab = expression(f(x)), ylim = c(0, 0.45))
    lines(ef$x, ghddens, col = 'red')
    lines(ef$x, hypdens, col = 'blue')
    lines(ef$x, nigdens, col = 'green')
    lines(ef$x, nordens, col = 'orange')
    legend('topleft', legend = c("empirical" , "GHD" , "HYP" , "NIG" , "NORM"),
           col = col.def,  lty = 1)
    name = paste("density_",Vniki[1,5], sep = '')}


  ## QQ-Plots
  if (QQImageSave==TRUE){
  jpeg(paste(Vniki[1,5],"_QQ.jpeg", sep = ''), )
  qqghyp(ghdfit, line = TRUE, ghyp.col = 'red',
         plot.legend = FALSE, gaussian = FALSE,
         main = '', cex = 0.8)
  qqghyp(hypfit, add = TRUE, ghyp.pch = 2, ghyp.col = 'blue',
         gaussian = FALSE, line = FALSE, cex = 0.8)
  qqghyp(nigfit, add = TRUE, ghyp.pch = 3, ghyp.col = 'green',
         gaussian = FALSE, line = FALSE, cex = 0.8)
  legend('topleft', legend = c("GHD" , "HYP" , "NIG"),
         col = col.def[-c(1,5)], pch = 1:3)
  dev.off()} else{
    qqghyp(ghdfit, line = TRUE, ghyp.col = 'red',
           plot.legend = FALSE, gaussian = FALSE,
           main = '', cex = 0.8)
    qqghyp(hypfit, add = TRUE, ghyp.pch = 2, ghyp.col = 'blue',
           gaussian = FALSE, line = FALSE, cex = 0.8)
    qqghyp(nigfit, add = TRUE, ghyp.pch = 3, ghyp.col = 'green',
           gaussian = FALSE, line = FALSE, cex = 0.8)
    legend('topleft', legend = c("GHD" , "HYP" , "NIG"),
           col = col.def[-c(1,5)], pch = 1:3)
  }
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
  if (VaRImageSave==TRUE){
  jpeg(paste(Vniki[1,5],"_VaR.jpeg", sep = ''), )
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
  dev.off()} else{
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
  }
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
  if (ESImageSave==TRUE){
  jpeg(paste(Vniki[1,5],"_ES.jpeg", sep = ''), )
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
  dev.off()} else{
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
  }
}

