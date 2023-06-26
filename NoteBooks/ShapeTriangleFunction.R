ShapeTriangle <- function(path, df, ShapeSave=FALSE ){
library( ghyp )
library (timeSeries )
library ( fBasics )
source("xichiFunction.R")
  tryCatch({Vniki <- read.csv(paste(path,df, sep = ''))
  ## Return calculation
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
  
  rd <- c(1, 5, 10, 20, 40)
  yrets <- na.omit(matrix(unlist(lapply(rd, function(x) diff(log(y),
                                                             lag = x))), ncol = 5))
  
  ## HYP fitting
  hypfits <- apply(yrets, 2, fit.hypuv, symmetric = FALSE)
  points <- matrix(unlist(lapply(hypfits, xichi)), ncol = 2, byrow = TRUE)
  
  ## Shape Triangle
  col.def <- c("black", "blue", "red", "green", "orange")
  leg.def <- paste(rd, rep("day return", 5))
  if (ShapeSave == TRUE){  jpeg(gsub("csv","jpeg",df), )
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
    text(x = 0.0, y = -0.5, labels = df, srt = 0)
    
    
    dev.off()} else{  plot(points, ylim = c(-0.5, 1.2), xlim = c(-1.5, 1.5), col = col.def, pch = 16,
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
      text(x = 0.0, y = -0.5, labels = df, srt = 0)
      
    }
  

  }, error = function(e){
    message("An error occurred")
    print(e)},
  warning = function(w){
    message("A Warning occurred")
    print(w)
  }
  )
  print(df)
return(points)
}