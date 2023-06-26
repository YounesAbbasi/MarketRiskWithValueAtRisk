
## Function for xi/chi coefficients
xichi <-  function(x){
  library(NormalLaplace)
  library(numDeriv)
  library(MASS)
  library(timeDate)
  library( ghyp )
  library (timeSeries )
  library ( fBasics )
  library(e1071)
  library(data.table)
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