#Funciones para analisis de residuales de objetos de clase nls.
#pruebas Durbin-Watson y Breusch-Pagan
#https://github.com/hugosal
#18/03/2020

#Estas funciones son adaptaciones de las funciones durbinWatsonTest (car version 3.0-6)
#Fox, J. (2016) Applied Regression Analysis and Generalized Linear Models, Third Edition. Sage.
#Y bptest (lmtest version 0.9-37) T.S. Breusch & A.R. Pagan (1979), A Simple Test for Heteroscedasticity 
#and Random Coefficient Variation. Econometrica 47, 1287-1294

#para la funcion durbinWatsonTest_nls el paquete car debe estar instalado
durbinWatsonTest_nls <- function(model_nls, alternative="positive",
                                 max.lag=1, method=c("resample","normal"),reps=100){
  if (!inherits(model_nls, "nls"))stop("method only for nls objects")
  residuals<-as.vector(resid(model_nls))
  if (any(is.na(residuals))) stop ("residuals include missing values")
  method <- match.arg(method)
  mu<-fitted(model_nls)
  if(!inherits(model_nls$data,"name")) stop("missign data argument on nls call")
  data_base<-eval(model_nls$data)
  X<-as.matrix(data_base[,colnames(data_base)%in%all.vars(formula(model_nls))[-1]])
  n <- length(residuals)
  r <- dw <-rep(0, max.lag)
  S <- sigma(model_nls)
  den <- sum(residuals^2)
  for (lag in 1:max.lag){
    dw[lag] <- (sum((residuals[(lag+1):n] - residuals[1:(n-lag)])^2))/den
    r[lag] <- (sum(residuals[(lag+1):n]*residuals[1:(n-lag)]))/den
  }
  Y <- if (method == "resample") 
    matrix(sample(residuals, n*reps, replace=TRUE), n, reps) + matrix(mu, n, reps)
  else matrix(rnorm(n*reps, 0, S), n, reps) + matrix(mu, n, reps)
  E <- residuals(lm(Y ~ X - 1))
  DW <- apply(E, 2, car::durbinWatsonTest, max.lag=max.lag)
  if (max.lag == 1) DW <- rbind(DW)
  p <- rep(0, max.lag)
  if (alternative == 'two.sided'){
    for (lag in 1:max.lag) {
      p[lag] <- (sum(dw[lag] < DW[lag,]))/reps
      p[lag] <- 2*(min(p[lag], 1 - p[lag]))
    }
  }
  else if (alternative == "positive"){
    for (lag in 1:max.lag) {
      p[lag] <- (sum(dw[lag] > DW[lag,]))/reps}}
  else {
    for (lag in 1:max.lag) {
      p[lag] <- (sum(dw[lag] < DW[lag,]))/reps}
  }
  result <- list(r=r, dw=dw, p=p, alternative=alternative)
  class(result)<-"durbinWatsonTest_nls"
  return(result)
}

print.durbinWatsonTest_nls <- function(x, ...){
  max.lag <- length(x$dw)
  result <- if (is.null(x$p)) cbind(lag=1:max.lag,Autocorrelation=x$r, "D-W Statistic"=x$dw)
  else cbind(lag=1:max.lag,Autocorrelation = x$r, "D-W Statistic" = x$dw, 
             "p-value"= x$p)
  rownames(result) <- rep("", max.lag)
  print(result)
  cat(paste("Alternative hypothesis: rho", if(max.lag > 1) "[lag]" else "",
            c(" != ", " > ", " < ")[which(x$alternative == c("two.sided", "positive", "negative"))],
            "0\n", sep=""))
  invisible(x)}

bptest_nls <- function(model_nls,studentize=F){
  if(!inherits(model_nls$data,"name")) stop("missign data argument on nls call")
  if (!inherits(model_nls, "nls"))stop("method only for nls objects")
  residuales<-matrix(resid(model_nls),ncol = 1)
  data_base<-eval(model_nls$data)
  Z<-cbind(matrix(1,nrow = length(residuales)),as.matrix(as.matrix(data_base[,colnames(data_base)%in%all.vars(formula(model_nls))[-1]])))
  k <- ncol(Z)
  n <- nrow(Z)
  sigma2 <- sum(residuales^2)/n
  f <- residuales^2/sigma2 -1
  if(studentize){ 
    w <- residuales^2 - sigma2
    aux <- lm.fit(Z, w)    
    bp <- n * sum(aux$fitted.values^2)/sum(w^2)
    method <- "studentized Breusch-Pagan test"
  }
  else{
    f <- residuales^2/sigma2 -1
    aux <- lm.fit(Z, f)
    bp <- 0.5 * sum(aux$fitted.values^2)
    method <- "Breusch-Pagan test"
  }
  names(bp) <- "BP"
  df <- c("df" = aux$rank-1)
  RVAL <- list(statistic = bp,
               data.name = model_nls$data,
               method = method,
               parameter = df,
               p.value= pchisq(bp, df, lower.tail = FALSE))
  class(RVAL) <- "htest"
  return(RVAL)
}
