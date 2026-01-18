
#Time series project 
#Jerome Allier, Axel Pinçon


#Importation of the required packages 
library(zoo)
library(tseries) 
library(fUnitRoots) 
library(ggplot2)
library(xtable)
library(forecast)

#_________________________________ DATA 

# Set working directory to project root (run from src/ directory)
data <- read.csv("../data/arms_truncated_index.csv", sep=';')
data <- data[, -c(3)]
data <- data[-c(1,2,3), ]
data <- data[seq(1,193),]
colnames(data) <- c("date","indice_prod")

dates <- as.yearmon(seq(from=2008, to=2024, by=1/12)) #Set of dates of the observations 


#___________________________________ TIME SERIES  

index_insee <- zoo(data$indice_prod, order.by=dates) #Main time series 
index_insee <- as.numeric(index_insee)

#Plot of the series 
par(mfrow=c(1,1))
ggplot(data = data.frame(time = as.numeric(dates), indice = index_insee), aes(x = time, y = indice)) +
  geom_line() +
  labs(title = 'French arms and ammunition production index from 2008 to 2024',x = 'Year',y = 'Index') +
  scale_x_continuous(breaks = seq(2008, 2024, by = 4),labels = as.character(seq(2008, 2024, by = 4))) +
  theme_minimal()  
#The series is not stationary, it seems that the series follow an increasing trend.

#Linear regression
dates_bis <- as.yearmon(seq(from=2008, to=2024, by=1/12)) 
summary(lm(index_insee ~ dates_bis))


#ACF and PACF plot 
par(mfrow=c(1,2))
pacf(index_insee);acf(index_insee)


#Make the series stationary :

#Differenciated series of order 1 
ind_spread <- diff(index_insee,1) 
ind_spread <- ind_spread - mean(ind_spread) 

par(mfrow=c(1,1))
ggplot(data = data.frame(time = as.numeric(dates)[1:192], indice = ind_spread), aes(x = time, y = indice)) +
  geom_line() +
  labs(title = 'Index series differenciated of order 1',x = 'Year',y = 'Index - lag=1') +
  scale_x_continuous(breaks = seq(2008, 2024, by = 4),labels = as.character(seq(2008, 2024, by = 4))) +
  theme_minimal()  

#ACF and PACF 
par(mfrow=c(1,2))
pacf(ind_spread);acf(ind_spread)


#Tests of stationarity 
pp.test(index_insee)
pp.test(ind_spread)

adf <- adfTest(index_insee,type='ct')
adf_1 <- adfTest(ind_spread)
adf
adf_1



#__________________________ FUNCTIONS 

#Definition of the function Qtests 
Qtests <- function(series, k, fitdf=0) {
  pvals <- apply(matrix(1:k), 1, FUN=function(l) {
    pval <- if (l<=fitdf) NA else Box.test(series, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value
    return(c("lag"=l,"pval"=pval))
  })
  return(t(pvals))
}

#Definition of a function to test every lag between 1 and 24  
#To choose the lag with no autocorrelation between the observations  
series <- ind_spread; kmax <- 24 #; adftype="ct"
adfTest_valid <- function(series, kmax){
  k <- 0
  noautocorr <- 0
  while (noautocorr==0){
    cat(paste0("ADF with ",k," lags: residuals OK? "))
    adf <- adfTest(series, lags=k)
    pvals <- Qtests(adf@test$lm$residuals, 24, fitdf = length(adf@test$lm$coefficients))[,2]
    if (sum(pvals<0.05,na.rm=T)==0) {
      noautocorr <- 1; cat("OK \n")
    } else cat("nope \n")
    k <- k+1
  }
  return(adf)
}


#______________ OTHER FUNCTIONS 

#test if the coefficients are significant 
signif <- function(estim){ 
  coef <- estim$coef
  se <- sqrt(diag(estim$var.coef))
  t <- coef/se
  pval <- (1-pnorm(abs(t)))*2
  return(rbind(coef,se,pval))
}

#To check the validity of models
arimafit <- function(estim){
  adjust <- round(signif(estim),3)
  pvals <- Qtests(estim$residuals,24,length(estim$coef)-1)
  pvals <- matrix(apply(matrix(1:24,nrow=6),2,function(c) round(pvals[c,],3)),nrow=6)
  colnames(pvals) <- rep(c("lag", "pval"),4)
  cat("tests de nullite des coefficients :\n")
  print(adjust)
  cat("\n tests d'absence d'autocorrelation des residus : \n")
  print(pvals)
}

modelchoice <- function(p,q,data=ind_spread, k=24){
  estim <- try(arima(data, c(p,0,q),optim.control=list(maxit=20000)))
  if (class(estim)=="try-error") return(c("p"=p,"q"=q,"arsignif"=NA,"masignif"=NA,"resnocorr"=NA, "ok"=NA))
  arsignif <- if (p==0) NA else signif(estim)[3,p]<=0.05
  masignif <- if (q==0) NA else signif(estim)[3,p+q]<=0.05
  resnocorr <- sum(Qtests(estim$residuals,24,length(estim$coef)-1)[,2]<=0.05,na.rm=T)==0
  checks <- c(arsignif,masignif,resnocorr)
  ok <- as.numeric(sum(checks,na.rm=T)==(3-sum(is.na(checks))))
  return(c("p"=p,"q"=q,"arsignif"=arsignif,"masignif"=masignif,"resnocorr"=resnocorr,"ok"=ok))
}

armamodelchoice <- function(pmax,qmax){
  pqs <- expand.grid(0:pmax,0:qmax)
  t(apply(matrix(1:dim(pqs)[1]),1,function(row) {
    p <- pqs[row,1]; q <- pqs[row,2]
    cat(paste0("Computing ARMA(",p,",",q,") \n"))
    modelchoice(p,q)
  }))
}

#For the minimization of the information criteria 
choice_arma <- function(valid_models){
  mat <- matrix(NA, nrow=length(valid_models), ncol=2)
  colnames(mat) <- paste0(c("AIC","BIC")) 
  for (i in 1:length(valid_models)){
    estim <- arima(ind_spread,valid_models[[i]], include.mean=F) #try to estimate the ARIMA model 
    mat[i,1] <- estim$aic
    mat[i,2] <- BIC(estim)
  }
  return(mat)
}

adj_r2 <- function(model){
  ss_res <- sum(model$residuals^2)
  p <- model$arma[1]
  q <- model$arma[2]
  ss_tot <- sum(index_insee[-c(1:max(p,q))]^2)
  n <- model$nobs-max(p,q)
  adj_r2 <- 1-(ss_res/(n-p-q-1))/(ss_tot/(n-1))
  return(adj_r2)
}



#_______________ VALID MODELS  

#Analyzing the autocorrelograms of the series ind_spread (series indice differenciated of order 1)
#PACF --> pmax=16 
#ACF --> qmax=16 

armamodelchoice(4,5) #For all the pairs (p,q) with p<=16 and q<=16, one can test if the ARMA(p,q) is valid 
#Thus we have all the possible valid model 


#________________________ MINIMIZATION OF THE INFORMATION CRITERIA


#To avoid to run the previous function, here is the result: 
#ARIMA models validated by the previous step of the study (significant coeff, no autocorrelation)
valid_models <- list(c(3,0,1),c(0,0,1),c(3,0,1),c(4,0,1),c(3,0,2),c(1,0,3),c(0,0,4),c(1,0,4),c(2,0,5),c(4,0,5))
p_valid <- list(0,0,1,1,2,3,3,3,4,4)
q_valid <- list(0,1,1,1,2,3,4,4,5,5)


#To choose the best ARIMA model for the series, one calculates the information criteria AIC and BIC for each model previously validated 
info_criteria <- choice_arma(valid_models) #matrix of AIC and BIC for each valid model 
info_criteria 


#ARIMA models that minimize the AIC and BIC 
best_arima_aic <- arima(ind_spread,c(4,0,1),include.mean=F) 
best_arima_bic <- arima(ind_spread,c(0,0,1),include.mean=F) 

best_arima_aic
best_arima_bic


#Test of the null hypothesis H0 of autocorrelation of the residuals 
#Thus one want to reject the null hyp (large p-value) 
Qtests(best_arima_aic$residuals, 24, fitdf=3)
Qtests(best_arima_bic$residuals, 24, fitdf=3)



#______________________ESTIMATION OF THE ARIMA MODELS FOR THE ORIGINAL SERIES 

arima401 <- arima(index_insee,c(4,1,1),include.mean=F)
arima001 <- arima(index_insee,c(0,1,1),include.mean=F)

adfTest(arima001$coef)

Qtests(arima401$residuals, 24, fitdf=3)
Qtests(arima001$residuals, 24, fitdf=3)

#Accuracy 
adj_r2(arima401)
adj_r2(arima001)

accuracy(arima401)


#Plot of the residuals 
dev.off()
par(mfrow=c(1,1))
ggplot(data = data.frame(time = as.numeric(dates)[1:193], indice = arima401$residuals), aes(x = time, y = indice)) +
  geom_line() +
  labs(title = 'ARIMA(4,0,1) Residuals',x = 'Year',y = 'Residuals') +
  scale_x_continuous(breaks = seq(2008, 2024, by = 4),labels = as.character(seq(2008, 2024, by = 4))) +
  theme_minimal()
par(mfrow=c(1,1))
ggplot(data = data.frame(time = as.numeric(dates)[1:193], indice = arima001$residuals), aes(x = time, y = indice)) +
  geom_line() +
  labs(title = 'ARIMA(0,0,1) Residuals',x = 'Year',y = 'Residuals') +
  scale_x_continuous(breaks = seq(2008, 2024, by = 4),labels = as.character(seq(2008, 2024, by = 4))) +
  theme_minimal()

#Tests on the residuals 
Box.test(arima401$residuals, type= "Ljung-Box")
shapiro.test(arima401$residuals) 

qqnorm(arima401$residuals,main = "Normal Q-Q Plot",
       xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
       plot.it = TRUE, datax = FALSE)


#____________________________PREDICTION  

arima1 <- predict(arima401, n.ahead=5) #Prediction of the 5 next values of the series based on the ARIMA(13,1,3) model
arima1$pred

future = forecast(arima401, h = 6,level = c(95))
plot(future) #plot of the predicted values with confidence interval at the 5% level 



