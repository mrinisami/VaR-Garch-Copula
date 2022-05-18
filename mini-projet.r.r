library(tseries)
library(rugarch)
library(copula)
library(e1071)
library(gsl)

#can <- get.hist.quote(instrument = "XIC.TO", start = "2013-04-15", end= "2022-03-11", quote = "AdjClose")
#us <- get.hist.quote(instrument = "XUS.TO", start = "2013-04-15", end= "2022-03-11", quote = "AdjClose")
#can_df <- as.data.frame(can)[1:length(can), 1]
#us_df <- as.data.frame(us)[1:length(us), 1]
#can_r <- diff(log(can_df))[1:2000]
#us_r <- diff(log(us_df))[1:2000]


usdata <- read.csv("xus.csv")
candata <- read.csv("xic.csv")

us_r <- usdata[1:2000, 4]
can_r <- candata[1:2000, 4]
#stats des distributions de rendements

print(c(mean(can_r), sd(can_r), skewness(can_r), kurtosis(can_r)))
print(c(mean(us_r), sd(us_r), skewness(us_r), kurtosis(us_r)))

#sortir les parametres initiaux
test <- ugarchspec(mean.model=list(armaOrder=c(0,0)))
test_ug <- ugarchfit(spec = test, data=can_r)
coef(test_ug)
par_init <- c(0, log(0.00001), log(0.16), log(0.803))
par_init
ctest <- pobs(residuals(test_ug, standardize=TRUE))
#Garch norm can

can_gnorm <- function(param){
  
  n <- length(can_r)
  mu <- param[1]
  alpha1 <- exp(param[2])
  alpha2 <- exp(param[3])
  beta <- exp(param[4])
  tot <- 1e-10
  
  if (alpha2 + beta < 1){
    
    vol_init <- alpha1 / (1 - alpha2 - beta)
    vol <- matrix(0, n, 1)
    vol[1] <- vol_init
    tot <- 0
    #like <- matrix(0, n, 1)
    for (i in 2:n){
      vol[i] <- alpha1 + (alpha2 * (can_r[i - 1] - mu) ^ 2) + beta * vol[i - 1]
      
    }
    like <- -0.5 * log(2 * pi) - 0.5 * log(vol) - (0.5 * (can_r - mu) ^ 2) / vol
    tot <- sum(like)
    #tot <- sum(like)
  }
  print(tot)
  return (-tot)
}


#On sort les parametres du garch pour CAN

val <- optim(par_init, can_gnorm)
cgarch_param <- c(val$par[1], exp(val$par[2:4]))
cgarch_mle <- val$value

#Garch norm US

#test pour param initiaux avec US data
test_us <- ugarchfit(spec=test, data=us_r)
coef(test_us)
us_init <- c(0.00007, log(0.000005), log(0.11), log(0.83))

us_gnorm <- function(param){
  
  n <- length(us_r)
  mu <- param[1]
  alpha1 <- exp(param[2])
  alpha2 <- exp(param[3])
  beta <- exp(param[4])
  
  tot <- 1e-10
  
  if (alpha2 + beta < 1){
    
    vol_init <- alpha1 / (1 - alpha2 - beta)
    vol <- matrix(0, n, 1)
    vol[1] <- vol_init
    tot <- 0
    for (i in 2:n){
      vol[i] <- alpha1 + (alpha2 * (us_r[i - 1] - mu) ^ 2) + beta * vol[i - 1]
    }
    like <- -0.5 * log(2 * pi) - 0.5 * log(vol) - (0.5 * (us_r - mu) ^ 2) / vol
    tot <- sum(like)
  }
  print(tot)
  return(-tot)
}

#Le GARCH avec erreures normales pour US
us_val <- optim(us_init, us_gnorm)
usgarch_param <- c(us_val$par[1], exp(us_val$par[2:4]))
ngarch_mus <- c(cgarch_param[1], usgarch_param[1])
usgarch_mle <- us_val$value

#Student GARCH avec rendements canadiens

can_gt <- function(param){
  
  n <- length(can_r)
  mu <- param[1]
  alpha1 <- (param[2])
  alpha2 <- (param[3])
  beta <- (param[4])
  df <- (param[5])
  tot <- -1e10
  
  if (alpha2 + beta < 1){
    
    vol_init <- alpha1 / (1 - alpha2 - beta)
    vol <- matrix(0, n, 1)
    vol[1] <- vol_init 
    tot <- 0
    e_sd <- matrix(0, n, 1)
    e_sd[1] <- (can_r[1] - mu) / (sqrt(vol_init) * sqrt((df - 2) / df))
    for (i in 2:n){
      vol[i] <- alpha1 + alpha2 * ((e_sd[i - 1] / (sqrt(vol[i - 1]) * sqrt((df - 2) / df))) ^ 2) + beta * vol[i - 1]
      e_sd[i] <- (can_r[i] - mu)
    }
    
    tot <- log(gamma((df + 1) / 2)) - log(gamma(df / 2)) - 0.5 * log(pi * (df - 2) * vol ^ 2) - ((df + 1) / 2) * log(1 + (e_sd ^ 2) / (vol ^ 2 * (df - 2)))
    tot <- sum(tot)
    
  }
  print(tot)
  return (-tot)
}

can_tinit <- c(0.00007, (0.00001), (1.581106e-01), (0.80), (4))
can_tgarchval <- optim(can_tinit, can_gt, method = "L-BFGS-B", lower=c(-1, 0.0001, 0.001, 0.1, 3), upper=c(1, 0.98, 0.98, 0.98, 10))
can_tparam <- can_tgarchval$par

#Student GARCH avec rendements US

us_gt <- function(param){
  
  n <- length(us_r)
  mu <- param[1]
  alpha1 <- (param[2])
  alpha2 <- (param[3])
  beta <- (param[4])
  df <- (param[5])
  tot <- -1e10
  
  if (alpha2 + beta < 1){
    
    vol_init <- alpha1 / (1 - alpha2 - beta)
    vol <- matrix(0, n, 1)
    vol[1] <- vol_init 
    tot <- 0
    e_sd <- matrix(0, n, 1)
    e_sd[1] <- (us_r[1] - mu) / (sqrt(vol_init) * sqrt((df - 2) / df))
    for (i in 2:n){
      vol[i] <- alpha1 + alpha2 * ((e_sd[i - 1] / (sqrt(vol[i - 1]) * sqrt((df - 2) / df))) ^ 2) + beta * vol[i - 1]
      e_sd[i] <- (us_r[i] - mu)
    }
    
    tot <- log(gamma((df + 1) / 2)) - log(gamma(df / 2)) - 0.5 * log(pi * (df - 2) * vol ^ 2) - ((df + 1) / 2) * log(1 + (e_sd ^ 2) / (vol ^ 2 * (df - 2)))
    tot <- sum(tot)
    
  }
  print(tot)
  return (-tot)
}

us_tinit <- c(0, 0, 0.15, 0.835, 6)
us_tgarchval <- optim(us_tinit, us_gt, method = "L-BFGS-B", lower=c(-1, 0.0001, 0.001, 0.1, 3), upper=c(1, 0.98, 0.98, 0.98, 10))
us_tparam <- us_tgarchval$par

#Il faut generer les residus des deux distributions

res <- function(param, data){
  
  n <- length(data)
  mu <- param[1]
  alpha1 <- param[2]
  alpha2 <- param[3]
  beta <- param[4]
  vol_init <- alpha1 / (1 - alpha2 - beta)
  vol <- matrix(0, n, 1)
  vol[1] <- vol_init
  res <- matrix(0, n, 1)
  res_init <- (data[1] - mu) / sqrt(vol_init)
  res[1] <- res_init
  for (i in 2:n){
    vol[i] <- alpha1 + (alpha2 * (data[i - 1] - mu) ^ 2) + beta * vol[i - 1]
    res[i] <- (data[i] - mu) / sqrt(vol[i])
  }
  return (c(res[2:n], vol[length(vol)]))
}

#On couple les residus GARCH normal

us_res <- pobs(res(usgarch_param, us_r))
can_res <- pobs(res(cgarch_param, can_r))
res_data <- cbind(us_res, can_res)
res_data[1:5, 2]

#fonction pour retourner aic et bic

aic_bic <- function(mle, n, p){
  aic <- -2 * mle + 2 * p
  bic <- -2 * mle + log(n) * p
  return(c(aic, bic))
}

#MES COPULES

#Copule Gausienne

my_ngc <- function(param){
  
  theta <- param[1]
  
  cop <- normalCopula(param=theta)
  tot <- 0
  for (i in 1:length(can_res)){
    coplike <- log(dCopula(c(res_data[i, 1], res_data[i, 2]), cop))
    tot <- tot + coplike
  }
  
  print(tot)
  return(-tot)
}

#Donnees copule gaussienne + garch normal
ngc_init <- cor(res_data)[1, 2]
ngc_data <- optim(ngc_init,my_ngc,method="L-BFGS-B",lower=c(-0.98),upper=c(0.98),hessian=TRUE)
ngc_mle <- ngc_data$value
ngc_mul <- aic_bic(ngc_mle, length(res_data), length(ngc_data$par))
ngc_aic <- ngc_mul[1]
ngc_bic <- ngc_mul[2]
ngc_se <- sqrt(diag(solve(ngc_data$hessian)))
ngc_par <- ngc_data$par

#Copule Student

my_nstc <- function(param){
  
  theta <- param[1]
  v <- param[2]
  cop <- tCopula(param=theta, df=v)
  tot <- 0
  for (i in 1:length(can_res)){
    coplike <- log(dCopula(c(res_data[i, 1], res_data[i, 2]), cop))
    tot <- tot + coplike
  }
  
  print(tot)
  return(-tot)
  
}

#Donnees copule student + garch normal
nstc_init <- c(cor(res_data)[1, 2], 4)
nstc_data <- optim(nstc_init,my_nstc,method="L-BFGS-B",lower=c(-0.98, 2.1),upper=c(0.98, 10),hessian=TRUE)
nstc_mle <- nstc_data$value
nstc_mul <- aic_bic(nstc_mle, length(res_data), length(nstc_data$par))
nstc_aic <- nstc_mul[1]
nstc_bic <- nstc_mul[2]
nstc_se <- sqrt(diag(solve(nstc_data$hessian)))
nstc_par <- nstc_data$par

#Copule Frank + GARCH normal

my_nfc <- function(param){
  
  theta <- param[1]

  cop <- frankCopula(param=theta)
  tot <- 0
  for (i in 1:length(can_res)){
    coplike <- log(dCopula(c(res_data[i, 1], res_data[i, 2]), cop))
    tot <- tot + coplike
  }
  
  print(tot)
  return(-tot)
  
}

#Donnees Frank Copula + GARCH normal
nfc_init <- 3
nfc_data <- optim(nfc_init,my_nfc,method="L-BFGS-B",lower=c(-2),upper=c(5) ,hessian=TRUE)
nfc_mle <- nfc_data$value
nfc_mul <- aic_bic(nfc_mle, length(res_data), length(nfc_data$par))
nfc_aic <- nfc_mul[1]
nfc_bic <- nfc_mul[2]
nfc_se <- sqrt(diag(solve(nfc_data$hessian)))
nfc_par <- nfc_data$par

#Copule Clayton + GARCH normal


my_ncc <- function(param){
  
  theta <- param[1]
  
  cop <- claytonCopula(param=theta)
  tot <- 0
  for (i in 1:length(can_res)){
    coplike <- log(dCopula(c(res_data[i, 1], res_data[i, 2]), cop))
    tot <- tot + coplike
  }
  
  print(tot)
  return(-tot)
  
}

#Donnees Clayton + GARCH normal
ncc_init <- 0.8
ncc_data <- optim(ncc_init,my_ncc,method="L-BFGS-B",lower=c(-0.98),upper=c(2) ,hessian=TRUE)
ncc_mle <- ncc_data$value
ncc_mul <- aic_bic(ncc_mle, length(res_data), length(ncc_data$par))
ncc_aic <- ncc_mul[1]
ncc_bic <- ncc_mul[2]
ncc_se <- sqrt(diag(solve(ncc_data$hessian)))
ncc_par <- ncc_data$par

#Gaussian copula VERIFICATION


gaus_c <- normalCopula(param=0.9, dim=2)
g_cop <- fitCopula(gaus_c, data=res_data, method="ml")
g_cop_aic <- -2 * g_cop@loglik + 2 * length(g_cop@estimate)
g_cop_bic <- -2 * g_cop@loglik + log(g_cop@nsample) * length(g_cop@estimate)
g_cop_std <- g_cop@var.est ^ 0.5

#Student copula VERIFICATION

t_c <- tCopula(param=0.1, dim=2, df=2.1)
t_cop <- fitCopula(t_c, data=res_data, method="ml")
t_cop_aic <- -2 * t_cop@loglik + 2 * length(t_cop@estimate)
t_cop_bic <- -2 * t_cop@loglik + log(t_cop@nsample) * length(t_cop@estimate)
t_cop_std <- t_cop@var.est[c(1, 4)] ^ 0.5
 
#Frank copula VERIFICATION

frank_c <- frankCopula(param=0.1, dim=2)
frank_cop <- fitCopula(frank_c, data=res_data, method="ml")
frank_cop_aic <- -2 * frank_cop@loglik + 2 * length(frank_cop@estimate)
frank_cop_bic <- -2 * frank_cop@loglik + log(frank_cop@nsample) * length(frank_cop@estimate)
frank_cop_std <- frank_cop@var.est ^ 0.5

#Clayton copula VERIFICATION

clay_c <- claytonCopula(param=2, dim=2)
clay_cop <- fitCopula(clay_c, data=res_data, method="ml")
clay_cop_aic <- -2 * clay_cop@loglik + 2 * length(clay_cop@estimate)
clay_cop_bic <- -2 * clay_cop@loglik + log(clay_cop@nsample) * length(clay_cop@estimate)
clay_cop_std <- clay_cop@var.est ^ 0.5

##############################
#Student distribution

one_one <- ugarchspec(mean.model=list(armaOrder=c(0,0)), distribution.model="std")
can_test <- ugarchfit(spec=one_one, data=can_r)
tcan_coef <- coef(can_test)

us_test <- ugarchfit(spec=one_one, data=us_r)
tus_coef <- coef(us_test)

c(cgarch_mle, usgarch_mle, 7271, 6845)
c(aic_bic(cgarch_mle, 2000, 4), aic_bic(7271, 2000, 5))
c(aic_bic(cgarch_))
#Couple les residus

tres <- cbind(pobs(can_test@fit$residuals), pobs(us_test@fit$residuals))

#Copule Gausienne

my_tgc <- function(param){
  
  theta <- param[1]
  
  cop <- normalCopula(param=theta)
  tot <- 0
  for (i in 1:length(can_r)){
    coplike <- log(dCopula(c(tres[i, 1], tres[i, 2]), cop))
    tot <- tot + coplike
  }
  
  print(tot)
  return(-tot)
}

#Donnees Gaussienne + GARCH student

tgc_init <- 0.5
tgc_data <- optim(tgc_init,my_tgc,method="L-BFGS-B",lower=c(-0.98),upper=c(0.98) ,hessian=TRUE)
tgc_mle <- tgc_data$value
tgc_mul <- aic_bic(tgc_mle, length(res_data), length(tgc_data$par))
tgc_aic <- tgc_mul[1]
tgc_bic <- tgc_mul[2]
tgc_se <- sqrt(diag(solve(tgc_data$hessian)))
tgc_param <- tgc_data$par

#Copule Student

my_tstc <- function(param){
  
  theta <- param[1]
  v <- param[2]
  cop <- tCopula(param=theta, df=v)
  tot <- 0
  for (i in 1:length(can_r)){
    coplike <- log(dCopula(c(tres[i, 1], tres[i, 2]), cop))
    tot <- tot + coplike
  }
  
  print(tot)
  return(-tot)
  
}

#Donnees copule student + garch student
tstc_init <- c(cor(tres)[1, 2], 4)
tstc_data <- optim(tstc_init,my_tstc,method="L-BFGS-B",lower=c(-0.98, 2.1),upper=c(0.98, 10),hessian=TRUE)
tstc_mle <- tstc_data$value
tstc_mul <- aic_bic(nstc_mle, length(res_data), length(tstc_data$par))
tstc_aic <- nstc_mul[1]
tstc_bic <- nstc_mul[2]
tstc_se <- sqrt(diag(solve(tstc_data$hessian)))
tstc_param <- tstc_data$par

#Copule Frank + GARCH student

my_tfc <- function(param){
  
  theta <- param[1]
  
  cop <- frankCopula(param=theta)
  tot <- 0
  for (i in 1:length(can_r)){
    coplike <- log(dCopula(c(tres[i, 1], tres[i, 2]), cop))
    tot <- tot + coplike
  }
  
  print(tot)
  return(-tot)
  
}

#Donnees Frank Copula + GARCH student
tfc_init <- 3
tfc_data <- optim(tfc_init,my_tfc,method="L-BFGS-B",lower=c(-2),upper=c(5) ,hessian=TRUE)
tfc_mle <- tfc_data$value
tfc_mul <- aic_bic(tfc_mle, length(res_data), length(tfc_data$par))
tfc_aic <- tfc_mul[1]
tfc_bic <- tfc_mul[2]
tfc_se <- sqrt(diag(solve(tfc_data$hessian)))
tfc_param <- tfc_data$par

#Donnees Clayton Copula + GARCH student

my_tcc <- function(param){
  
  theta <- param[1]
  
  cop <- claytonCopula(param=theta)
  tot <- 0
  for (i in 1:length(can_r)){
    coplike <- log(dCopula(c(tres[i, 1], tres[i, 2]), cop))
    tot <- tot + coplike
  }
  
  print(tot)
  return(-tot)
  
}

#Donnees Clayton + GARCH student
tcc_init <- 0.8
tcc_data <- optim(tcc_init,my_tcc,method="L-BFGS-B",lower=c(-0.98),upper=c(2), hessian=TRUE)
tcc_mle <- tcc_data$value
tcc_mul <- aic_bic(tcc_mle, length(res_data), length(tcc_data$par))
tcc_aic <- tcc_mul[1]
tcc_bic <- tcc_mul[2]
tcc_se <- sqrt(diag(solve(tcc_data$hessian)))
tcc_param <- tcc_data$par

#Creation of data frame for summary of models

truelog_like <- c(ngc_mle, nstc_mle, nfc_mle, ncc_mle, tgc_mle, tstc_mle, tfc_mle, tcc_mle)
truest_err <- c(ngc_se, nstc_se[1], nfc_se, ncc_se, tgc_se, tstc_se[1], tfc_se, tcc_se)
truestu_err2 <- c("N.A", nstc_se[2], "N.A", "N.A", "N.A",  tstc_se[2], "N.A", "N.A")
trueaic <- c(ngc_aic, nstc_aic, nfc_aic, ncc_aic, tgc_aic, tstc_aic, tfc_aic, tcc_aic)
truebic <- c(ngc_bic, nstc_bic, nfc_bic, ncc_bic, tgc_bic, tstc_bic, tfc_bic, tcc_bic)
truenames <- c('Gaus - Normal', "Student - Normal", "Frank - Normal", "Clayton - Normal", "Gaus - Student", "Student - Student", 'Frank - Student', "Clayton - Student")

#Le tableau en question
cbind(truenames, truelog_like, truest_err, truestu_err2, trueaic, truebic)

#Generer 634 predictions avec GARCH norm

nforecast <- function(param, data, n){
  
  vols <- res(param, data)
  last_vol <- vols[length(vols)]
  mu <- param[1]
  alpha1 <- param[2]
  alpha2 <- param[3]
  beta <- param[4]
  last_ret <- data[length(data)] - mu
  preds <- matrix(0, n, 1)
  preds[1] <- alpha1 + alpha2 * last_ret ^ 2 + beta * last_vol
  for (i in 2:n){
    preds[i] <- alpha1 + (alpha2 + beta) * preds[i - 1]
  }
  return(preds[2:n])
  
}

ncan_preds <- nforecast(cgarch_param, can_r, 235)
nus_preds <- nforecast(usgarch_param, us_r, 235)

#Generer pred avec GARCH student

tforecast <- function(param, data, n){
  
  vols <- res(param, data)
  last_vol <- vols[length(vols)]
  mu <- param[1]
  alpha1 <- param[2]
  alpha2 <- param[3]
  beta <- param[4]
  last_ret <- data[length(data)] - mu
  preds <- matrix(0, n, 1)
  preds[1] <- alpha1 + alpha2 * last_ret ^ 2 + beta * last_vol
  for (i in 2:n){
    preds[i] <- alpha1 + (alpha2 + beta) * preds[i - 1]
  }
  return(preds[2:n])
  
}

tcgarch_param <- tcan_coef
tusgarch_param <- tus_coef
tcan_preds <- tforecast(tcgarch_param, can_r, 235)
tus_preds <- tforecast(tusgarch_param, us_r, 235)
tgarch_mus <- c(tcgarch_param[1], tusgarch_param[1])

#Simulate data from Gaussian Copula + norm GARCH

ncop_forecast <- function(param, n){
  
  norm <- param[1]
  std <- param[2]
  k <- param[3]
  frank <- param[4]
  clayton <- param[5]
  ncop <- normalCopula(norm)
  rncop <- rCopula(n, ncop)
  tcop <- tCopula(std, df=k)
  rtcop <- rCopula(n, tcop)
  fcop <- frankCopula(frank)
  rfcop <- rCopula(n, fcop)
  ccop <- claytonCopula(clayton)
  rccop <- rCopula(n, ccop)
  
  return(cbind(rncop, rtcop, rfcop, rccop))
  
}

#VaRs computed for all NORMAL GARCH models at once

cmp_var <- function(param, n){+
  
  ncan_preds <- nforecast(cgarch_param, can_r, 235)
  nus_preds <- nforecast(usgarch_param, us_r, 235)
  ptf_mus <- 0.5 * cgarch_param[1] + 0.5 * cgarch_param[2]
  quant <- matrix(0, 234, 8)
  vars <- matrix(0, 234, 8)
  for (i in 1:234){
    data <- ncop_forecast(param, n)
    gaus <- data[1:n, 1:2]
    gaus <- qnorm(gaus)
    std <- data[1:n, 3:4]
    std <- qnorm(std)
    frank <- data[1:n, 5:6]
    frank <- qnorm(frank)
    clayton <- data[1:n, 7:8]
    clayton <- qnorm(clayton)
    
    #gaussian copula simulated  errors, their quantiles and the VaR
    cov <- (2 / pi * asin(param[1])) * sqrt(ncan_preds) * sqrt(nus_preds)
    var_ptf <- ncan_preds * 0.5 ^ 2 + nus_preds * 0.5 ^ 2 + 2 * (0.5 ^ 2 * cov)
    comb_err <- gaus[1:n, 1] * 0.5 + gaus[1:n, 2] * 0.5
    quant[i, 1:2] <- quantile(comb_err, probs=c(0.95, 0.99))
    vars[i, 1:2] <- ptf_mus - quant[i, 1:2] * sqrt(var_ptf[i])
    
    #student copula simulated  errors, their quantiles and the VaR
    corr <- (2 / pi * asin(param[2]))  * sqrt(ncan_preds) * sqrt(nus_preds)
    var_ptf <- ncan_preds * 0.5 ^ 2 + nus_preds * 0.5 ^ 2 + 2 * (0.5 ^ 2 * cov)
    comb_err <- std[1:n, 1] * 0.5 + std[1:n, 2] * 0.5
    quant[i, 3:4] <- quantile(comb_err, probs=c(0.95, 0.99))
    vars[i, 3:4] <- ptf_mus - quant[i, 3:4] * sqrt(var_ptf[i])
    
    #Frank copula simulated  errors, their quantiles and the VaR
    corr <- 1 - (4 / param[4]) * (1 - debye_1(param[4]))  * sqrt(ncan_preds) * sqrt(nus_preds)
    var_ptf <- ncan_preds * 0.5 ^ 2 + nus_preds * 0.5 ^ 2 + 2 * (0.5 ^ 2 * cov)
    comb_err <- frank[1:n, 1] * 0.5 + frank[1:n, 2] * 0.5
    quant[i, 5:6] <- quantile(comb_err, probs=c(0.95, 0.99))
    vars[i, 5:6] <- ptf_mus - quant[i, 5:6] * sqrt(var_ptf[i])
    
    #Clayton copula simulated errors, their quantiles and the VaR
    corr <- param[5] / (2 + param[5]) * sqrt(ncan_preds) * sqrt(nus_preds)
    var_ptf <- ncan_preds * 0.5 ^ 2 + nus_preds * 0.5 ^ 2 + 2 * (0.5 ^ 2 * cov)
    comb_err <- clayton[1:n, 1] * 0.5 + clayton[1:n, 2] * 0.5
    quant[i, 7:8] <- quantile(comb_err, probs=c(0.95, 0.99))
    vars[i, 7:8] <- ptf_mus - quant[i, 7:8] * sqrt(var_ptf[i]) 
  }
  return (vars)
}

#Data for the VaRs using NORMAL GARCH
var_norm <- cmp_var(c(ngc_par, nstc_par, nfc_par, ncc_par), 1000)

#VaRs function for student distributed marginals

cmp_vart <- function(param, n){
  
  tcan_preds <- tforecast(tcgarch_param, can_r, 235)
  tus_preds <- tforecast(tusgarch_param, us_r, 235)
  us_v <- tusgarch_param[5]
  can_v <- tcgarch_param[5]
  ptf_mus <- 0.5 * cgarch_param[1] + 0.5 * cgarch_param[2]
  vars <- matrix(0, 234,8)
  quant <- matrix(0, 234, 8)
  for (i in 1:234){
    data <- ncop_forecast(param, n)
    gaus <- data[1:n, 1:2]
    gaus_c <- qt(gaus[1:n, 1], can_v)
    gaus_u <- qt(gaus[1:n, 2], us_v)
    std <- data[1:n, 3:4]
    std_c <- qt(std[1:n, 1], can_v)
    std_u <- qt(std[1:n, 2], us_v)
    frank <- data[1:n, 5:6]
    frank_c <- qt(frank[1:n, 1], can_v)
    frank_u <- qt(frank[1:n, 2], us_v)
    clayton <- data[1:n, 7:8]
    clayton_c <- qt(clayton[1:n, 1], can_v)
    clayton_u <- qt(clayton[1:n, 2], us_v)
    
    #gaussian copula simulated  errors, their quantiles and the VaR
    cov <- (2 / pi * asin(param[1])) * sqrt(tcan_preds) * sqrt(tus_preds)
    var_ptf <- tcan_preds * 0.5 ^ 2 + tus_preds * 0.5 ^ 2 + 2 * (0.5 ^ 2 * cov)
    comb_err <- gaus_c * 0.5 + gaus_u * 0.5
    quant[i, 1:2] <- quantile(comb_err, probs=c(0.95, 0.99))
    vars[i, 1:2] <- ptf_mus - quant[i, 1:2] * sqrt(var_ptf[i])
    
    #student copula simulated  errors, the quantiles and the VaR
    corr <- (2 / pi * asin(param[2]))  * sqrt(tcan_preds) * sqrt(tus_preds)
    var_ptf <- tcan_preds * 0.5 ^ 2 + tus_preds * 0.5 ^ 2 + 2 * (0.5 ^ 2 * cov)
    comb_err <- std_c * 0.5 + std_u * 0.5
    quant[i, 3:4] <- quantile(comb_err, probs=c(0.95, 0.99))
    vars[i, 3:4] <- ptf_mus - quant[i, 3:4] * sqrt(var_ptf[i])
    
    #Frank copula simulated  errors, their quantiles and the VaR
    corr <- 1 - (4 / param[4]) * (1 - debye_1(param[4]))  * sqrt(tcan_preds) * sqrt(tus_preds)
    var_ptf <- tcan_preds * 0.5 ^ 2 + tus_preds * 0.5 ^ 2 + 2 * (0.5 ^ 2 * cov)
    comb_err <- frank_c * 0.5 + frank_u * 0.5
    quant[i, 5:6] <- quantile(comb_err, probs=c(0.95, 0.99))
    vars[i, 5:6] <- ptf_mus - quant[i, 5:6] * sqrt(var_ptf[i])
    
    #Clayton copula simulated errors, their quantiles and the VaR
    corr <- param[5] / (2 + param[5]) * sqrt(tcan_preds) * sqrt(tus_preds)
    var_ptf <- tcan_preds * 0.5 ^ 2 + tus_preds * 0.5 ^ 2 + 2 * (0.5 ^ 2 * cov)
    comb_err <- clayton_c * 0.5 + clayton_u * 0.5
    quant[i, 7:8] <- quantile(comb_err, probs=c(0.95, 0.99))
    vars[i, 7:8] <- ptf_mus - quant[i, 7:8] * sqrt(var_ptf[i]) 
  }
  return (vars)
}

var_t <- cmp_vart(c(tgc_param, tstc_param, tfc_param, tcc_param), 1000)

#Graphs + hits


graph_var <- function(vars, title){
  us_realized <- usdata[2001:2234, 4]
  can_realized <- candata[2001:2234, 4]
  ptf <- us_realized * 0.5 + can_realized * 0.5
  nb_hits <- length(which(vars > ptf))
  dates <- seq(1, 234)
  plot(dates, ptf, xlab="Acting as dates", ylab="Returns")
  lines(vars)
  text(35, 0.015, paste0("Number of hits :",nb_hits))
  title(title)

  
  return(nb_hits)
}

#Generating the graphs : order : garch normal and garch student alongside same copula
#95% GARCH NORMAL THEN GARCH STUDENT BOTH W/ GAUSS COP
graph_var(var_norm[1:234, 1], "VaR 95% w/ Garch normal and gaussian copula")
graph_var(var_t[1:234, 1], "VaR 95% w/ Garch student and gaussian copula")

#99% W/ GAUSS COP
graph_var(var_norm[1:234, 2], "VaR 99% w/ Garch normal and gaussian copula")
graph_var(var_t[1:234, 2], "VaR 99% w/ Garch student and gaussian copula")

#95% W/ STUDENT COP
graph_var(var_norm[1:234, 3], "VaR 95% w/ Garch normal and student copula")
graph_var(var_t[1:234, 3], "VaR 95% w/ Garch student and student copula")

#99% W/ STUDENT COP
graph_var(var_norm[1:234, 4], "VaR 99% w/ Garch normal and student copula")
graph_var(var_t[1:234, 4], "VaR 99% w/ Garch student and student copula")

#95% W/ FRANK COP
graph_var(var_norm[1:234, 5], "VaR 95% w/ Garch normal and frank copula")
graph_var(var_t[1:234, 5], "VaR 95% w/ Garch student and frank copula")

#99% W/ FRANK COP
graph_var(var_norm[1:234, 6], "VaR 99% w/ Garch normal and frank copula")
graph_var(var_t[1:234, 6], "VaR 99% w/ Garch student and frank copula")

#95% W/ CLAYTON COP
graph_var(var_norm[1:234, 7], "VaR 95% w/ Garch normal and clayton copula")
graph_var(var_t[1:234, 7], "VaR 95% w/ Garch student and clayton copula")

 #99% W/ CLAYTON COP
graph_var(var_norm[1:234, 8], "VaR 99% w/ Garch normal and clayton copula")
graph_var(var_t[1:234, 8], "VaR 99% w/ Garch student and clayton copula")





