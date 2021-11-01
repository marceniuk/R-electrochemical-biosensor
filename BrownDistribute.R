# ATTENTION!!! R version 3.6.1 + Julia 1.5.3
# gamma distributed delays model of electrochemical biosensor


# !!!!! ATTENTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!! IT IS NECESSARY TO CREATE DATAFRAME USING CORRECT LOCALIZATION OF THE DATA!!!!!
#!!!! in the example below we use XLSX data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!! in case of CSV use read_csv !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
require(xlsx)
#install.packages('xlsx', type = 'source', INSTALL_opts='--no-multiarch')
require(dplyr)
library(ggplot2)
library(reshape2) # library for reshaping data (tall-narrow <-> short-wide)
library(deSolve) # library for solving differential equations
library(minpack.lm) # library for least squares fit using levenberg-marquart algorithm

#!!!!!!!!! it should be correct localization of input data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
dname<-"C:/My_doc/KlosWitkowska/ElectroChemicalBiosensorBSA/Plots/gammadistribution_multiple_volumes"
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
dname1<-"CLEA_substract"
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
dir.create(file.path(dname,dname1), showWarnings = FALSE)
setwd(file.path(dname,dname1))

my_theme <- theme_bw() +  theme(axis.text=element_text(size=14),axis.title=element_text(size=16),legend.title=element_text(size=16),legend.text = element_text(size=16)) 


# ingest experimental data of conductivity


df <- read.xlsx("C:/My_doc/KlosWitkowska/ElectroChemicalBiosensorBSA/TimeSerieses/CLEAandSubstrate.xlsx",sheetName = "Sheet1",stringsAsFactors=F)
df <- df[which(df$time_sec>=250),]
tmpvec <- df$substrate_0dot1
df$substrate_0dot1 <- df$substrate_0dot1 - c(rep(tmpvec[1],length(tmpvec))) 
tmpvec <- df$substrate_0dot3
df$substrate_0dot3 <- df$substrate_0dot3 - c(rep(tmpvec[1],length(tmpvec))) 
tmpvec <- df$substrate_0dot9
df$substrate_0dot9 <- df$substrate_0dot9 - c(rep(tmpvec[1],length(tmpvec))) 
tmpvec <- df$substrate_1dot5
df$substrate_1dot5 <- df$substrate_1dot5 - c(rep(tmpvec[1],length(tmpvec))) 

# plot data
tmp=melt(df,id.vars=c("time_sec"),variable.name="conductivity",value.name="muS/cm")
png(paste(paste("substrate_0dot1"),"png", sep = "."), width = 800, height = 600)
p<-ggplot(data=tmp[which(tmp$conductivity=="substrate_0dot1"),],aes(x=time_sec,y=`muS/cm`,color=conductivity))+geom_point(size=1)+my_theme
plot(p)
dev.off()
png(paste(paste("substrate_0dot3"),"png", sep = "."), width = 800, height = 600)
p<-ggplot(data=tmp[which(tmp$conductivity=="substrate_0dot3"),],aes(x=time_sec,y=`muS/cm`,color=conductivity))+geom_point(size=1)+my_theme
plot(p)
dev.off()
png(paste(paste("substrate_0dot9"),"png", sep = "."), width = 800, height = 600)
p<-ggplot(data=tmp[which(tmp$conductivity=="substrate_0dot9"),],aes(x=time_sec,y=`muS/cm`,color=conductivity))+geom_point(size=1)+my_theme
plot(p)
dev.off()
png(paste(paste("substrate_1dot5"),"png", sep = "."), width = 800, height = 600)
p<-ggplot(data=tmp[which(tmp$conductivity=="substrate_1dot5"),],aes(x=time_sec,y=`muS/cm`,color=conductivity))+geom_point(size=1)+my_theme
plot(p)
dev.off()


################################
# Julia installation required ##
################################

# install.packages("diffeqr")
# devtools::install_github('SciML/diffeqr', build_vignettes=T)

# diffeqr::diffeq_setup()  # installation of Julia

library(JuliaCall)
library(diffeqr)

de <- diffeqr::diffeq_setup()

f_dct <- JuliaCall::julia_eval("function f_dct(du, u, h, p, t)
  k_d=p[1]
  a=p[2]
  m=p[3]
  tau_m=p[4]
  n_S_0=p[5]
  n_E_0=p[6]
  n_P_0=p[7]
  tau=tau_m+(m+1)/a
  du[1]= -k_d*u[2]*u[1]
  du[2]=-k_d*u[2]*u[1] + k_d*h(p,t-tau)[2]*h(p,t-tau)[1]
  du[3]=k_d*h(p,t-tau)[2]*h(p,t-tau)[1]
end")  
h <- JuliaCall::julia_eval("function h(p, t)
  [p[5],p[6],p[7]]
end")
psi <- JuliaCall::julia_eval("function psi(a, m, tau_m, tau)
  if tau <= tau_m
    res = 0
  else 
    res = (a^(m+1)/gamma(m+1))*(tau-tau_m)^m * exp(-a*(tau-tau_m))
  end
  res
end")

# S_0 = 0.3

u0 <- c(0.1,1,0)
tspan <- c(max(0,min(df$time_sec)),max(0,max(df$time_sec)))
abstol <- 1e-2
reltol <- 1e-2
saveat <- round(seq(from=min(df$time),to=max(df$time),length.out=500))
saveat <- unique(saveat)
p <- c(0.05,1,20,5,u0[1:3])
constant_lags <- c(100.0)
JuliaCall::julia_assign("u0", u0)
JuliaCall::julia_assign("tspan", tspan)
JuliaCall::julia_assign("constant_lags", constant_lags)
JuliaCall::julia_assign("p", p)
prob <- JuliaCall::julia_eval("using DifferentialEquations; using SpecialFunctions; DDEProblem(f_dct, u0, h, tspan, p)")
system.time(sol <- de$solve(prob,abstol=abstol,reltol=reltol, saveat=saveat, de$MethodOfSteps(de$Tsit5())))
udf <- as.data.frame(cbind(sol$t,t(sapply(sol$u,identity))))
colnames(udf) <- c("t","S","E","P")

library(plotly)

font <- list(
  family = "arial",
  size = 14,
  color = 'black')

fig <- plotly::plot_ly(udf, x = sol$t, y = ~S, type = 'scatter', mode = 'lines') 
fig <- fig %>% layout(title=list(text="Example",y=0.95), font=font)
orca(fig, file = paste(paste("Sexample"),"png", sep = "."))
fig <- plotly::plot_ly(udf, x = sol$t, y = ~E, type = 'scatter', mode = 'lines') 
fig <- fig %>% layout(title=list(text="Example",y=0.95), font=font)
orca(fig, file = paste(paste("Eexample"),"png", sep = "."))
fig <- plotly::plot_ly(udf, x = sol$t, y = ~P, type = 'scatter', mode = 'lines') 
fig <- fig %>% layout(title=list(text="Example",y=0.95), font=font)
orca(fig, file = paste(paste("Pexample"),"png", sep = "."))
fig <- plotly::plot_ly(udf, x = sol$t, y = ~S, type = 'scatter', mode = 'lines') %>% plotly::add_trace(y = ~E) %>% plotly::add_trace(y = ~P)
fig <- fig %>% layout(title=list(text="Example",y=0.95), font=font)
orca(fig, file = paste(paste("SEPexample"),"png", sep = "."))


#### Parameter identification

ssq=function(parms){
  de <- diffeqr::diffeq_setup()
  print(parms)
  
  f <- JuliaCall::julia_eval("function f(du, u, h, p, t)
  k_d=p[1]
  a=p[2]
  m=p[3]
  tau_m=p[4]
  n_S_0=p[5]
  n_E_0=p[6]
  n_P_0=p[7]
  n=750
  rmse=sqrt((m+1)/a^2)
  lags=range(1/n, length=n, stop=(tau_m+(m+1)/a+3*rmse)*(1-1/n))
  du[1]= -k_d*u[2]*u[1]
  du[2]=-k_d*u[2]*u[1] + k_d*(tau_m+(m+1)/a+3*rmse)*(1/n)*sum([psi(a,m,tau_m,tau)*h(p,t-tau)[2]*h(p,t-tau)[1] for tau in lags])
  du[3]=k_d*(tau_m+(m+1)/a+3*rmse)*(1/n)*sum([psi(a,m,tau_m,tau)*h(p,t-tau)[2]*h(p,t-tau)[1] for tau in lags])
end")  
  h <- JuliaCall::julia_eval("function h(p, t)
  [p[5],p[6],p[7]]
end")
  psi <- JuliaCall::julia_eval("function psi(a, m, tau_m, tau)
  if tau <= tau_m
    res = 0
  else 
    res = (a^(m+1)/gamma(m+1))*(tau-tau_m)^m * exp(-a*(tau-tau_m))
  end
  res
end")
  
  # solve DDEs for a given set of parameters
  ################ for 0.1 #################################
  u0 <- c(0.1,parms[6],df$substrate_0dot1[1])
  tspan <- c(max(0,min(df$time_sec)),max(0,max(df$time_sec)))
  abstol <- 1e-2
  reltol <- 1e-2
  saveat <- round(seq(from=min(df$time),to=max(df$time),length.out=500))
  saveat <- unique(saveat)
  parms[5]<-0.1
  p <- parms
  constant_lags <- c(100.0)
  JuliaCall::julia_assign("u0", u0)
  JuliaCall::julia_assign("tspan", tspan)
  JuliaCall::julia_assign("constant_lags", constant_lags)
  JuliaCall::julia_assign("p", p)
  #prob <- JuliaCall::julia_eval("using SpecialFunctions; DDEProblem(f, u0, h, tspan, constant_lags = constant_lags, p)")
  prob <- JuliaCall::julia_eval("using SpecialFunctions; DDEProblem(f, u0, h, tspan, p)")
  system.time(sol <- de$solve(prob,abstol=abstol,reltol=reltol, saveat=saveat, de$MethodOfSteps(de$Tsit5())))
  udf <- as.data.frame(cbind(sol$t,t(sapply(sol$u,identity))))
  colnames(udf) <- c("t","S","E","P")
  
  # Filter data that contains time points where data is available
  
  outdf=udf[udf$t %in% df$time_sec,]
  
  # Evaluate predicted vs experimental residual
  preddf=outdf$P
  expdf=df$substrate_0dot1
  
  #print(parms[8]*preddf - expdf)
  #ssqres_0dot1=parms[8]*sqrt(preddf)-expdf
  ssqres_0dot1=parms[8]*preddf - parms[9]*'^'(preddf,3/2) - expdf
  ################ for 0.3 #################################
  u0 <- c(0.3,parms[6],df$substrate_0dot3[1])
  tspan <- c(max(0,min(df$time_sec)),max(0,max(df$time_sec)))
  abstol <- 1e-2
  reltol <- 1e-2
  saveat <- round(seq(from=min(df$time),to=max(df$time),length.out=500))
  saveat <- unique(saveat)
  parms[5]<-0.3
  p <- parms
  constant_lags <- c(100.0)
  JuliaCall::julia_assign("u0", u0)
  JuliaCall::julia_assign("tspan", tspan)
  JuliaCall::julia_assign("constant_lags", constant_lags)
  JuliaCall::julia_assign("p", p)
  #prob <- JuliaCall::julia_eval("using SpecialFunctions; DDEProblem(f, u0, h, tspan, constant_lags = constant_lags, p)")
  prob <- JuliaCall::julia_eval("using SpecialFunctions; DDEProblem(f, u0, h, tspan, p)")
  system.time(sol <- de$solve(prob,abstol=abstol,reltol=reltol, saveat=saveat, de$MethodOfSteps(de$Tsit5())))
  udf <- as.data.frame(cbind(sol$t,t(sapply(sol$u,identity))))
  colnames(udf) <- c("t","S","E","P")
  
  # Filter data that contains time points where data is available
  
  outdf=udf[udf$t %in% df$time_sec,]
  
  # Evaluate predicted vs experimental residual
  preddf=outdf$P
  expdf=df$substrate_0dot3
  
  #print(parms[8]*preddf - expdf)
  
  #ssqres_0dot3=parms[8]*sqrt(preddf)-expdf
  ssqres_0dot3=parms[8]*preddf - parms[9]*'^'(preddf,3/2) - expdf
  ################ for 0.9 #################################
  u0 <- c(0.9,parms[6],df$substrate_0dot3[1])
  tspan <- c(max(0,min(df$time_sec)),max(0,max(df$time_sec)))
  abstol <- 1e-2
  reltol <- 1e-2
  saveat <- round(seq(from=min(df$time),to=max(df$time),length.out=500))
  saveat <- unique(saveat)
  parms[5]<-0.9
  p <- parms
  constant_lags <- c(100.0)
  JuliaCall::julia_assign("u0", u0)
  JuliaCall::julia_assign("tspan", tspan)
  JuliaCall::julia_assign("constant_lags", constant_lags)
  JuliaCall::julia_assign("p", p)
  #prob <- JuliaCall::julia_eval("using SpecialFunctions; DDEProblem(f, u0, h, tspan, constant_lags = constant_lags, p)")
  prob <- JuliaCall::julia_eval("using SpecialFunctions; DDEProblem(f, u0, h, tspan, p)")
  system.time(sol <- de$solve(prob,abstol=abstol,reltol=reltol, saveat=saveat, de$MethodOfSteps(de$Tsit5())))
  udf <- as.data.frame(cbind(sol$t,t(sapply(sol$u,identity))))
  colnames(udf) <- c("t","S","E","P")
  
  # Filter data that contains time points where data is available
  
  outdf=udf[udf$t %in% df$time_sec,]
  
  # Evaluate predicted vs experimental residual
  preddf=outdf$P
  expdf=df$substrate_0dot9
  
  #print(parms[8]*preddf - expdf)
  
  #ssqres_0dot9=parms[8]*sqrt(preddf)-expdf
  ssqres_0dot9=parms[8]*preddf - parms[9]*'^'(preddf,3/2) - expdf
  
  ################ for 1.5 #################################
  u0 <- c(1.5,parms[6],df$substrate_1dot5[1])
  tspan <- c(max(0,min(df$time_sec)),max(0,max(df$time_sec)))
  abstol <- 1e-2
  reltol <- 1e-2
  saveat <- round(seq(from=min(df$time),to=max(df$time),length.out=500))
  saveat <- unique(saveat)
  parms[5]<-1.5
  p <- parms
  constant_lags <- c(100.0)
  JuliaCall::julia_assign("u0", u0)
  JuliaCall::julia_assign("tspan", tspan)
  JuliaCall::julia_assign("constant_lags", constant_lags)
  JuliaCall::julia_assign("p", p)
  #prob <- JuliaCall::julia_eval("using SpecialFunctions; DDEProblem(f, u0, h, tspan, constant_lags = constant_lags, p)")
  prob <- JuliaCall::julia_eval("using SpecialFunctions; DDEProblem(f, u0, h, tspan, p)")
  system.time(sol <- de$solve(prob,abstol=abstol,reltol=reltol, saveat=saveat, de$MethodOfSteps(de$Tsit5())))
  udf <- as.data.frame(cbind(sol$t,t(sapply(sol$u,identity))))
  colnames(udf) <- c("t","S","E","P")
  
  # Filter data that contains time points where data is available
  
  outdf=udf[udf$t %in% df$time_sec,]
  
  # Evaluate predicted vs experimental residual
  preddf=outdf$P
  expdf=df$substrate1dot5
  
  #print(parms[8]*preddf - expdf)
  
  #ssqres_1dot5=parms[8]*sqrt(preddf)-expdf
  ssqres_1dot5=parms[8]*preddf - parms[9]*'^'(preddf,3/2) - expdf
  
  ##########################################
  # return predicted vs experimental residual
  ssqres <- c(ssqres_0dot1,ssqres_0dot3,ssqres_0dot9,ssqres_1dot5)
  #print(ssqres)
  
  print(mean(abs(ssqres)))
  
  plot(x=udf$t,y=udf$P)
  
  return(ssqres)
}

# parameter fitting 
# initial guess for parameters

parms <- c(0.04, 1,20,5,0.5,4.1,0,6000,50)
lower <- c(rep(1.e-10,4),0.5,4.1,0, 1.e-10,1.e-10)
upper <- c(1,rep(1000,3),0.5,4.1,0,1.e6,1.e6)

################ fitting using levenberg marquart algorithm
if(FALSE){
  fitval=nls.lm(par=parms,lower=lower,upper=upper,fn=ssq, control = nls.lm.control(maxiter = 5,nprint = 2))
  #fitval=nls.lm(par=parms,lower=lower,upper=upper,fn=ssq)
  
  # Summary of fit
  #summary(fitval)
  
  # Estimated parameter
  parest=as.list(coef(fitval))
  parest
  
  
  # degrees of freedom: # data points - # parameters
  dof=3*nrow(df)-2
  dof
  
  # mean error
  ms=sqrt(deviance(fitval)/dof)
  ms
}

########## # parameter fitting using COBYLA #######################
# M. J. D. Powell, ``A direct search optimization method that models the objective and constraint functions by linear interpolation,'' in Advances in Optimization and Numerical Analysis, eds. S. Gomez and J.-P. Hennart (Kluwer Academic: Dordrecht, 1994), p. 51-67.
if(TRUE){
  library(nloptr)
  fn.sep <- function(x) {
    sqrt(sum(ssq(x)^2))
  }
  hin.sep <- function(x) {
    h <- numeric(18)
    h[1:9] <- upper - x
    h[10:18] <- x - lower
    return(h)
  }
  fitval.cobyla <- cobyla(parms, fn.sep, hin = hin.sep,
                          nl.info = TRUE, control = list(xtol_rel = 1e-3,maxeval = 100))
  
  parest <- fitval.cobyla$par
  
  fitval.cobyla$value
  
  # degrees of freedom: # data points - # parameters
  dof=3*nrow(df)-2
  dof
  
  # mean error
  ms=sqrt(deviance(fitval.cobyla$value)/dof)
  ms
  
}
##########################################
# plot estimated model ###################
##########################################
saveat <- round(seq(from=min(df$time_sec),to=max(df$time_sec),length.out=1000))
saveat <- unique(saveat)
tspan <- c(max(0,min(df$time_sec)),max(0,max(df$time_sec)))#c(max(0,min(df$time_sec)),max(0,max(df$time_sec)))
##########################################
# S_0=0.9 ################################
##########################################
u0 <- unlist(parest[5:7]) 
u0[1] <- 0.9
#tspan <- c(max(0,min(df$time_sec)),max(0,max(df$time_sec)))
p <- parest
constant_lags <- c(10.0)
JuliaCall::julia_assign("u0", u0)
JuliaCall::julia_assign("tspan", tspan)
JuliaCall::julia_assign("constant_lags", constant_lags)
JuliaCall::julia_assign("p", p)
prob <- JuliaCall::julia_eval("DDEProblem(f, u0, h, tspan, constant_lags = constant_lags, p)")
prob_dct <- JuliaCall::julia_eval("DDEProblem(f_dct, u0, h, tspan, constant_lags = constant_lags, p)")
system.time(sol <- de$solve(prob, abstol=abstol,reltol=reltol, saveat=saveat, de$MethodOfSteps(de$Tsit5())))
system.time(sol_dct <- de$solve(prob_dct, abstol=abstol,reltol=reltol, saveat=saveat, de$MethodOfSteps(de$Tsit5())))
udf <- as.data.frame(cbind(sol$t,t(sapply(sol$u,identity)),t(sapply(sol_dct$u,identity))))
colnames(udf) <- c("t","S","E","P","Sdct","Edct","Pdct")

library(plotly)
fig <- plotly::plot_ly(udf, x = sol$t, y = ~S, name="S",type = 'scatter', mode = 'lines') %>% plotly::add_trace(y = ~E, name="E") %>% plotly::add_trace(y = ~P,name="P") %>%
  layout(
    title = list(text="Prediction for substrate volume 0.9 ml",y=0.95),
    xaxis = list(title = "time, secs"),
    yaxis = list(title = "volume, ml"),
    font = font
  )
orca(fig, file = paste(paste("Prediction_0dot9"),"png", sep = "."))


# Comparisons

library(plotly)

fig <- plotly::plot_ly(udf, x = sol$t, y = array(outer(array(outer(udf$P,unlist(parest[8]))),array(outer(array(outer(udf$P,3/2,FUN='^')),unlist(parest[9]))),FUN='-')), type = 'scatter', name="&#954;<sub>pred</sub>", mode = 'lines') %>% 
  plotly::add_trace(df, x = df$time_sec, y = df$substrate_0dot9, name="<sub>exp</sub>", mode='markers',marker = list(size = 3)) %>% 
  layout(
    title = list(text="Conductivity",y=0.95),
    xaxis = list(title = "time, secs"),
    yaxis = list(title = "muS/cm"),
    font=font
  ) %>%
  config(mathjax = "cdn")
orca(fig, file = paste(paste("P_pred_exp_0dot9"),"png", sep = "."))

fig <- plotly::plot_ly(udf, x = sol$t, y = array(outer(array(outer(udf$P,unlist(parest[8]))),array(outer(array(outer(udf$P,3/2,FUN='^')),unlist(parest[9]))),FUN='-')), type = 'scatter', name="&#954;<sub>pred,distributed</sub>", mode = 'lines') %>% 
  plotly::add_trace(udf, x = sol$t, y = array(outer(array(outer(udf$Pdct,unlist(parest[8]))),array(outer(array(outer(udf$Pdct,3/2,FUN='^')),unlist(parest[9]))),FUN='-')), type = 'scatter', name="&#954;<sub>pred,discrete</sub>", mode = 'lines') %>% 
  plotly::add_trace(df, x = df$time_sec, y = df$substrate_0dot9, name="&#954;<sub>exp</sub>", mode='markers',marker = list(size = 3)) %>% 
  layout(
    title = list(text="Conductivity for 0.9 ml of substrate",y=0.95,yref="container"),
    xaxis = list(title = "time, secs"),
    yaxis = list(title = "&mu; S/cm"),
    font=font
  ) %>%
  config(mathjax = "cdn")
orca(fig, file = paste(paste("P_pred_exp_0dot9_discrete_distributed"),"png", sep = "."))

##########################################
# S_0=0.1 ################################
##########################################
u0 <- unlist(parest[5:7]) 
u0[1] <- 0.1
#tspan <- c(max(0,min(df$time_sec)),max(0,max(df$time_sec)))
p <- parest
p[5] <- 0.1
constant_lags <- c(10.0)
JuliaCall::julia_assign("u0", u0)
JuliaCall::julia_assign("tspan", tspan)
JuliaCall::julia_assign("constant_lags", constant_lags)
JuliaCall::julia_assign("p", p)
prob <- JuliaCall::julia_eval("DDEProblem(f, u0, h, tspan, constant_lags = constant_lags, p)")
prob_dct <- JuliaCall::julia_eval("DDEProblem(f_dct, u0, h, tspan, constant_lags = constant_lags, p)")
system.time(sol <- de$solve(prob, abstol=abstol,reltol=reltol, saveat=saveat, de$MethodOfSteps(de$Tsit5())))
system.time(sol_dct <- de$solve(prob_dct, abstol=abstol,reltol=reltol, saveat=saveat, de$MethodOfSteps(de$Tsit5())))
udf <- as.data.frame(cbind(sol$t,t(sapply(sol$u,identity)),t(sapply(sol_dct$u,identity))))
colnames(udf) <- c("t","S","E","P","Sdct","Edct","Pdct")


library(plotly)
fig <- plotly::plot_ly(udf, x = sol$t, y = ~S, name="S",type = 'scatter', mode = 'lines') %>% plotly::add_trace(y = ~E, name="E") %>% plotly::add_trace(y = ~P,name="P") %>%
  layout(
    title = list(text="Prediction for substrate volume 0.1 ml",y=0.95),
    xaxis = list(title = "time, secs"),
    yaxis = list(title = "volume, ml"),
    font = font
  )
orca(fig, file = paste(paste("Prediction_0dot1"),"png", sep = "."))


# Comparisons

library(plotly)

fig <- plotly::plot_ly(udf, x = sol$t, y = array(outer(array(outer(udf$P,unlist(parest[8]))),array(outer(array(outer(udf$P,3/2,FUN='^')),unlist(parest[9]))),FUN='-')), type = 'scatter', name="&#954;<sub>pred</sub>", mode = 'lines') %>% 
  plotly::add_trace(df, x = df$time_sec, y = df$substrate_0dot1, name="&#954;<sub>exp</sub>", mode='markers',marker = list(size = 3)) %>% 
  layout(
    title = list(text="Conductivity",y=0.95),
    xaxis = list(title = "time, secs"),
    yaxis = list(title = "&mu; S/cm"),
    font=font
  ) %>%
  config(mathjax = "cdn")
orca(fig, file = paste(paste("P_pred_exp_0dot1"),"png", sep = "."))

fig <- plotly::plot_ly(udf, x = sol$t, y = array(outer(array(outer(udf$P,unlist(parest[8]))),array(outer(array(outer(udf$P,3/2,FUN='^')),unlist(parest[9]))),FUN='-')), type = 'scatter', name="&#954;<sub>pred,distributed</sub>", mode = 'lines') %>% 
  plotly::add_trace(udf, x = sol$t, y = array(outer(array(outer(udf$Pdct,unlist(parest[8]))),array(outer(array(outer(udf$Pdct,3/2,FUN='^')),unlist(parest[9]))),FUN='-')), type = 'scatter', name="&#954;<sub>pred,discrete</sub>", mode = 'lines') %>% 
  plotly::add_trace(df, x = df$time_sec, y = df$substrate_0dot1, name="&#954;<sub>exp</sub>", mode='markers',marker = list(size = 3)) %>% 
  layout(
    title = list(text="Conductivity for 0.1 ml of substrate",y=0.95),
    xaxis = list(title = "time, secs"),
    yaxis = list(title = "&mu; S/cm"),
    font=font
  ) %>%
  config(mathjax = "cdn")
orca(fig, file = paste(paste("P_pred_exp_0dot1_discrete_distributed"),"png", sep = "."))


##########################################
# S_0=0.3 ################################
##########################################
u0 <- unlist(parest[5:7]) 
u0[1] <- 0.3
#tspan <- c(max(0,min(df$time_sec)),max(0,max(df$time_sec)))
p <- parest
p[5] <- 0.3
constant_lags <- c(10.0)
JuliaCall::julia_assign("u0", u0)
JuliaCall::julia_assign("tspan", tspan)
JuliaCall::julia_assign("constant_lags", constant_lags)
JuliaCall::julia_assign("p", p)
prob <- JuliaCall::julia_eval("DDEProblem(f, u0, h, tspan, constant_lags = constant_lags, p)")
prob_dct <- JuliaCall::julia_eval("DDEProblem(f_dct, u0, h, tspan, constant_lags = constant_lags, p)")
system.time(sol <- de$solve(prob, abstol=abstol,reltol=reltol, saveat=saveat, de$MethodOfSteps(de$Tsit5())))
system.time(sol_dct <- de$solve(prob_dct, abstol=abstol,reltol=reltol, saveat=saveat, de$MethodOfSteps(de$Tsit5())))
udf <- as.data.frame(cbind(sol$t,t(sapply(sol$u,identity)),t(sapply(sol_dct$u,identity))))
colnames(udf) <- c("t","S","E","P","Sdct","Edct","Pdct")


library(plotly)
fig <- plotly::plot_ly(udf, x = sol$t, y = ~S, name="S",type = 'scatter', mode = 'lines') %>% plotly::add_trace(y = ~E, name="E") %>% plotly::add_trace(y = ~P,name="P") %>%
  layout(
    title = list(text="Prediction for substrate volume 0.3 ml",y=0.95),
    xaxis = list(title = "time, secs"),
    yaxis = list(title = "volume, ml"),
    font = font
  )
orca(fig, file = paste(paste("Prediction_0dot3"),"png", sep = "."))


# Comparisons

library(plotly)

fig <- plotly::plot_ly(udf, x = sol$t, y = array(outer(array(outer(udf$P,unlist(parest[8]))),array(outer(array(outer(udf$P,3/2,FUN='^')),unlist(parest[9]))),FUN='-')), type = 'scatter', name="&#954;<sub>pred</sub>", mode = 'lines') %>% 
  plotly::add_trace(df, x = df$time_sec, y = df$substrate_0dot3, name="P<sub>exp</sub>", mode='markers',marker = list(size = 3)) %>% 
  layout(
    title = list(text="Conductivity",y=0.95),
    xaxis = list(title = "time, secs"),
    yaxis = list(title = "&mu; S/cm"),
    font=font
  ) %>%
  config(mathjax = "cdn")
orca(fig, file = paste(paste("P_pred_exp_0dot3"),"png", sep = "."))

fig <- plotly::plot_ly(udf, x = sol$t, y = array(outer(array(outer(udf$P,unlist(parest[8]))),array(outer(array(outer(udf$P,3/2,FUN='^')),unlist(parest[9]))),FUN='-')), type = 'scatter', name="&#954;<sub>pred,distributed</sub>", mode = 'lines') %>% 
  plotly::add_trace(udf, x = sol$t, y = array(outer(array(outer(udf$Pdct,unlist(parest[8]))),array(outer(array(outer(udf$Pdct,3/2,FUN='^')),unlist(parest[9]))),FUN='-')), type = 'scatter', name="&#954;<sub>pred,discrete</sub>", mode = 'lines') %>% 
  plotly::add_trace(df, x = df$time_sec, y = df$substrate_0dot3, name="&#954;<sub>exp</sub>", mode='markers',marker = list(size = 3)) %>% 
  layout(
    title = list(text="Conductivity  for 0.3 ml of substrate",y=0.95),
    xaxis = list(title = "time, secs"),
    yaxis = list(title = "&mu; S/cm"),
    font=font
  ) %>%
  config(mathjax = "cdn")
orca(fig, file = paste(paste("P_pred_exp_0dot3_discrete_distributed"),"png", sep = "."))

##########################################
# S_0=1.5 ################################
##########################################
u0 <- unlist(parest[5:7]) 
u0[1] <- 1.5
#tspan <- c(max(0,min(df$time_sec)),max(0,max(df$time_sec)))
p <- parest
p[5] <- 1.5
constant_lags <- c(10.0)
JuliaCall::julia_assign("u0", u0)
JuliaCall::julia_assign("tspan", tspan)
JuliaCall::julia_assign("constant_lags", constant_lags)
JuliaCall::julia_assign("p", p)
prob <- JuliaCall::julia_eval("DDEProblem(f, u0, h, tspan, constant_lags = constant_lags, p)")
prob_dct <- JuliaCall::julia_eval("DDEProblem(f_dct, u0, h, tspan, constant_lags = constant_lags, p)")
system.time(sol <- de$solve(prob, abstol=abstol,reltol=reltol, saveat=saveat, de$MethodOfSteps(de$Tsit5())))
system.time(sol_dct <- de$solve(prob_dct, abstol=abstol,reltol=reltol, saveat=saveat, de$MethodOfSteps(de$Tsit5())))
udf <- as.data.frame(cbind(sol$t,t(sapply(sol$u,identity)),t(sapply(sol_dct$u,identity))))
colnames(udf) <- c("t","S","E","P","Sdct","Edct","Pdct")

library(plotly)
fig <- plotly::plot_ly(udf, x = sol$t, y = ~S, name="S",type = 'scatter', mode = 'lines') %>% plotly::add_trace(y = ~E, name="E") %>% plotly::add_trace(y = ~P,name="P") %>%
  layout(
    title = list(text="Prediction for substrate volume 1.5 ml",y=0.95),
    xaxis = list(title = "time, secs"),
    yaxis = list(title = "volume, ml"),
    font = font
  )
orca(fig, file = paste(paste("Prediction_1dot5"),"png", sep = "."))


# Comparisons

library(plotly)

fig <- plotly::plot_ly(udf, x = sol$t, y = array(outer(array(outer(udf$P,unlist(parest[8]))),array(outer(array(outer(udf$P,3/2,FUN='^')),unlist(parest[9]))),FUN='-')), type = 'scatter', name="&#954;<sub>pred</sub>", mode = 'lines') %>% 
  plotly::add_trace(df, x = df$time_sec, y = df$substrate_1dot5, name="&#954;<sub>exp</sub>", mode='markers',marker = list(size = 4)) %>% 
  layout(
    title = list(text="Conductivity",y=0.95),
    xaxis = list(title = "time, secs"),
    yaxis = list(title = "&mu; S/cm"),
    font=font
  ) %>%
  config(mathjax = "cdn")
orca(fig, file = paste(paste("P_pred_exp_1dot5"),"png", sep = "."))

fig <- plotly::plot_ly(udf, x = sol$t, y = array(outer(array(outer(udf$P,unlist(parest[8]))),array(outer(array(outer(udf$P,3/2,FUN='^')),unlist(parest[9]))),FUN='-')), type = 'scatter', name="&#954;<sub>pred,distributed</sub>", mode = 'lines') %>% 
  plotly::add_trace(udf, x = sol$t, y = array(outer(array(outer(udf$Pdct,unlist(parest[8]))),array(outer(array(outer(udf$Pdct,3/2,FUN='^')),unlist(parest[9]))),FUN='-')), type = 'scatter', name="&#954;<sub>pred,discrete</sub>", mode = 'lines') %>% 
  plotly::add_trace(df, x = df$time_sec, y = df$substrate_1dot5, name="&#954;<sub>exp</sub>", mode='markers',marker = list(size = 3)) %>% 
  layout(
    title = list(text="Conductivity for 1.5 ml of substrate",y=0.95),
    xaxis = list(title = "time, secs"),
    yaxis = list(title = "&mu; S/cm"),
    font=font
  ) %>%
  config(mathjax = "cdn")
orca(fig, file = paste(paste("P_pred_exp_1dot5_discrete_distributed"),"png", sep = "."))


########################################################
# delay distribution functions #########################
########################################################
a=as.numeric(parest[2])
m=as.numeric(parest[3])
tau_m=as.numeric(parest[4])
rmse=sqrt((m+1)/a^2)
psi_func = function(x,a,m,tau_m) {(a^(m+1)/gamma(m+1))*(x-tau_m)^m * exp(-a*(x-tau_m))}
arguments = seq(tau_m,tau_m+(m+1)/a+3*rmse,0.05)

data<-data.frame(time=c(arguments),psi=psi_func(c(arguments),a,m,tau_m))


png(paste(paste("delay_densities_0dot3"),"png", sep = "."), width = 800, height = 600)
densities <- ggplot()+geom_line(data=data,aes(x=time,y=psi),size=1.5)+
  xlab("t, secs") + ylab(expression(italic(psi))) + 
  my_theme
plot(densities)
dev.off()

save(list=c("parest","parms","lower","upper"),file = paste("S_0dot3",".Rdata",sep = "."))

save.image(file="BrownModelGammaDist_multiple_volumes.RData")