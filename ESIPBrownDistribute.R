# ATTENTION!!! R version 3.6.1 + Julia 1.5.3
# gamma distributed delays model of electrochemical biosensor: Enzyme+Substrate+Inhibitor
# ver 4 of the script - error analysis was added
# ver 5 k_w coefficient was added for the product of reaction
# video on Michaelis-Menten equilibrium
# https://www.khanacademy.org/test-prep/mcat/biomolecules/enzyme-kinetics/e/enzyme-kinetics-questions


# https://www.r-bloggers.com/learning-r-parameter-fitting-for-models-involving-differential-equations/

require(xlsx)
#install.packages('xlsx', type = 'source', INSTALL_opts='--no-multiarch')
require(dplyr)
library(ggplot2)
library(reshape2) # library for reshaping data (tall-narrow <-> short-wide)
library(deSolve) # library for solving differential equations
library(minpack.lm) # library for least squares fit using levenberg-marquart algorithm

dname<-"C:/My_doc/KlosWitkowska/ElectroChemicalBiosensorBSA/Plots/gammadistribution_multiple_volumes"

dname1<-"ESIP_ver5"

dir.create(file.path(dname,dname1), showWarnings = FALSE)
setwd(file.path(dname,dname1))

my_theme <- theme_bw() +  theme(axis.text=element_text(size=14),axis.title=element_text(size=16),legend.title=element_text(size=16),legend.text = element_text(size=16)) 


# ingest experimental data of conductivity

df <- read.csv("C:/My_doc/Sverstiuk/Dzyadevych/Fw__1_mM_Bu_+_1,3,5_mkM_GA/data1_3_5_10.csv", header = TRUE, sep="\t")
colnames(df)[1] <- "time_sec"

df1 <- df[which(df$time_sec>=601&df$time_sec<=1332),]
inhibitor_1 <- df1$chan_2
df1 <- df[which(df$time_sec>=2954.5&df$time_sec<=3680.5),]
inhibitor_3 <- df1$chan_2
df1 <- df[which(df$time_sec>=6656.0&df$time_sec<=7210.5),]
inhibitor_5 <- df1$chan_2
df1 <- df[which(df$time_sec>=9619.0&df$time_sec<=9977.0),]
inhibitor_10 <- df1$chan_2

min_dim <- min(length(inhibitor_1),length(inhibitor_3),length(inhibitor_5),length(inhibitor_10))

df <- data.frame(df$time_sec[1:min_dim],inhibitor_1[1:min_dim],inhibitor_3[1:min_dim],inhibitor_5[1:min_dim],inhibitor_10[1:min_dim])
colnames(df) <- c("time_sec","inhibitor_1","inhibitor_3","inhibitor_5","inhibitor_10")
tmpval <- inhibitor_3[1]-inhibitor_1[1]
df$inhibitor_3 <- df$inhibitor_3 - c(rep(tmpval,min_dim)) 
tmpval <- inhibitor_5[1]-inhibitor_1[1]
df$inhibitor_5 <- df$inhibitor_5 - c(rep(tmpval,min_dim)) 
tmpval <- inhibitor_10[1]-inhibitor_1[1]
df$inhibitor_10 <- df$inhibitor_10 - c(rep(tmpval,min_dim)) 

# plot data
tmp=melt(df,id.vars=c("time_sec"),variable.name="conductivity",value.name="muS/cm")

png(paste(paste("inhibitor_1mM"),"png", sep = "."), width = 800, height = 600)
p<-ggplot(data=tmp[which(tmp$conductivity=="inhibitor_1"),],aes(x=time_sec,y=`muS/cm`,color=conductivity))+geom_point(size=1)+my_theme
plot(p)
dev.off()

png(paste(paste("inhibitor_3mM"),"png", sep = "."), width = 800, height = 600)
p<-ggplot(data=tmp[which(tmp$conductivity=="inhibitor_3"),],aes(x=time_sec,y=`muS/cm`,color=conductivity))+geom_point(size=1)+my_theme
plot(p)
dev.off()

png(paste(paste("inhibitor_5mM"),"png", sep = "."), width = 800, height = 600)
p<-ggplot(data=tmp[which(tmp$conductivity=="inhibitor_5"),],aes(x=time_sec,y=`muS/cm`,color=conductivity))+geom_point(size=1)+my_theme
plot(p)
dev.off()

png(paste(paste("inhibitor_10mM"),"png", sep = "."), width = 800, height = 600)
p<-ggplot(data=tmp[which(tmp$conductivity=="inhibitor_10"),],aes(x=time_sec,y=`muS/cm`,color=conductivity))+geom_point(size=1)+my_theme
plot(p)
dev.off()

###################
# Julia examples 
#  https://github.com/SciML/diffeqr
####################

# install.packages("diffeqr")
# devtools::install_github('SciML/diffeqr', build_vignettes=T)

# diffeqr::diffeq_setup()  # installation of Julia

library(JuliaCall)
library(diffeqr)

de <- diffeqr::diffeq_setup()

f_dct <- JuliaCall::julia_eval("function f_dct(du, u, h, p, t)
  alpha_01d=p[1]
  alpha_11d=p[2]
  beta_10d=p[3]
  beta_11d=p[4]
  a_01=p[5]
  m_01=p[6]
  tau_m01=p[7]
  a_11=p[8]
  m_11=p[9]
  tau_m11=p[10]
  b_10=p[11]
  k_10=p[12]
  h_m10=p[13]
  b_11=p[14]
  k_11=p[15]
  h_m11=p[16]
  n_P_0=p[17]
  n_S1_0=p[18]
  n_E_0=p[19]
  n_I1_0=p[20]
  k_w=p[21]
  tau_01=tau_m01+(m_01+1)/a_01
  tau_11=tau_m11+(m_11+1)/a_11
  h_10=h_m10+(k_10+1)/b_10
  h_11=h_m11+(k_11+1)/b_11
  du[1]=alpha_01d*h(p,t-tau_01)[3]*h(p,t-tau_01)[2]+alpha_11d*h(p,t-tau_11)[3]*h(p,t-tau_11)[2]-k_w*u[1]
  du[2]=-(alpha_01d+alpha_11d)*u[2]*u[3]
  du[3]=alpha_01d*h(p,t-tau_01)[3]*h(p,t-tau_01)[2]+alpha_11d*h(p,t-tau_11)[3]*h(p,t-tau_11)[2]-(alpha_01d+alpha_11d)*u[2]*u[3]+beta_10d*h(p,t-h_10)[4]*h(p,t-h_10)[3]+beta_11d*h(p,t-h_11)[4]*h(p,t-h_11)[3]-(beta_10d+beta_11d)*u[4]*u[3]
  du[4]=-(beta_10d+beta_11d)*u[4]*u[3]
end")  
h <- JuliaCall::julia_eval("function h(p, t)
  [p[17],p[18],p[19],p[20]]
end")
psi <- JuliaCall::julia_eval("function psi(a, m, tau_m, tau)
  if tau <= tau_m
    res = 0
  else 
    res = (a^(m+1)/gamma(m+1))*(tau-tau_m)^m * exp(-a*(tau-tau_m))
  end
  res
end")

### for distributed delay############################

f <- JuliaCall::julia_eval("function f(du, u, h, p, t)
  alpha_01d=p[1]
  alpha_11d=p[2]
  beta_10d=p[3]
  beta_11d=p[4]
  a_01=p[5]
  m_01=p[6]
  tau_m01=p[7]
  a_11=p[8]
  m_11=p[9]
  tau_m11=p[10]
  b_10=p[11]
  k_10=p[12]
  h_m10=p[13]
  b_11=p[14]
  k_11=p[15]
  h_m11=p[16]
  n_P_0=p[17]
  n_S1_0=p[18]
  n_E_0=p[19]
  n_I1_0=p[20]
  k_w=p[21]
  n=100
  lags_alpha_01=range(1/n, length=n, stop=(tau_m01+(m_01+1)/a_01+3*rmse(a_01,m_01))*(1-1/n))
  lags_alpha_11=range(1/n, length=n, stop=(tau_m11+(m_11+1)/a_11+3*rmse(a_11,m_11))*(1-1/n))
  lags_beta_10=range(1/n, length=n, stop=(h_m10+(k_10+1)/b_10+3*rmse(b_10,k_10))*(1-1/n)) 
  lags_beta_11=range(1/n, length=n, stop=(h_m11+(k_11+1)/b_11+3*rmse(b_11,k_11))*(1-1/n))
  #lags_alpha_01=lags(n,a_01,m_01,tau_m01)
  #lags_alpha_11=lags(n,a_11,m_11,tau_m11)
  #lags_beta_10=lags(n,b_10,k_10,h_m10)
  #lags_beta_11=lags(n,b_11,k_11,h_m11)
  du[1]=alpha_01d*tau_M(a_01,m_01,tau_m01)*(1/n)*sum([psi(a_01,m_01,tau_m01,s)*h(p,t-s)[3]*h(p,t-s)[2] for s in lags_alpha_01])+alpha_11d*tau_M(a_11,m_11,tau_m11)*(1/n)*sum([psi(a_11,m_11,tau_m11,s)*h(p,t-s)[3]*h(p,t-s)[2] for s in lags_alpha_11])-k_w*u[1]
  du[2]=-(alpha_01d+alpha_11d)*u[2]*u[3]
  du[3]=alpha_01d*tau_M(a_01,m_01,tau_m01)*(1/n)*sum([psi(a_01,m_01,tau_m01,s)*h(p,t-s)[3]*h(p,t-s)[2] for s in lags_alpha_01])+alpha_11d*tau_M(a_11,m_11,tau_m11)*(1/n)*sum([psi(a_11,m_11,tau_m11,s)*h(p,t-s)[3]*h(p,t-s)[2] for s in lags_alpha_11])-(alpha_01d+alpha_11d)*u[2]*u[3]+beta_10d*tau_M(b_10,k_10,h_m10)*(1/n)*sum([psi(b_10,k_10,h_m10,s)*h(p,t-s)[3]*h(p,t-s)[4] for s in lags_beta_10])+beta_11d*tau_M(b_11,k_11,h_m11)*(1/n)*sum([psi(b_11,k_11,h_m11,s)*h(p,t-s)[3]*h(p,t-s)[4] for s in lags_beta_11])-(beta_10d+beta_11d)*u[4]*u[3]
  du[4]=-(beta_10d+beta_11d)*u[4]*u[3]
end")  
h <- JuliaCall::julia_eval("function h(p, t)
  [p[17],p[18],p[19],p[20]]
end")
rmse <- JuliaCall::julia_eval("function rmse(a, m)
  sqrt((m+1)/a^2)
end")
lags <- JuliaCall::julia_eval("function lags(n,a,m,tau_m)
  return range(1/n, length=n, stop=(tau_m+(m+1)/a+3*rmse(a,m))*(1-1/n))
end")
tau_M <- JuliaCall::julia_eval("function tau_M(a,m,tau_m)
  tau_m+(m+1)/a+3*rmse(a,m)
end")
psi <- JuliaCall::julia_eval("function psi(a, m, tau_m, tau)
  if tau <= tau_m
    res = 0
  else 
    res = (a^(m+1)/gamma(m+1))*(tau-tau_m)^m * exp(-a*(tau-tau_m))
  end
  res
end")

#####################################################

# I_0 = 1

u0 <- c(30,0.1,1,1)
tspan <- c(max(0,min(df$time_sec)),max(0,max(df$time_sec)))
abstol <- 1e-2
reltol <- 1e-2
saveat <- round(seq(from=min(df$time),to=max(df$time),length.out=500))
saveat <- unique(saveat)
p <- c(0.05,0.05,0.05,0.05,1,20,5,1,20,5,1,20,5,1,20,5,u0[1:4],0.168)
constant_lags <- c(100.0)
JuliaCall::julia_assign("u0", u0)
JuliaCall::julia_assign("tspan", tspan)
JuliaCall::julia_assign("constant_lags", constant_lags)
JuliaCall::julia_assign("p", p)
prob_dct <- JuliaCall::julia_eval("using DifferentialEquations; using SpecialFunctions; DDEProblem(f_dct, u0, h, tspan, p)")
prob <- JuliaCall::julia_eval("using DifferentialEquations; using SpecialFunctions; DDEProblem(f, u0, h, tspan, p)")
system.time(sol <- de$solve(prob,abstol=abstol,reltol=reltol, saveat=saveat, de$MethodOfSteps(de$Tsit5())))
udf <- as.data.frame(cbind(sol$t,t(sapply(sol$u,identity))))
colnames(udf) <- c("t","P","S_1","E","I_1")

library(plotly)

font <- list(
  family = "arial",
  size = 14,
  color = 'black')

fig <- plotly::plot_ly(udf, x = sol$t, y = ~P, type = 'scatter', mode = 'lines') 
fig <- fig %>% layout(title=list(text="Example",y=0.95), font=font)
orca(fig, file = paste(paste("Pexample"),"png", sep = "."))
fig <- plotly::plot_ly(udf, x = sol$t, y = ~S_1, type = 'scatter', mode = 'lines') 
fig <- fig %>% layout(title=list(text="Example",y=0.95), font=font)
orca(fig, file = paste(paste("S_1example"),"png", sep = "."))
fig <- plotly::plot_ly(udf, x = sol$t, y = ~E, type = 'scatter', mode = 'lines') 
fig <- fig %>% layout(title=list(text="Example",y=0.95), font=font)
orca(fig, file = paste(paste("Eexample"),"png", sep = "."))
fig <- plotly::plot_ly(udf, x = sol$t, y = ~I_1, type = 'scatter', mode = 'lines') 
fig <- fig %>% layout(title=list(text="Example",y=0.95), font=font)
orca(fig, file = paste(paste("I_1example"),"png", sep = "."))
fig <- plotly::plot_ly(udf, x = sol$t, y = ~S_1, type = 'scatter', mode = 'lines') %>% plotly::add_trace(y = ~E) %>% plotly::add_trace(y = ~I_1) %>% plotly::add_trace(y = ~P)
fig <- fig %>% layout(title=list(text="Example",y=0.95), font=font)
orca(fig, file = paste(paste("PSEIexample"),"png", sep = "."))


#### Parameter identification

ssq=function(parms){
  de <- diffeqr::diffeq_setup()
  print(parms)
  
  f <- JuliaCall::julia_eval("function f(du, u, h, p, t)
  alpha_01d=p[1]
  alpha_11d=p[2]
  beta_10d=p[3]
  beta_11d=p[4]
  a_01=p[5]
  m_01=p[6]
  tau_m01=p[7]
  a_11=p[8]
  m_11=p[9]
  tau_m11=p[10]
  b_10=p[11]
  k_10=p[12]
  h_m10=p[13]
  b_11=p[14]
  k_11=p[15]
  h_m11=p[16]
  n_P_0=p[17]
  n_S1_0=p[18]
  n_E_0=p[19]
  n_I1_0=p[20]
  k_w=p[21]
  n=100
  lags_alpha_01=range(1/n, length=n, stop=(tau_m01+(m_01+1)/a_01+3*rmse(a_01,m_01))*(1-1/n))
  lags_alpha_11=range(1/n, length=n, stop=(tau_m11+(m_11+1)/a_11+3*rmse(a_11,m_11))*(1-1/n))
  lags_beta_10=range(1/n, length=n, stop=(h_m10+(k_10+1)/b_10+3*rmse(b_10,k_10))*(1-1/n)) 
  lags_beta_11=range(1/n, length=n, stop=(h_m11+(k_11+1)/b_11+3*rmse(b_11,k_11))*(1-1/n))
  #lags_alpha_01=lags(n,a_01,m_01,tau_m01)
  #lags_alpha_11=lags(n,a_11,m_11,tau_m11)
  #lags_beta_10=lags(n,b_10,k_10,h_m10)
  #lags_beta_11=lags(n,b_11,k_11,h_m11)
  du[1]=alpha_01d*tau_M(a_01,m_01,tau_m01)*(1/n)*sum([psi(a_01,m_01,tau_m01,s)*h(p,t-s)[3]*h(p,t-s)[2] for s in lags_alpha_01])+alpha_11d*tau_M(a_11,m_11,tau_m11)*(1/n)*sum([psi(a_11,m_11,tau_m11,s)*h(p,t-s)[3]*h(p,t-s)[2] for s in lags_alpha_11])-k_w*u[1]
  du[2]=-(alpha_01d+alpha_11d)*u[2]*u[3]
  du[3]=alpha_01d*tau_M(a_01,m_01,tau_m01)*(1/n)*sum([psi(a_01,m_01,tau_m01,s)*h(p,t-s)[3]*h(p,t-s)[2] for s in lags_alpha_01])+alpha_11d*tau_M(a_11,m_11,tau_m11)*(1/n)*sum([psi(a_11,m_11,tau_m11,s)*h(p,t-s)[3]*h(p,t-s)[2] for s in lags_alpha_11])-(alpha_01d+alpha_11d)*u[2]*u[3]+beta_10d*tau_M(b_10,k_10,h_m10)*(1/n)*sum([psi(b_10,k_10,h_m10,s)*h(p,t-s)[3]*h(p,t-s)[4] for s in lags_beta_10])+beta_11d*tau_M(b_11,k_11,h_m11)*(1/n)*sum([psi(b_11,k_11,h_m11,s)*h(p,t-s)[3]*h(p,t-s)[4] for s in lags_beta_11])-(beta_10d+beta_11d)*u[4]*u[3]
  du[4]=-(beta_10d+beta_11d)*u[4]*u[3]
end")  
  h <- JuliaCall::julia_eval("function h(p, t)
  [p[17],p[18],p[19],p[20]]
end")
  rmse <- JuliaCall::julia_eval("function rmse(a, m)
  sqrt((m+1)/a^2)
end")
  lags <- JuliaCall::julia_eval("function lags(n,a,m,tau_m)
  return range(1/n, length=n, stop=(tau_m+(m+1)/a+3*rmse(a,m))*(1-1/n))
end")
  tau_M <- JuliaCall::julia_eval("function tau_M(a,m,tau_m)
  tau_m+(m+1)/a+3*rmse(a,m)
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
  ################ for 1 #################################
  u0 <- c(df$inhibitor_1[1],0.1,1,1)
  tspan <- c(max(0,min(df$time_sec)),max(0,max(df$time_sec)))
  abstol <- 1e-2
  reltol <- 1e-2
  saveat <- round(seq(from=min(df$time),to=max(df$time),length.out=500))
  saveat <- unique(saveat)
  parms[20]<-1
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
  colnames(udf) <- c("t","P","S_1","E","I_1")
  
  # Filter data that contains time points where data is available
  
  outdf=udf[udf$t %in% df$time_sec,]
  
  # Evaluate predicted vs experimental residual
  preddf=outdf$P[which(outdf$t %in% df$time_sec)]
  expdf=df$inhibitor_1[which(outdf$t %in% df$time_sec)]
  ssqres_1=unlist(parms)[22]*preddf - unlist(parms)[23]*'^'(preddf,3/2) - expdf
  
  ################ for 3 #################################
  u0 <- c(df$inhibitor_3[1],0.1,1,3)
  tspan <- c(max(0,min(df$time_sec)),max(0,max(df$time_sec)))
  abstol <- 1e-2
  reltol <- 1e-2
  saveat <- round(seq(from=min(df$time),to=max(df$time),length.out=500))
  saveat <- unique(saveat)
  parms[20]<-3
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
  colnames(udf) <- c("t","P","S_1","E","I_1")
  
  # Filter data that contains time points where data is available
  
  outdf=udf[udf$t %in% df$time_sec,]
  
  # Evaluate predicted vs experimental residual
  preddf=outdf$P[which(outdf$t %in% df$time_sec)]
  expdf=df$inhibitor_3[which(outdf$t %in% df$time_sec)]
  ssqres_3=unlist(parms)[22]*preddf - unlist(parms)[23]*'^'(preddf,3/2) - expdf
  ################ for 5 #################################
  u0 <- c(df$inhibitor_5[1],0.1,1,5)
  tspan <- c(max(0,min(df$time_sec)),max(0,max(df$time_sec)))
  abstol <- 1e-2
  reltol <- 1e-2
  saveat <- round(seq(from=min(df$time),to=max(df$time),length.out=500))
  saveat <- unique(saveat)
  parms[20]<-5
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
  colnames(udf) <- c("t","P","S_1","E","I_1")
  
  # Filter data that contains time points where data is available
  
  outdf=udf[udf$t %in% df$time_sec,]
  
  # Evaluate predicted vs experimental residual
  preddf=outdf$P[which(outdf$t %in% df$time_sec)]
  expdf=df$inhibitor_5[which(outdf$t %in% df$time_sec)]
  ssqres_5=unlist(parms)[22]*preddf - unlist(parms)[23]*'^'(preddf,3/2) - expdf
  
  ################ for 10 #################################
  u0 <- c(df$inhibitor_10[1],0.1,1,10)
  tspan <- c(max(0,min(df$time_sec)),max(0,max(df$time_sec)))
  abstol <- 1e-2
  reltol <- 1e-2
  saveat <- round(seq(from=min(df$time),to=max(df$time),length.out=500))
  saveat <- unique(saveat)
  parms[20]<-10
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
  colnames(udf) <- c("t","P","S_1","E","I_1")
  
  # Filter data that contains time points where data is available
  
  outdf=udf[udf$t %in% df$time_sec,]
  
  # Evaluate predicted vs experimental residual
  preddf=outdf$P[which(outdf$t %in% df$time_sec)]
  expdf=df$inhibitor_10[which(outdf$t %in% df$time_sec)]
  ssqres_10=unlist(parms)[22]*preddf - unlist(parms)[23]*'^'(preddf,3/2) - expdf
  
  ##########################################
  # return predicted vs experimental residual
  ssqres <- c(ssqres_1,ssqres_3,ssqres_5,ssqres_10)
  #print(ssqres)
  
  print(mean(abs(ssqres)))
  
  plot(x=udf$t,y=udf$P)
  
  return(ssqres)
}

# parameter fitting 
# initial guess for parameters

#parms <- c( 1e-3,0.128817,1e-3,0.156118,1,20,5,1,20,5,1,20,5,1,20,5,3,1e-10,1.9,1,0.168,4.31,8.085)
#parms <- c( 1e-3,0.128817,1e-3,0.156118,1,20,5,1,20,5,1,20,5,1,20,5,3,1e-10,1.9,1,0.168,4.31,8.085)
#parms <- c( 1e-10,1.28817,1e-10,95611.8,1e-10,1e-10,412.122,23.2164,933.821,1e-10,1,20,5,1,20,5,3,1e-10,1834.9,1,43759.1,8015.85)
#parms <- c(0.05,0.05,0.05,0.05,1,20,5,1,20,5,1,20,5,1,20,5,u0[1:4],729,246)
##parms <- c(1e-10,0.185906,1e-10,105376,6.91667,253.517,1e-10,1e-10,1e-10,114.871,1,20,5,1,20,5,3,219.578,1e-10,1,729,246)
lower <- c(rep(1.e-4,23))
upper <- c(rep(1.e4,23))  
parms <- c( rep(1.e-2,4), rep(c(1,7,5),4),3,1e-10,1.9,1,0.168,4.31,8.085)
#parms <- c( rep(5.e-2,4), rep(c(1,20,5),4),30,1e-10,1.9,1,0.168,4.31,8.085)

######################################################
################ fitting using levenberg marquart algorithm
if(TRUE){
  fitval=nls.lm(par=parms,lower=lower,upper=upper,fn=ssq, control = nls.lm.control(maxiter = 10,nprint = 2))
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

########## # parameter fitting using NLopt
if(FALSE){
  # Steven G. Johnson, The NLopt nonlinear-optimization package, http://ab-initio.mit.edu/nlopt
  
  library(nloptr)
  fn.sir <- function(x) {
    sqrt(sum(ssq(x)^2))
  }
  
  # equalities
  eval_g_eq <- function( x ) {
    constr <- c( x[14]-x[15]-x[16]-x[17]-x[18]-x[19]-x[20]-x[21] )
    grad <- c( 0,
               0,
               0,
               0,
               0,
               0,
               0,
               0,
               0,
               0,
               0,
               0,
               0,
               1,
               -1,
               -1,
               -1,
               -1,
               -1,
               -1,
               -1,
               0,
               0,
               0,
               0,
               0,
               0,
               0,
               0)
    return( list( "constraints"=constr, "jacobian"=grad ) )
  }
  
  
  local_opts <- list( "algorithm" = "NLOPT_LD_MMA",
                      "xtol_rel" = 1.0e-7 )
  opts <- list( "algorithm" = "NLOPT_LD_MMA",
                "xtol_rel" = 1.0e-7,
                "maxeval" = 1000,
                "local_opts" = local_opts )
  res <- nloptr( x0=parms,
                 eval_f=fn.sir,
                 lb=lower,
                 ub=upper,
                 eval_g_eq=eval_g_eq,
                 opts=opts)
  print( res )
  
}
########## # parameter fitting using COBYLA
# M. J. D. Powell, ``A direct search optimization method that models the objective and constraint functions by linear interpolation,'' in Advances in Optimization and Numerical Analysis, eds. S. Gomez and J.-P. Hennart (Kluwer Academic: Dordrecht, 1994), p. 51-67.
if(FALSE){
  library(nloptr)
  fn.sep <- function(x) {
    sqrt(sum(ssq(x)^2))
  }
  hin.sep <- function(x) {
    h <- numeric(18)
    h[1:9] <- upper - x
    h[10:18] <- x - lower
    #h[60] <- x[14]-x[15]-x[16]-x[17]-x[18]-x[19]-x[20]-x[21]
    #h[61] <- -x[14]+x[15]+x[16]+x[17]+x[18]+x[19]+x[20]+x[21]
    #h[2] <- -x[14]+x[15]+x[16]+x[17]+x[18]+x[19]+x[20]+x[21]
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
# I_1_0=1 ################################
##########################################
u0 <- unlist(parest[17:20]) 
u0[4] <- 1
#tspan <- c(max(0,min(df$time_sec)),max(0,max(df$time_sec)))
p <- parest
constant_lags <- c(100.0)
JuliaCall::julia_assign("u0", u0)
JuliaCall::julia_assign("tspan", tspan)
JuliaCall::julia_assign("constant_lags", constant_lags)
JuliaCall::julia_assign("p", p)
prob <- JuliaCall::julia_eval("DDEProblem(f, u0, h, tspan, constant_lags = constant_lags, p)")
prob_dct <- JuliaCall::julia_eval("DDEProblem(f_dct, u0, h, tspan, constant_lags = constant_lags, p)")
system.time(sol <- de$solve(prob, abstol=abstol,reltol=reltol, saveat=saveat, de$MethodOfSteps(de$Tsit5())))
system.time(sol_dct <- de$solve(prob_dct, abstol=abstol,reltol=reltol, saveat=saveat, de$MethodOfSteps(de$Tsit5())))
udf <- as.data.frame(cbind(sol$t,t(sapply(sol$u,identity)),t(sapply(sol_dct$u,identity))))
colnames(udf) <- c("t","P","S_1","E","I_1","Pdct","S_1dct","Edct","I_1dct")

library(plotly)
fig <- plotly::plot_ly(udf, x = sol$t, y = ~S_1, name="S<sub>1</sub>",type = 'scatter', mode = 'lines') %>% plotly::add_trace(y = ~E, name="E") %>% plotly::add_trace(y = ~P,name="P") %>% plotly::add_trace(y = ~I_1,name="I<sub>1</sub>") %>%
  layout(
    title = list(text="Prediction for inhibitor volume 1 mkM",y=0.95),
    xaxis = list(title = "time, secs"),
    yaxis = list(title = "volume, ml"),
    font = font,
    legend = list(x = 0.7, y = 0.9,font = list(size = 20))
  )
orca(fig, file = paste(paste("Prediction_1"),"png", sep = "."))


# Comparisons

library(plotly)

P<-udf$P[which(udf$t %in% df$time_sec)]
Pdct<-udf$Pdct[which(udf$t %in% df$time_sec)]
udf_final<-cbind(df,P)
udf_final<-cbind(udf_final,Pdct)
udf_final$K=unlist(parest)[22]*udf_final$P - unlist(parest)[23]*'^'(udf_final$P,3/2)
udf_final$Kdct=unlist(parest)[22]*udf_final$Pdct - unlist(parest)[23]*'^'(udf_final$Pdct,3/2)

fig <- plotly::plot_ly(udf_final, x = udf_final$time_sec, y = udf_final$K, type = 'scatter', name="&#954;<sub>pred</sub>", mode = 'lines') %>% 
  plotly::add_trace(udf_final, x = udf_final$time_sec, y = udf_final$inhibitor_1, name="&#954;<sub>exp</sub>", mode='markers',marker = list(size = 3)) %>% 
  layout(
    title = list(text="Conductivity",y=0.95),
    xaxis = list(title = "time, secs"),
    yaxis = list(title = "muS/cm"),
    font=font,
    legend = list(x = 0.7, y = 0.9,font = list(size = 20))
  ) %>%
  config(mathjax = "cdn")
orca(fig, file = paste(paste("P_pred_exp_1"),"png", sep = "."))

fig <- plotly::plot_ly(udf_final, x = udf_final$time_sec, y = udf_final$K, type = 'scatter', name="&#954;<sub>pred,distributed</sub>", mode = 'lines') %>% 
  plotly::add_trace(udf_final, x = udf_final$time_sec, y = udf_final$Kdct, type = 'scatter', name="&#954;<sub>pred,discrete</sub>", mode = 'lines') %>% 
  plotly::add_trace(udf_final, x = udf_final$time_sec, y = udf_final$inhibitor_1, name="&#954;<sub>exp</sub>", mode='markers',marker = list(size = 3)) %>% 
  layout(
    title = list(text="Conductivity for 1 mkM of inhibitor",y=0.95,yref="container"),
    xaxis = list(title = "time, secs"),
    yaxis = list(title = "&mu; S/cm"),
    font=font,
    legend = list(x = 0.7, y = 0.9,font = list(size = 20))
  ) %>%
  config(mathjax = "cdn")
orca(fig, file = paste(paste("P_pred_exp_1_discrete_distributed"),"png", sep = "."))


fig <- plotly::plot_ly(udf_final,x = udf_final$time_sec, y = abs(udf_final$K - udf_final$inhibitor_1), type = 'scatter', name="&#954;<sub>pred,distributed</sub>", mode = 'lines') %>% 
  plotly::add_trace(udf_final,x = udf_final$time_sec, y = abs(udf_final$Kdct - udf_final$inhibitor_1), type = 'scatter', name="&#954;<sub>pred,discrete</sub>", mode = 'lines') %>% 
  layout(
    title = list(text="Errors of models for 1 mkM of inhibitor",y=0.95,yref="container"),
    xaxis = list(title = "time, secs"),
    yaxis = list(title = "&mu; S/cm"),
    font=font,
    legend = list(x = 0.7, y = 0.9,font = list(size = 20))
  ) %>%
  config(mathjax = "cdn")
orca(fig, file = paste(paste("P_pred_exp_1_error_discrete_distributed"),"png", sep = "."))


##########################################
# I_1_0=3 ################################
##########################################
u0 <- unlist(parest[17:20]) 
u0[4] <- 3
#tspan <- c(max(0,min(df$time_sec)),max(0,max(df$time_sec)))
p <- parest
constant_lags <- c(100.0)
JuliaCall::julia_assign("u0", u0)
JuliaCall::julia_assign("tspan", tspan)
JuliaCall::julia_assign("constant_lags", constant_lags)
JuliaCall::julia_assign("p", p)
prob <- JuliaCall::julia_eval("DDEProblem(f, u0, h, tspan, constant_lags = constant_lags, p)")
prob_dct <- JuliaCall::julia_eval("DDEProblem(f_dct, u0, h, tspan, constant_lags = constant_lags, p)")
system.time(sol <- de$solve(prob, abstol=abstol,reltol=reltol, saveat=saveat, de$MethodOfSteps(de$Tsit5())))
system.time(sol_dct <- de$solve(prob_dct, abstol=abstol,reltol=reltol, saveat=saveat, de$MethodOfSteps(de$Tsit5())))
udf <- as.data.frame(cbind(sol$t,t(sapply(sol$u,identity)),t(sapply(sol_dct$u,identity))))
colnames(udf) <- c("t","P","S_1","E","I_1","Pdct","S_1dct","Edct","I_1dct")

library(plotly)
fig <- plotly::plot_ly(udf, x = sol$t, y = ~S_1, name="S<sub>1</sub>",type = 'scatter', mode = 'lines') %>% plotly::add_trace(y = ~E, name="E") %>% plotly::add_trace(y = ~P,name="P") %>% plotly::add_trace(y = ~I_1,name="I<sub>1</sub>") %>%
  layout(
    title = list(text="Prediction for inhibitor volume 3 mkM",y=0.95),
    xaxis = list(title = "time, secs"),
    yaxis = list(title = "volume, ml"),
    font = font,
    legend = list(x = 0.7, y = 0.9,font = list(size = 20))
  )
orca(fig, file = paste(paste("Prediction_3"),"png", sep = "."))


# Comparisons

library(plotly)

P<-udf$P[which(udf$t %in% df$time_sec)]
Pdct<-udf$Pdct[which(udf$t %in% df$time_sec)]
udf_final<-cbind(df,P)
udf_final<-cbind(udf_final,Pdct)
udf_final$K=unlist(parest)[22]*udf_final$P - unlist(parest)[23]*'^'(udf_final$P,3/2)
udf_final$Kdct=unlist(parest)[22]*udf_final$Pdct - unlist(parest)[23]*'^'(udf_final$Pdct,3/2)

fig <- plotly::plot_ly(udf_final, x = udf_final$time_sec, y = udf_final$K, type = 'scatter', name="&#954;<sub>pred</sub>", mode = 'lines') %>% 
  plotly::add_trace(udf_final, x = udf_final$time_sec, y = udf_final$inhibitor_1, name="&#954;<sub>exp</sub>", mode='markers',marker = list(size = 3)) %>% 
  layout(
    title = list(text="Conductivity",y=0.95),
    xaxis = list(title = "time, secs"),
    yaxis = list(title = "muS/cm"),
    font=font,
    legend = list(x = 0.7, y = 0.9,font = list(size = 20))
  ) %>%
  config(mathjax = "cdn")
orca(fig, file = paste(paste("P_pred_exp_3"),"png", sep = "."))

fig <- plotly::plot_ly(udf_final, x = udf_final$time_sec, y = udf_final$K, type = 'scatter', name="&#954;<sub>pred,distributed</sub>", mode = 'lines') %>% 
  plotly::add_trace(udf_final, x = udf_final$time_sec, y = udf_final$Kdct, type = 'scatter', name="&#954;<sub>pred,discrete</sub>", mode = 'lines') %>% 
  plotly::add_trace(udf_final, x = udf_final$time_sec, y = udf_final$inhibitor_3, name="&#954;<sub>exp</sub>", mode='markers',marker = list(size = 3)) %>% 
  layout(
    title = list(text="Conductivity for 3 mkM of inhibitor",y=0.95,yref="container"),
    xaxis = list(title = "time, secs"),
    yaxis = list(title = "&mu; S/cm"),
    font=font,
    legend = list(x = 0.7, y = 0.9,font = list(size = 20))
  ) %>%
  config(mathjax = "cdn")
orca(fig, file = paste(paste("P_pred_exp_3_discrete_distributed"),"png", sep = "."))


fig <- plotly::plot_ly(udf_final,x = udf_final$time_sec, y = abs(udf_final$K - udf_final$inhibitor_3), type = 'scatter', name="&#954;<sub>pred,distributed</sub>", mode = 'lines') %>% 
  plotly::add_trace(udf_final,x = udf_final$time_sec, y = abs(udf_final$Kdct - udf_final$inhibitor_3), type = 'scatter', name="&#954;<sub>pred,discrete</sub>", mode = 'lines') %>% 
  layout(
    title = list(text="Errors of models for 3 mkM of inhibitor",y=0.95,yref="container"),
    xaxis = list(title = "time, secs"),
    yaxis = list(title = "&mu; S/cm"),
    font=font,
    legend = list(x = 0.7, y = 0.9,font = list(size = 20))
  ) %>%
  config(mathjax = "cdn")
orca(fig, file = paste(paste("P_pred_exp_3_error_discrete_distributed"),"png", sep = "."))

##########################################
# I_1_0=5 ################################
##########################################
u0 <- unlist(parest[17:20]) 
u0[4] <- 5
#tspan <- c(max(0,min(df$time_sec)),max(0,max(df$time_sec)))
p <- parest
constant_lags <- c(100.0)
JuliaCall::julia_assign("u0", u0)
JuliaCall::julia_assign("tspan", tspan)
JuliaCall::julia_assign("constant_lags", constant_lags)
JuliaCall::julia_assign("p", p)
prob <- JuliaCall::julia_eval("DDEProblem(f, u0, h, tspan, constant_lags = constant_lags, p)")
prob_dct <- JuliaCall::julia_eval("DDEProblem(f_dct, u0, h, tspan, constant_lags = constant_lags, p)")
system.time(sol <- de$solve(prob, abstol=abstol,reltol=reltol, saveat=saveat, de$MethodOfSteps(de$Tsit5())))
system.time(sol_dct <- de$solve(prob_dct, abstol=abstol,reltol=reltol, saveat=saveat, de$MethodOfSteps(de$Tsit5())))
udf <- as.data.frame(cbind(sol$t,t(sapply(sol$u,identity)),t(sapply(sol_dct$u,identity))))
colnames(udf) <- c("t","P","S_1","E","I_1","Pdct","S_1dct","Edct","I_1dct")

library(plotly)
fig <- plotly::plot_ly(udf, x = sol$t, y = ~S_1, name="S<sub>1</sub>",type = 'scatter', mode = 'lines') %>% plotly::add_trace(y = ~E, name="E") %>% plotly::add_trace(y = ~P,name="P") %>% plotly::add_trace(y = ~I_1,name="I<sub>1</sub>") %>%
  layout(
    title = list(text="Prediction for inhibitor volume 5 mkM",y=0.95),
    xaxis = list(title = "time, secs"),
    yaxis = list(title = "volume, ml"),
    font = font,
    legend = list(x = 0.7, y = 0.9,font = list(size = 20))
  )
orca(fig, file = paste(paste("Prediction_5"),"png", sep = "."))


# Comparisons

library(plotly)

P<-udf$P[which(udf$t %in% df$time_sec)]
Pdct<-udf$Pdct[which(udf$t %in% df$time_sec)]
udf_final<-cbind(df,P)
udf_final<-cbind(udf_final,Pdct)
udf_final$K=unlist(parest)[22]*udf_final$P - unlist(parest)[23]*'^'(udf_final$P,3/2)
udf_final$Kdct=unlist(parest)[22]*udf_final$Pdct - unlist(parest)[23]*'^'(udf_final$Pdct,3/2)

fig <- plotly::plot_ly(udf_final, x = udf_final$time_sec, y = udf_final$K, type = 'scatter', name="&#954;<sub>pred</sub>", mode = 'lines') %>% 
  plotly::add_trace(udf_final, x = udf_final$time_sec, y = udf_final$inhibitor_1, name="&#954;<sub>exp</sub>", mode='markers',marker = list(size = 3)) %>% 
  layout(
    title = list(text="Conductivity",y=0.95),
    xaxis = list(title = "time, secs"),
    yaxis = list(title = "muS/cm"),
    font=font,
    legend = list(x = 0.7, y = 0.9,font = list(size = 20))
  ) %>%
  config(mathjax = "cdn")
orca(fig, file = paste(paste("P_pred_exp_5"),"png", sep = "."))

fig <- plotly::plot_ly(udf_final, x = udf_final$time_sec, y = udf_final$K, type = 'scatter', name="&#954;<sub>pred,distributed</sub>", mode = 'lines') %>% 
  plotly::add_trace(udf_final, x = udf_final$time_sec, y = udf_final$Kdct, type = 'scatter', name="&#954;<sub>pred,discrete</sub>", mode = 'lines') %>% 
  plotly::add_trace(udf_final, x = udf_final$time_sec, y = udf_final$inhibitor_5, name="&#954;<sub>exp</sub>", mode='markers',marker = list(size = 3)) %>% 
  layout(
    title = list(text="Conductivity for 5 mkM of inhibitor",y=0.95,yref="container"),
    xaxis = list(title = "time, secs"),
    yaxis = list(title = "&mu; S/cm"),
    font=font,
    legend = list(x = 0.7, y = 0.9,font = list(size = 20))
  ) %>%
  config(mathjax = "cdn")
orca(fig, file = paste(paste("P_pred_exp_5_discrete_distributed"),"png", sep = "."))


fig <- plotly::plot_ly(udf_final,x = udf_final$time_sec, y = abs(udf_final$K - udf_final$inhibitor_5), type = 'scatter', name="&#954;<sub>pred,distributed</sub>", mode = 'lines') %>% 
  plotly::add_trace(udf_final,x = udf_final$time_sec, y = abs(udf_final$Kdct - udf_final$inhibitor_5), type = 'scatter', name="&#954;<sub>pred,discrete</sub>", mode = 'lines') %>% 
  layout(
    title = list(text="Errors of models for 5 mkM of inhibitor",y=0.95,yref="container"),
    xaxis = list(title = "time, secs"),
    yaxis = list(title = "&mu; S/cm"),
    font=font,
    legend = list(x = 0.7, y = 0.9,font = list(size = 20))
  ) %>%
  config(mathjax = "cdn")
orca(fig, file = paste(paste("P_pred_exp_5_error_discrete_distributed"),"png", sep = "."))

##########################################
# I_1_0=10 ################################
##########################################
u0 <- unlist(parest[17:20]) 
u0[4] <- 10
#tspan <- c(max(0,min(df$time_sec)),max(0,max(df$time_sec)))
p <- parest
constant_lags <- c(100.0)
JuliaCall::julia_assign("u0", u0)
JuliaCall::julia_assign("tspan", tspan)
JuliaCall::julia_assign("constant_lags", constant_lags)
JuliaCall::julia_assign("p", p)
prob <- JuliaCall::julia_eval("DDEProblem(f, u0, h, tspan, constant_lags = constant_lags, p)")
prob_dct <- JuliaCall::julia_eval("DDEProblem(f_dct, u0, h, tspan, constant_lags = constant_lags, p)")
system.time(sol <- de$solve(prob, abstol=abstol,reltol=reltol, saveat=saveat, de$MethodOfSteps(de$Tsit5())))
system.time(sol_dct <- de$solve(prob_dct, abstol=abstol,reltol=reltol, saveat=saveat, de$MethodOfSteps(de$Tsit5())))
udf <- as.data.frame(cbind(sol$t,t(sapply(sol$u,identity)),t(sapply(sol_dct$u,identity))))
colnames(udf) <- c("t","P","S_1","E","I_1","Pdct","S_1dct","Edct","I_1dct")

library(plotly)
fig <- plotly::plot_ly(udf, x = sol$t, y = ~S_1, name="S<sub>1</sub>",type = 'scatter', mode = 'lines') %>% plotly::add_trace(y = ~E, name="E") %>% plotly::add_trace(y = ~P,name="P") %>% plotly::add_trace(y = ~I_1,name="I<sub>1</sub>") %>%
  layout(
    title = list(text="Prediction for inhibitor volume 10 mkM",y=0.95),
    xaxis = list(title = "time, secs"),
    yaxis = list(title = "volume, ml"),
    font = font,
    legend = list(x = 0.7, y = 0.9,font = list(size = 20))
  )
orca(fig, file = paste(paste("Prediction_10"),"png", sep = "."))


# Comparisons

library(plotly)

P<-udf$P[which(udf$t %in% df$time_sec)]
Pdct<-udf$Pdct[which(udf$t %in% df$time_sec)]
udf_final<-cbind(df,P)
udf_final<-cbind(udf_final,Pdct)
udf_final$K=unlist(parest)[22]*udf_final$P - unlist(parest)[23]*'^'(udf_final$P,3/2)
udf_final$Kdct=unlist(parest)[22]*udf_final$Pdct - unlist(parest)[23]*'^'(udf_final$Pdct,3/2)

fig <- plotly::plot_ly(udf_final, x = udf_final$time_sec, y = udf_final$K, type = 'scatter', name="&#954;<sub>pred</sub>", mode = 'lines') %>% 
  plotly::add_trace(udf_final, x = udf_final$time_sec, y = udf_final$inhibitor_10, name="&#954;<sub>exp</sub>", mode='markers',marker = list(size = 3)) %>% 
  layout(
    title = list(text="Conductivity",y=0.95),
    xaxis = list(title = "time, secs"),
    yaxis = list(title = "muS/cm"),
    font=font,
    legend = list(x = 0.7, y = 0.9,font = list(size = 20))
  ) %>%
  config(mathjax = "cdn")
orca(fig, file = paste(paste("P_pred_exp_10"),"png", sep = "."))

fig <- plotly::plot_ly(udf_final, x = udf_final$time_sec, y = udf_final$K, type = 'scatter', name="&#954;<sub>pred,distributed</sub>", mode = 'lines') %>% 
  plotly::add_trace(udf_final, x = udf_final$time_sec, y = udf_final$Kdct, type = 'scatter', name="&#954;<sub>pred,discrete</sub>", mode = 'lines') %>% 
  plotly::add_trace(udf_final, x = udf_final$time_sec, y = udf_final$inhibitor_10, name="&#954;<sub>exp</sub>", mode='markers',marker = list(size = 3)) %>% 
  layout(
    title = list(text="Conductivity for 10 mkM of inhibitor",y=0.95,yref="container"),
    xaxis = list(title = "time, secs"),
    yaxis = list(title = "&mu; S/cm"),
    font=font,
    legend = list(x = 0.7, y = 0.9,font = list(size = 20))
  ) %>%
  config(mathjax = "cdn")
orca(fig, file = paste(paste("P_pred_exp_10_discrete_distributed"),"png", sep = "."))


fig <- plotly::plot_ly(udf_final,x = udf_final$time_sec, y = abs(udf_final$K - udf_final$inhibitor_10), type = 'scatter', name="&#954;<sub>pred,distributed</sub>", mode = 'lines') %>% 
  plotly::add_trace(udf_final,x = udf_final$time_sec, y = abs(udf_final$Kdct - udf_final$inhibitor_10), type = 'scatter', name="&#954;<sub>pred,discrete</sub>", mode = 'lines') %>% 
  layout(
    title = list(text="Errors of models for 10 mkM of inhibitor",y=0.95,yref="container"),
    xaxis = list(title = "time, secs"),
    yaxis = list(title = "&mu; S/cm"),
    font=font,
    legend = list(x = 0.7, y = 0.9,font = list(size = 20))
  ) %>%
  config(mathjax = "cdn")
orca(fig, file = paste(paste("P_pred_exp_10_error_discrete_distributed"),"png", sep = "."))


########################################################
# delay distribution functions #########################
########################################################
a=as.numeric(parest[5])
m=as.numeric(parest[6])
tau_m=as.numeric(parest[7])
rmse=sqrt((m+1)/a^2)
psi_func = function(x,a,m,tau_m) {(a^(m+1)/gamma(m+1))*(x-tau_m)^m * exp(-a*(x-tau_m))}
arguments = seq(tau_m,tau_m+(m+1)/a+3*rmse,0.05)

data<-data.frame(time=c(arguments),psi=psi_func(c(arguments),a,m,tau_m))


png(paste(paste("ESIP_delay_densities_1"),"png", sep = "."), width = 800, height = 600)
densities <- ggplot()+geom_line(data=data,aes(x=time,y=psi),size=1.5)+
  xlab("t, secs") + ylab(expression(italic(psi))) + 
  my_theme
plot(densities)
dev.off()

save(list=c("parest","parms","lower","upper"),file = paste("S_0dot3",".Rdata",sep = "."))

save.image(file="ESIP_BrownModelGammaDist_multiple_volumes.RData")