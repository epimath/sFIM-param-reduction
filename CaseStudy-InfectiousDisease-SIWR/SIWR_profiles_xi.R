#############################
## SIWR Example - Profiles ##
## Marisa Eisenberg        ##
## 1-29-18                 ##
#############################

## To Do's:
# - make a function to do some of the profile & FIM stuff that we're repeating everywhere
# - clean up notation throughout and make consistent with paper

#### Load libraries ####
library(deSolve)
library(Matrix)
library(ggplot2)
library(gridExtra)
# library(plotly)
library(reshape2)
library(lhs)
library(MASS)
source('sFIM.R')
source('ProfParam.R')

#### Model equations ####

SIWRode = function(t, x, params){
  S = x[1]
  I = x[2]
  W = x[3]
  R = x[4]
  
  betaI = params[1]
  betaW = params[2]
  gamma = 0.25
  xi = params[3]
  
  dSdt = - betaI*I*S - betaW*W*S
  dIdt = betaI*I*S + betaW*W*S - gamma*I
  dWdt = xi*(I-W)
  dRdt = gamma*I
  
  list(c(dSdt, dIdt, dWdt, dRdt))
}

#### Parameters & setup ####
params = c('betaI'=0.256, 'betaW'=1.21, 'xi' = 0.00756, 'k' = 1.1212e-5) # from Eisenberg 2013


# Function to calculate the initial conditions from wherever we are in parameter space
x0fun = function(cases,params) {
  x0 = c(1-(cases[1]*params[4]), cases[1]*params[4], 0, 0)
  names(x0) = c('S0','I0','W0', 'R0')
  x0}


# Function to calculate the model output y = k*I, but params[4] = 1/k
yfun = function(odeSim, params){odeSim[,3]/params[4]} 


# Least squares cost function
SIWRLS=function(params,times,data){
  params = abs(params)
  # Simulate model
  xcurr = ode(x0fun(data,params), times, SIWRode, params, method='ode45')
  
  # Measurement equation
  y = yfun(xcurr,params)
  
  # Sum of squares
  RSS =  sum((y - data)^2)
}

#### Plot fit with default parameters ####
# Data
times = c(0, 7, 14, 21, 28, 35, 42, 49, 56, 63, 70, 77, 84, 91, 98, 105, 112, 119, 126, 133, 140, 147, 154, 161)
cases  = c(113, 60, 75, 148, 379, 2911, 4572, 5361, 5300, 6348, 5346, 4412, 3558, 2271, 1931, 2251, 1692, 1184, 816, 748, 770, 522, 553, 379)

xfit = ode(x0fun(cases,params) + c(0,0,0.00005,0), times, SIWRode, params, method='ode45')

# Measurement equation
yfit = yfun(xfit,params)

modelfit = ggplot() + 
  geom_line(aes(x=times, y=yfit), color = 'black') + 
  geom_point(aes(x=times, y=cases)) + labs(x="Time (days)", y="Cases") + 
  theme_bw(base_size = 12) +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color = "black",size=0.6),
        axis.line.y = element_line(color = "black",size=0.6))

modelfit
ggsave(file="SIWRfit.pdf", width = 3, height = 3)

#### Profiles as xi increases ####

# this section is full of hard-coded junk that needs to get cleaned up but we can do it later.

# Parameter to profile
j = 2

# xifactors = c(1,1000)
# timefactors = c(150,20)
xifactors = c(1,1000,5000)
timefactors = c(150,20,17)
# xifactors = c(1,1000,5000,10000) # runs slow if you do profile, but useful for showing that rank decays to 2
# timefactors = c(150,20,17,17)
# xifactors = c(1,10,100,1000)
# timefactors = c(150,50,30,20)


profiles.xi = list()
profileplotdata = c()
FIMs = list()

for(i in 1:length(xifactors)){
  timestemp = seq(0,timefactors[i],length.out=20)
  
  parameststemp = params
  # parameststemp = paramests
  
  parameststemp[3] = parameststemp[3]*xifactors[i]
  
  xtemp = ode(x0fun(113,parameststemp), timestemp, SIWRode, parameststemp, method='ode45')
  
  # Measurement equation
  yhattemp = yfun(xtemp,parameststemp)
  
  FIMs[[i]] = sFIM(timestemp,parameststemp,x0fun,SIWRode,yfun,113)
  
  # Calculate the rank
  print(rankMatrix(FIMs[[i]])[1])
  
  # profiles.xi[[i]] = ProfParam(paramests,2,SIWRML,timestemp,cases,perrange=0.2)
  profiles.xi[[i]] = ProfParam(parameststemp,j,SIWRLS,timestemp,yhattemp,perrange=0.05,numpoints = 15)
  
  newdata = cbind(profiles.xi[[i]]$profparvals, profiles.xi[[i]]$fnvals)
  colnames(newdata) = c(paste("profparvals",as.character(xifactors[i]),sep=""), paste("fnvals",as.character(xifactors[i]),sep=""))
  
  profileplotdata = cbind(profileplotdata, newdata)
  
  # plot(profiles.xi[[i]]$profparvals, profiles.xi[[i]]$fnvals)
}

threshold = qchisq(0.95,length(parameststemp))/2 + SIWRLS(parameststemp,timestemp,yhattemp)

# profileplot.xi = ggplot(data = as.data.frame(profileplotdata)) +
#   geom_line(aes(x=profparvals1/(parameststemp[j]), y=fnvals1), color='black', linetype=2, size=0.7) +
#   # geom_line(aes(x=profparvals10/(parameststemp[j]), y=fnvals10), color='palegoldenrod') +
#   # geom_line(aes(x=profparvals100/(parameststemp[j]), y=fnvals100), color='dodgerblue2') +
#   # geom_line(aes(x=profparvals1000/(parameststemp[j]), y=fnvals1000), color='black', size=0.7, linetype=3) +
#   geom_line(aes(x=profparvals5000/(parameststemp[j]), y=fnvals5000), color='black', size=0.7,, linetype=1) +
#   geom_line(aes(x=profparvals1000/(parameststemp[j]), y=threshold), color='firebrick2', linetype=3, size=0.7) +
#   labs(x=expression(beta[W]/hat(beta)[W]), y="NLL") + #ylim(0,20) + 
#   theme_bw(base_size = 10) +
#   theme(panel.border = element_blank(),
#         panel.background = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.line.x = element_line(color = "black",size=0.6),
#         axis.line.y = element_line(color = "black",size=0.6))

profileplot.xi = ggplot(data = as.data.frame(profileplotdata)) +
  geom_line(aes(x=profparvals1, y=fnvals1), color='black', linetype=2, size=0.8) +
  # geom_line(aes(x=profparvals10, y=fnvals10), color='palegoldenrod') +
  # geom_line(aes(x=profparvals100, y=fnvals100), color='dodgerblue2') +
  # geom_line(aes(x=profparvals1000, y=fnvals1000), color='black', size=0.7, linetype=3) +
  geom_line(aes(x=profparvals5000, y=fnvals5000), color='black', size=0.8,, linetype=1) +
  geom_line(aes(x=profparvals1000, y=threshold), color='firebrick2', linetype=3, size=0.8) +
  labs(x=expression(beta[W]), y="NLL") + #ylim(0,20) + 
  theme_bw(base_size = 12) +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color = "black",size=0.6),
        axis.line.y = element_line(color = "black",size=0.6))

# profileplot.xi = ggplot(data = as.data.frame(profileplotdata)) +
#   geom_line(aes(x=profparvals1/(paramests[j]*xifactors[1]), y=fnvals1), color='firebrick2') +
#   geom_line(aes(x=profparvals10/(paramests[j]*xifactors[2]), y=fnvals10), color='palegoldenrod') +
#   geom_line(aes(x=profparvals100/(paramests[j]*xifactors[3]), y=fnvals100), color='dodgerblue2') +
#   geom_line(aes(x=profparvals1000/(paramests[j]*xifactors[4]), y=fnvals1000), color='black') +
#   geom_line(aes(x=profparvals1000/(paramests[j]*xifactors[4]), y=threshold), color='black', linetype=2) +
#   labs(x=names(params)[j], y="NLL") + ylim(-5,50) + 
#   theme_bw(base_size = 8) +
#   theme(panel.border = element_blank(),
#         panel.background = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.line.x = element_line(color = "black",size=0.4),
#         axis.line.y = element_line(color = "black",size=0.4))

profileplot.xi
ggsave(file="SIWRprofiles_xilimit.pdf", width = 3, height = 3, profileplot.xi)

#### betaW vs betaI plot for fast xi #####
#redo-ing the profile to get more x-axis range---actually it's fine, nvm

# plot(betaWbetaIprofile$profparvals, betaWbetaIprofile$paramestvals[,1])

i=length(xifactors)
timestemp = seq(0,timefactors[i],length.out=20)
parameststemp = params
parameststemp[3] = parameststemp[3]*xifactors[i]
xtemp = ode(x0fun(113,parameststemp), timestemp, SIWRode, parameststemp, method='ode45')
yhattemp = yfun(xtemp,parameststemp)

# betaWbetaIprofile = ProfParam(parameststemp,2,SIWRLS,timestemp,yhattemp,perrange=0.9,numpoints = 10)
betaWbetaIprofile = profiles.xi[[length(profiles.xi)]]

bWbIplot.xi = ggplot() +
  geom_line(aes(x=betaWbetaIprofile$profparvals, y=betaWbetaIprofile$paramestvals[,1]), color='black', size=1.1) +
  labs(x=expression(beta[W]), y=expression(beta[I])) + #xlim(0,1.5) + ylim(0,1.5) + 
  theme_bw(base_size = 12) +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color = "black",size=0.6),
        axis.line.y = element_line(color = "black",size=0.6))

bWbIplot.xi
ggsave(file="bWbI_xilimit.pdf", width = 3, height = 3, bWbIplot.xi)



