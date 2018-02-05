#############################
## SIWR Example - FIM & AS ##
## Marisa Eisenberg        ##
## 2-3-18                  ##
#############################

# Notes
# - Generates local and global active subspace metrics using both the vector and scalar QOIs (either time series or cost function)
# - Uses the scaled form of the parameters
# - Global measures sample over a ±percent range



#### Load libraries ####
library(deSolve)
library(Matrix)
library(ggplot2)
library(gridExtra)
library(plotly)
library(reshape2)
library(lhs)
library(MASS)
source('sFIMscaledzeroparams.R')
source('ProfParam.R')

#### Parameters & setup ####
params = c('betaI'=0, 'betaW'=0, 'xi' = 0, 'k' = 0) # base parameters from Eisenberg 2013

i=1 # set to determine which value of xi to use
xifactors = c(1,5000)
timefactors = c(150,17)
times = seq(0,timefactors[i],length.out=20)

paramstemp = params
paramstemp[3] = paramstemp[3]*xifactors[i] # set xi to current value

#### Model functions and simulate data for scalar QOI ####

# Model equations
SIWRode = function(t, x, params){
  
  S = x[1]
  I = x[2]
  W = x[3]
  R = x[4]
  
  betaI = (1+params[1])*0.256
  betaW = (1+params[2])*1.21
  gamma = 0.25
  xi = (1+params[3])*0.00756
  
  dSdt = - betaI*I*S - betaW*W*S
  dIdt = betaI*I*S + betaW*W*S - gamma*I
  dWdt = xi*(I-W)
  dRdt = gamma*I
  
  list(c(dSdt, dIdt, dWdt, dRdt))
}

# Function to calculate the initial conditions from wherever we are in parameter space
x0fun = function(cases,params) {
  x0 = c(1-(cases[1]*(1+params[4])*1.1212e-5), cases[1]*(1+params[4])*1.1212e-5, 0, 0)
  names(x0) = c('S0','I0','W0', 'R0')
  x0}

# Function to calculate the model output y = k*I
yfun = function(odeSim, params){odeSim[,3]/((1+params[4])*1.1212e-5)} 


# Simulate data
xtemp = ode(x0fun(113,paramstemp), times, SIWRode, paramstemp, method='ode45')
ytemp = yfun(xtemp,paramstemp)

# Least squares cost function using above data
SIWRLS=function(p){
  # print(p)
  
  # Simulate model
  xcurr = ode(x0fun(113,p), times, SIWRode, p, method='ode45')
  
  # Measurement equation
  y = yfun(xcurr,p)
  
  # Sum of squares
  RSS =  sum((y - ytemp)^2)
}

# R0 function (in case we want that as our QOI instead)
R0 = function(p){((1+p[1])*0.256+(1+p[2])*1.21)/0.25}

#### Active Subspace Calculations ####
fmapscalar = SIWRLS
set.seed(8)
numsamples = 500
# paramsample = (randomLHS(numsamples,length(params)) - 0.5)*2
paramsample = (randomLHS(numsamples,length(params)) - 0.5)
# paramsample = (randomLHS(numsamples,length(params)) - 0.5)*0.5
# paramsample = (randomLHS(numsamples,length(params)) - 0.5)*0.2


# Local FIMs (F) - vector & scalar forms
FIM.loc.vec = sFIM(times,paramstemp,x0fun,SIWRode,yfun,ytemp)
FIM.loc.sca = grad(fmapscalar,paramstemp,method="Richardson") %o% grad(fmapscalar,paramstemp,method="Richardson")


# Global FIMs (C) - vector & scalar forms
FIMsample.vec = list()
FIMsample.sca = list()
fscalarsample = c()

C.vec = matrix(0,length(params),length(params))
C.sca = matrix(0,length(params),length(params))

for(i in 1:numsamples){
  fscalarsample = c(fscalarsample,fmapscalar(paramsample[i,]))
  
  FIMsample.vec[[i]] = sFIM(tspan = times,params = paramsample[i,],x0fcn = x0fun,xfcn = SIWRode,yfcn = yfun,data = ytemp)
  FIMsample.sca[[i]] = grad(fmapscalar,paramsample[i,],method="Richardson") %o% grad(fmapscalar,paramsample[i,],method="Richardson")
  
  C.vec = C.vec + FIMsample.vec[[i]]/numsamples
  C.sca = C.sca + FIMsample.sca[[i]]/numsamples
}

print(rankMatrix(FIM.loc.vec)[1])
print(rankMatrix(FIM.loc.sca)[1])

print(rankMatrix(C.vec)[1])
print(rankMatrix(C.sca)[1])

F.vec.decomp = eigen(FIM.loc.vec,symmetric=TRUE)
F.vec.eigenvalues = F.vec.decomp$values
F.vec.eigenvectors = F.vec.decomp$vectors

C.vec.decomp = eigen(C.vec,symmetric=TRUE)
C.vec.eigenvalues = C.vec.decomp$values
C.vec.eigenvectors = C.vec.decomp$vectors

F.sca.decomp = eigen(FIM.loc.sca,symmetric=TRUE)
F.sca.eigenvalues = F.sca.decomp$values
F.sca.eigenvectors = F.sca.decomp$vectors

C.sca.decomp = eigen(C.sca,symmetric=TRUE)
C.sca.eigenvalues = C.sca.decomp$values
C.sca.eigenvectors = C.sca.decomp$vectors

# Takehome from eigenvectors: same sort of basic trend in all four, but different ordering for the eigenvalues

#### Plots ####
eigenvalues = list(F.vec.eigenvalues, C.vec.eigenvalues, F.sca.eigenvalues, C.sca.eigenvalues)
eigenvectors = list(F.vec.eigenvectors, C.vec.eigenvectors, F.sca.eigenvectors, C.sca.eigenvectors)
names = list("-Fvec","-Cvec","-Fsca","-Csca")

# eigenvalue plots
for(i in 1:length(eigenvalues)){
  df_eig=as.data.frame(cbind(rep(0,4),eigenvalues[[i]]))
  names(df_eig) = c("x","y")
  eigplot = ggplot(df_eig,aes(x,y,xend=x+1,yend=y))+
    geom_segment()+
    scale_shape_identity()+
    theme_bw(base_size = 8) + labs(x="", y="Eigenvalues")+
    theme(panel.border = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line.x = element_line(color = "black",size=0.4),
          axis.line.y = element_line(color = "black",size=0.4))+
    scale_x_discrete()+
    # scale_y_continuous(trans="log",breaks=c(1E-25,1E-20,1E-15,1E-10,1E-5,1E3,1E4,1E5,1E6,1E7,1E8,1E9,1E10,1E11,1E12,1E13,1E14,1E15,1E16,1E17,1E18,1E19,1E20,1E21,1E22))
    # scale_y_continuous(trans="log",breaks=c(1E-27,1E-24,1E-21,1E-18,1E-15,1E-12,1E-9,1E-6,1E-3,1E2,1E4,1E6,1E8,1E10,1E12,1E14,1E16,1E18,1E20,1E22,1E24,1E26,1E28,1E30))
    scale_y_continuous(trans="log",breaks=c(1E-30,1E-25,1E-20,1E-15,1E-10,1E-5,1E3,1E4,1E5,1E6,1E7,1E8,1E9,1E10,1E11,1E12,1E13,1E14,1E15,1E16,1E17,1E18,1E19,1E20,1E21,1E22))
  eigplot
  ggsave(file=paste("eigplot",names[[i]],".pdf",sep=""), width = 1, height = 3.25, eigplot)
}


# sufficient summary plots
labels = list(expression(paste(Q["a,1"]^T,theta)), expression(paste(Q["a,2"]^T,theta)), expression(paste(Q["a,3"]^T,theta)), expression(paste(Q["a,4"]^T,theta)))

for(i in 1:length(eigenvalues)){
  for(j in 1:4){
    component=t(t(eigenvectors[[i]][,j])%*%t(paramsample))
    df=as.data.frame(cbind(component,fscalarsample))
    names(df) = c("x","y")
    eigsummary = ggplot(data=df,aes(x,y,color=y))+
      geom_point()+
      #scale_color_gradient2(low="dodgerblue2", high="firebrick2", mid="palegoldenrod", midpoint=25,space="Lab",guide =guide_colorbar(title="Cost function value"))+
      scale_color_gradientn(guide =guide_colorbar(title="Cost function value"),colors=c("dodgerblue1","lightgreen","gold","firebrick1"),values=c(0,0.25,0.5,1),limits=c(min(fscalarsample),max(fscalarsample)))+
      theme_bw(base_size = 8) + labs(x=labels[[j]], y='Cost function value')+
      theme(axis.title.y=element_text(vjust=1.5))+
      #theme(legend.position = c(0.75, 0.9),legend.key = element_blank(),legend.background=element_blank())+
      theme(panel.border = element_blank(),
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line.x = element_line(color = "black",size=0.4),
            axis.line.y = element_line(color = "black",size=0.4))+
      scale_x_continuous()+
      scale_y_continuous()
    
    eigsummary
    ggsave(file=paste("eig",j,"summary",names[[i]],".pdf",sep=""), width = 3.25, height = 3.25, eigsummary)
  }
}


# qvq plots, for eig 1 vs eig 2

for(i in 1:length(eigenvalues)){
  component1=t(t(eigenvectors[[i]][,1])%*%t(paramsample))
  component2=t(t(eigenvectors[[i]][,2])%*%t(paramsample))
  df=as.data.frame(cbind(component1,component2,fscalarsample))
  names(df) = c("x","y","z")
  
  eig1eig2 = ggplot(data=df,aes(x,y,color=z))+
    geom_point()+
    #scale_color_gradientn(colors=terrain.colors(10))+
    #scale_color_gradient2(low="dodgerblue2", high="firebrick2", mid="palegoldenrod", midpoint=25,space="Lab",guide =guide_colorbar(title="Cost function value"))+
    scale_color_gradientn(guide =guide_colorbar(title="Cost function value"),colors=c("dodgerblue1","lightgreen","gold","firebrick1"),values=c(0,0.25,0.5,1),limits=c(min(fscalarsample),max(fscalarsample)))+
    theme_bw(base_size = 8) + labs(x=expression(paste(Q["a,1"]^T,theta)), y=expression(paste(Q["a,2"]^T,theta)))+
    theme(axis.title.y=element_text(vjust=1.5))+
    theme(panel.border = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line.x = element_line(color = "black",size=0.4),
          axis.line.y = element_line(color = "black",size=0.4))+
    scale_x_continuous()+
    scale_y_continuous() 
  
  eig1eig2
  ggsave(file=paste("eig1eig2",names[[i]],".pdf",sep=""), width = 3.25, height = 3.25, eig1eig2)
}



# eigenvector loadings
labels = list("first","second","third","fourth")

for(i in 1:length(eigenvalues)){
  for(j in 1:4){
    df_eigvec = data.frame(cbind(1:4,eigenvectors[[i]][,j]))
    names(df_eigvec) = c("x","y")
    eigvec = ggplot(data=df_eigvec,aes(x,y=abs(y),fill=abs(y)))+
      geom_bar(stat = "identity")+
      scale_fill_gradient(low="dodgerblue1", high="firebrick1",space="Lab",limits=c(0, 1), guide =FALSE)+
      theme_bw(base_size = 14) + labs(x="Parameter", y=paste("Magnitude in the", labels[[j]], "eigenvector"))+
      theme(axis.title.y=element_text(vjust=1.5))+
      theme(panel.border = element_blank(),
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line.x = element_line(color = "black",size=0.4),
            axis.line.y = element_line(color = "black",size=0.4))+
      scale_x_continuous(labels=c(expression(beta[I]),expression(beta[W]),expression(xi),'k'), breaks = c(1,2,3,4))+
      scale_y_continuous(limits=c(0,1))
    eigvec
    ggsave(file=paste("eig",j,"load", names[[i]],".pdf",sep=""), width = 3.25, height = 3.25, eigvec)
  }
}
