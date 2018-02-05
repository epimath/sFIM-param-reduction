library(ggplot2)
library(ggpubr)
library(lhs)
library(deSolve)
library(MASS)
library(Bhat)
library(numDeriv)
library(compiler)
library(reshape2)

enableJIT(3)


dq=function(u,y,params) #function(time, state, parameters)
{
  Cdc20tot=params[1]
  E2Ftot=params[2]
  GF=params[3]
  Kda=params[4]
  Kdb=params[5]
  Kdd=params[6]
  Kde=params[7]
  Kgf=params[8]
  K1cdc20=params[9]
  K2cdc20=params[10]
  K1e2f=params[11]
  K2e2f=params[12]
  Vda=params[13]
  Vdb=params[14]
  Vdd=params[15]
  Vde=params[16]
  vsa=params[17]
  vsb=params[18]
  vsd=params[19]
  vse=params[20]
  V1cdc20=params[21]
  V2cdc20=params[22]
  V1e2f=params[23]
  V2e2f=params[24]
  
  Mddot = vsd * GF/(Kgf+GF) - Vdd * y[1]/(Kdd+y[1])
  e2fdot = V1e2f * ( E2Ftot - y[2])/(K1e2f + ( E2Ftot - y[2])) * (y[1] + y[3]) - V2e2f * y[2]/(K2e2f+y[2])*y[4]
  Medot = vse*y[2] - Vde*y[4]*y[3]/(Kde+y[3])
  Madot = vsa*y[2] - Vda*y[6]*y[4]/(Kda+y[4])
  Mbdot = vsb*y[4] - Vdb*y[6]*y[5]/(Kdb+y[5])
  Cdc20dot= V1cdc20*y[5] * (Cdc20tot-y[6])/(K1cdc20+(Cdc20tot-y[6])) - V2cdc20 * y[6]/(K2cdc20+y[6])
  
  
  return(list(c(Mddot,e2fdot,Medot,Madot,Mbdot,Cdc20dot)))
}


params=c(5,3,1,0.1,0.005,0.1,0.1,0.1,1,1,0.01,0.01,0.245,0.28,0.245,0.35,0.175,0.21,0.175,0.21,0.21, 0.35,0.805,0.7)

dt=0.05

# out=ode(y=c(Y1=0.1851852,Y2=2.159831,Y3=1.431198,Y4=3.093866,Y5=1.633533,y6=0.630261),seq(0,100,dt),dq, params, method='lsoda',atol=1e-15)
# plot(out[,1],out[,4],col="blue4",type="l",ylim=c(0,3.5))
# lines(out[,1],out[,5],col="green4",type="l")
# lines(out[,1],out[,6],col="red4",type="l")

#df_period=as.data.frame(cbind(rep(out[,1],3),c(out[,4],out[,5],out[,6]),c(rep("Cyclin E-Cdk2",length(out[,1])),rep("Cyclin A-Cdk2",length(out[,1])),rep("Cyclin B-Cdk1",length(out[,1])))))
df_period=as.data.frame(cbind(out[,1],out[,4],out[,5],out[,6]))
names(df_period) = c("x","Cyclin E-Cdk2","Cyclin A-Cdk2","Cyclin B-Cdk1")
df_period_melt=melt(df_period,id="x")


p0 = ggplot(data=df_period_melt,aes(x=x,y=value,color=variable))+
  geom_line()+
  theme_bw(base_size = 8) + labs(x="Time (hours)", y=expression(paste("Concentration (",mu,"M)")))+
  theme(axis.title.y=element_text(vjust=1.5))+
  scale_color_discrete(guide=guide_legend(direction="horizontal",title=""))+
  theme(legend.position = "top",legend.key = element_blank(),legend.background=element_blank())+
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color = "black",size=0.4),
        axis.line.y = element_line(color = "black",size=0.4))+
  scale_x_continuous()+
  scale_y_continuous()

p0
ggsave(file="Fig3e.pdf", width = 3.25, height = 3.25)


# diff(tail(which(diff(sign(diff(out[,4])))==-2)+1,n=2))*dt


f = function(par){
  out=ode(y=c(Y1=0.1,Y2=0.1,Y3=0.1,Y4=0.1,Y5=0.1,y6=0.1),seq(0,200,dt),dq, par, method='lsoda',atol=1e-15)
  diff(tail(which(diff(sign(diff(out[,6])))==-2)+1,n=2))*dt
}

# f(params)

N=10^3
set.seed(0)
params_matrix=matrix(params,nrow=N,ncol=24,byrow=TRUE)
LHS=randomLHS(n=N,k=24)
sample=params_matrix+ 0.5*(LHS-0.5)*params_matrix
# saveRDS(sample,"sample.rds")

all_gradients=matrix(NA,N,24)
period=matrix(NA,N,1)

for (i in 1:N){
  print(i)
  all_gradients[i,]= grad(f,sample[i,],method="Richardson")
  period[i,]= f(sample[i,])
}

# saveRDS(all_gradients,"all_gradients.rds")
# saveRDS(period,"period.rds")

##########

sample=readRDS("sample.rds")
all_gradients=readRDS("all_gradients.rds")
period=readRDS("period.rds")
N=10^3

C=matrix(0,24,24)

for (i in 1:N){
  C=C+1/N * all_gradients[i,] %o% all_gradients[i,]
}

decomp=eigen(C,symmetric=TRUE)
eigenvalues=decomp$values
eigenvectors=decomp$vectors

component1=t(t(eigenvectors[,1])%*%t(sample))
component2=t(t(eigenvectors[,2])%*%t(sample))

#plot(component1,period)
#plot(component2,period)

df_eig=as.data.frame(cbind(rep(0,24),eigenvalues))
names(df_eig) = c("x","y")


p1 = ggplot(df_eig,aes(x,y,xend=x+1,yend=y))+
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
  scale_y_continuous(trans="log",breaks=c(1E-5,1E-4,1E-3,1E-2,1E-1,1E0,1E1,1E2,1E3,1E4,1E5,1E6,1E7,1E8))
p1
ggsave(file="Fig3a.pdf", width = 1, height = 3.25)


df=as.data.frame(cbind(component1,component2,period))
names(df) = c("x","y","z")

p2 = ggplot(data=df,aes(x,z,color=z))+
  geom_point()+
  #scale_color_gradient2(low="dodgerblue2", high="firebrick2", mid="palegoldenrod", midpoint=25,space="Lab",guide =guide_colorbar(title="Period"))+
  scale_color_gradientn(guide =guide_colorbar(title="Period"),colors=c("dodgerblue1","lightgreen","gold","firebrick1"),values=c(0,0.25,0.5,1),limits=c(min(period),max(period)))+
  theme_bw(base_size = 8) + labs(x=expression(paste(Q["a,1"]^T,theta)), y='Period (hours)')+
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

p2
ggsave(file="Fig3b.pdf", width = 3.25, height = 3.25)



p3 = ggplot(data=df,aes(x,y,color=z))+
  geom_point()+
  #scale_color_gradientn(colors=terrain.colors(10))+
  #scale_color_gradient2(low="dodgerblue2", high="firebrick2", mid="palegoldenrod", midpoint=25,space="Lab",guide =guide_colorbar(title="Period"))+
  scale_color_gradientn(guide =guide_colorbar(title="Period"),colors=c("dodgerblue1","lightgreen","gold","firebrick1"),values=c(0,0.25,0.5,1),limits=c(min(period),max(period)))+
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
  
p3
ggsave(file="Fig3c.pdf", width = 3.25, height = 3.25)


df_eigvec1= data.frame(cbind(1:24,eigenvectors[,1]))
names(df_eigvec1)=c("x","y")

p4 = ggplot(data=df_eigvec1,aes(x,y=abs(y),fill=abs(y)))+
  geom_bar(stat="identity")+
  scale_fill_gradient(low="dodgerblue1", high="firebrick1",space="Lab",limits=c(0, 0.63), guide =FALSE)+
  theme_bw(base_size = 8) + labs(x="Parameter", y="Magnitude in first eigenvector")+
  theme(axis.title.y=element_text(vjust=1.5))+
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color = "black",size=0.4),
        axis.line.y = element_line(color = "black",size=0.4))+
  scale_x_continuous(labels=c())+
  scale_y_continuous(limits=c(0,0.7))+
  annotate("text",x=8,y=0.45,label="Degredation of CycB-Cdk1",size = 1.75)+
  annotate("text",x=10,y=0.55,label="Degredation of CycE-Cdk2",size = 1.75)+
  annotate("text",x=10,y=0.65,label="Synthesis of CycA-Cdk2",size = 1.75)+
  annotate("text",x=20,y=0.7,label="Synthesis of CycE-Cdk2",size = 1.75)+
  annotate("text",x=22.5,y=0.6,label="Activation of Cdc20",size = 1.75)+
  annotate("segment",x=21,xend=24,y=0.35,yend=0.57)+
  annotate("segment",x=20,xend=18,y=0.53,yend=0.67)+
  annotate("segment",x=17,xend=14.5,y=0.64,yend=0.65)+
  annotate("segment",x=16,xend=14,y=0.28,yend=0.52)+
  annotate("segment",x=14,xend=10,y=0.29,yend=0.42)
  

p4

df_eigvec2= data.frame(cbind(1:24,eigenvectors[,2]))
names(df_eigvec2)=c("x","y")

p5 = ggplot(data=df_eigvec2,aes(x,y=abs(y),fill=abs(y)))+
  geom_bar(stat="identity")+
  scale_fill_gradient(low="dodgerblue1", high="firebrick1",space="Lab",limits=c(0, 0.65),guide =FALSE)+
  theme_bw(base_size = 8) + labs(x="Parameter", y="Magnitude in second eigenvector")+
  theme(axis.title.y=element_text(vjust=1.5))+
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color = "black",size=0.4),
        axis.line.y = element_line(color = "black",size=0.4))+
  scale_x_continuous(labels=c())+
  scale_y_continuous(limits=c(0,0.7))+
  annotate("text",x=9,y=0.4,label="Degredation of CycA-Cdk2",size = 1.75)+
  annotate("text",x=12,y=0.5,label="Degredation of CycD-Cdk4/6",size = 1.75)+
  annotate("text",x=14,y=0.6,label="Synthesis of CycB-Cdk1",size = 1.75)+
  annotate("text",x=19,y=0.7,label="Synthesis of CycD-Cdk4/6",size = 1.75)+
  annotate("segment",x=19,xend=19,y=0.64,yend=0.67)+
  annotate("segment",x=18,xend=17,y=0.45,yend=0.57)
  #annotate("segment",x=15,xend=14,y=0.49,yend=0.47)+
  #annotate("segment",x=13,xend=13,y=0.36,yend=0.37)
  

p5

p6=ggarrange(p4, p5,ncol = 1, nrow = 2)
p6

ggsave(file="Fig3d.pdf", width = 3.25, height = 3.25)

