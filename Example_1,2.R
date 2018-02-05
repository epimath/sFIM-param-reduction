library(ggplot2)
library(lhs)

#################

data=matrix(NA,100*200,4)
data=as.data.frame(data)
colnames(data)=c("x","y","z1","z2")
data$x=rep(seq(0,1,length.out=100),each=200)
data$y=rep(seq(0,2,length.out=200),100)
data$z1=exp(data$x+data$y)
data$z2=exp(data$x+data$y)
max=tail(data$z1,n=1)

vec1=c(1,1)/sqrt(2)*1/2
vec2=c(1,-1)/sqrt(2)*1/2
vec3=c(1,1)/sqrt(2) *1/2
vec4=c(1,-1)/sqrt(2) * 1/2

# p=ggplot(data,aes(x=x,y=y,fill=z1))+
#   coord_fixed()+
#   geom_raster()+
#  geom_contour(breaks=exp(0.3+0.924057),aes(z=z1),color="black",size=2)+
#   scale_fill_gradientn(guide =guide_colorbar(title=expression(f(theta[1],theta[2]))),colors=c("dodgerblue1","lightgreen","gold","firebrick1"),values=c(0,0.05,0.25,1))+
#   annotate("point",x=0.30,y=0.924057,size=5)+
#   annotate("segment", x = 0.30, xend = (0.3+vec1[1]), y = 0.924057, yend = (0.924057+vec1[2]), colour="black", size=2, arrow=arrow())+
#   annotate("segment", x = 0.30, xend = (0.3+vec2[1]), y = 0.924057, yend = (0.924057+vec2[2]), colour="black", size=2, arrow=arrow())+
#     theme_bw(base_size = 12) + labs(x=expression(theta[1]), y=expression(theta[2]))+
#   theme(axis.title.y=element_text(vjust=1.5))+
#   theme(legend.key = element_blank(),legend.background=element_blank())+
#   theme(panel.border = element_blank(),
#         panel.background = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.line.x = element_line(color = "black",size=0.4),
#         axis.line.y = element_line(color = "black",size=0.4))+
#   scale_x_continuous(limits=c(0,1),expand=(c(0,0)))+
#   scale_y_continuous(limits=c(0,2),expand=c(0,0))

p=ggplot(data,aes(x=x,y=y,fill=z1))+
  coord_fixed()+
  geom_raster()+
  geom_contour(breaks=exp(0.3+1),aes(z=z1),color="grey50",size=2)+
  scale_fill_gradientn(guide =guide_colorbar(title=expression(f(theta[1],theta[2]))),colors=c("dodgerblue1","lightgreen","gold","firebrick1"),values=c(0,0.05,0.25,1),limits=c(1,max))+
  annotate("point",x=0.30,y=1,size=5)+
  annotate("segment", x = 0.30, xend = (0.3+vec1[1]), y = 1, yend = (1+vec1[2]), colour="black", size=2, arrow=arrow())+
  annotate("segment", x = 0.30, xend = (0.3+vec2[1]), y = 1, yend = (1+vec2[2]), colour="black", size=2, arrow=arrow())+
  theme_bw(base_size = 12) + labs(x=expression(theta[1]), y=expression(theta[2]))+
  theme(axis.title.y=element_text(vjust=1.5))+
  theme(legend.key = element_blank(),legend.background=element_blank())+
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color = "black",size=0.4),
        axis.line.y = element_line(color = "black",size=0.4))+
  scale_x_continuous(limits=c(0,1),expand=(c(0,0)))+
  scale_y_continuous(limits=c(0,2),expand=c(0,0))


p
ggsave(file="Exp1a.pdf", width = 3.25, height = 3.25)


# p=ggplot(data,aes(x=x,y=y,fill=z2))+
#   coord_fixed()+
#   geom_raster()+
#   geom_contour(bins=6,aes(z=log(z1)),color="black",size=1)+
#   scale_fill_gradientn(guide =guide_colorbar(title=expression(g(theta[1],theta[2]))),colors=c("dodgerblue1","lightgreen","gold","firebrick1"),values=c(0,0.05,0.25,1))+
#   #annotate("point",x=0.30,y=0.924057,size=5)+
#   #annotate("segment", x = 0.30, xend = (0.3+vec3[1]), y = 0.924057, yend = (0.924057+vec3[2]), colour="black", size=2, arrow=arrow())+
#   #annotate("segment", x = 0.30, xend = (0.3+vec4[1]), y = 0.924057, yend = (0.924057+vec4[2]), colour="black", size=2, arrow=arrow())+
#   theme_bw(base_size = 12) + labs(x=expression(theta[1]), y=expression(theta[2]))+
#   theme(axis.title.y=element_text(vjust=1.5))+
#   theme(legend.key = element_blank(),legend.background=element_blank())+
#   theme(panel.border = element_blank(),
#         panel.background = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.line.x = element_line(color = "black",size=0.4),
#         axis.line.y = element_line(color = "black",size=0.4))+
#   scale_x_continuous(limits=c(0,1),expand=(c(0,0)))+
#   scale_y_continuous(limits=c(0,2),expand=c(0,0))
# 
# p
# 

set.seed(0)
sample=as.data.frame(randomLHS(10^4,2))
sample[,2]=2*sample[,2]
sample[,3]=exp(sample[,1]+sample[,2])
sample[,4]=exp(sample[,1]+sample[,2])
colnames(sample)=c("t1","t2","y","x")


p=ggplot(sample,aes(x=x,y=y,color=y))+
  geom_point()+
  theme_bw(base_size = 12) + labs(x=expression(g(theta[1]+theta[2])), y=expression(f(theta[1]+theta[2])))+
  scale_color_gradientn(guide =FALSE,colors=c("dodgerblue1","lightgreen","gold","firebrick1"),values=c(0,0.05,0.25,1),limits=c(1,max))+
  theme(axis.title.y=element_text(vjust=1.5))+
  theme(legend.key = element_blank(),legend.background=element_blank())+
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color = "black",size=0.4),
        axis.line.y = element_line(color = "black",size=0.4))+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))

p
ggsave(file="Exp1b.pdf", width = 3.25, height = 3.25)



####################################################################

data=matrix(NA,100*200,4)
data=as.data.frame(data)
colnames(data)=c("x","y","z1","z2")
data$x=rep(seq(0,1,length.out=100),each=200)
data$y=rep(seq(0,2,length.out=200),100)
data$z1=exp(data$x*data$y)
data$z2=exp((0.183374* data$x + 0.386973* data$y)^2* 2.11029)
max=tail(data$z1,n=1)

# vec1=c(0.924057,0.3)/sqrt(2) 
# vec2=c(0.3,-0.924057)/sqrt(2) 
vec1=c(1,0.3)/sqrt(1^2+0.3^2) * 1/2
vec2=c(0.3,-1)/sqrt(1^2+0.3^2)* 1/2
vec3=c(0.428222,0.903673)/sqrt(2) * 1/2
vec4=c(0.903673,-0.428222)/sqrt(2) * 1/2

# p=ggplot(data,aes(x=x,y=y,fill=z1))+
#   coord_fixed()+
#   geom_raster()+
#   geom_contour(breaks=1+(max-1)*0.05,aes(z=z1),color="black",size=2)+
#   scale_fill_gradientn(guide =guide_colorbar(title=expression(f(theta[1],theta[2]))),colors=c("dodgerblue1","lightgreen","gold","firebrick1"),values=c(0,0.05,0.25,1))+
#   annotate("point",x=0.30,y=0.924057,size=5)+
#   annotate("segment", x = 0.30, xend = (0.3+vec1[1]), y = 0.924057, yend = (0.924057+vec1[2]), colour="black", size=2, arrow=arrow())+
#   annotate("segment", x = 0.30, xend = (0.3+vec2[1]), y = 0.924057, yend = (0.924057+vec2[2]), colour="black", size=2, arrow=arrow())+
#   theme_bw(base_size = 12) + labs(x=expression(theta[1]), y=expression(theta[2]))+
#   theme(axis.title.y=element_text(vjust=1.5))+
#   theme(legend.key = element_blank(),legend.background=element_blank())+
#   theme(panel.border = element_blank(),
#         panel.background = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.line.x = element_line(color = "black",size=0.4),
#         axis.line.y = element_line(color = "black",size=0.4))+
#   scale_x_continuous(limits=c(0,1),expand=(c(0,0)))+
#   scale_y_continuous(limits=c(0,2),expand=c(0,0))

p=ggplot(data,aes(x=x,y=y,fill=z1))+
  coord_fixed()+
  geom_raster()+
  geom_contour(breaks=exp(0.3),aes(z=z1),color="grey50",size=2)+
  scale_fill_gradientn(guide =guide_colorbar(title=expression(f(theta[1],theta[2]))),colors=c("dodgerblue1","lightgreen","gold","firebrick1"),values=c(0,0.05,0.25,1),limits=c(1,max))+
  annotate("point",x=0.30,y=1,size=5)+
  annotate("segment", x = 0.30, xend = (0.3+vec1[1]), y = 1, yend = (1+vec1[2]), colour="black", size=2, arrow=arrow())+
  annotate("segment", x = 0.30, xend = (0.3+vec2[1]), y = 1, yend = (1+vec2[2]), colour="black", size=2, arrow=arrow())+
  theme_bw(base_size = 12) + labs(x=expression(theta[1]), y=expression(theta[2]))+
  theme(axis.title.y=element_text(vjust=1.5))+
  theme(legend.key = element_blank(),legend.background=element_blank())+
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color = "black",size=0.4),
        axis.line.y = element_line(color = "black",size=0.4))+
  scale_x_continuous(limits=c(0,1),expand=(c(0,0)))+
  scale_y_continuous(limits=c(0,2),expand=c(0,0))



p
ggsave(file="Exp2a.pdf", width = 3.25, height = 3.25)



data_profile=matrix(NA,100,3)
data_profile=as.data.frame(data_profile)
colnames(data_profile)=c("x","y","z")
data_profile$x=seq(0.01,1,length.out=100)
data_profile$y=0.3/seq(0.01,1,length.out=100)
data_profile$z=exp(0.3)


p=ggplot(data_profile,aes(x=x,y=y,color=z))+
  geom_line(size=2)+
  annotate("point",x=0.30,y=1,size=5)+
  scale_color_gradientn(guide =FALSE,colors=c("dodgerblue1","lightgreen","gold","firebrick1"),values=c(0,0.05,0.25,1),limits=c(1,max))+
  theme_bw(base_size = 12) + labs(x=expression(paste(theta[1]^"*",", log-scale")), y=expression(paste("argmin c",(theta[2]),", log-scale")))+
  theme(axis.title.y=element_text(vjust=1.5))+
  theme(legend.key = element_blank(),legend.background=element_blank())+
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color = "black",size=0.4),
        axis.line.y = element_line(color = "black",size=0.4))+
  scale_x_log10(limits=c(0.1,1),expand=(c(0,0)))+
  scale_y_log10(limits=c(0.1,2),expand=c(0,0))

p
ggsave(file="Exp2b.pdf", width = 3.25, height = 3.25)

set.seed(0)
sample=as.data.frame(randomLHS(10^4,2))
sample[,2]=2*sample[,2]
sample[,3]=exp(sample[,1]*sample[,2])
sample[,4]=0.183374* sample[,1] + 0.386973* sample[,2]
colnames(sample)=c("t1","t2","y","x")
#plot(exp((0.183374* sample[,1] + 0.386973* sample[,2])^2*2.11029),exp(sample[,1]*sample[,2]),pch=16,cex=0.1)


p=ggplot(sample,aes(x=x,y=y,color=y))+
  geom_point()+
  theme_bw(base_size = 12) + labs(x=expression(paste(Q[a]^"T",theta)), y=expression(f(theta[1],theta[2])))+
  scale_color_gradientn(guide =FALSE,colors=c("dodgerblue1","lightgreen","gold","firebrick1"),values=c(0,0.05,0.25,1),limits=c(1,max))+
  theme(axis.title.y=element_text(vjust=1.5))+
  theme(legend.key = element_blank(),legend.background=element_blank())+
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color = "black",size=0.4),
        axis.line.y = element_line(color = "black",size=0.4))+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))

p
ggsave(file="Exp2c.pdf", width = 3.25, height = 3.25)



set.seed(0)
sample=as.data.frame(randomLHS(10^4,2))
sample[,2]=2*sample[,2]
sample[,3]=exp(sample[,1]*sample[,2])
sample[,4]=exp((0.183374* sample[,1] + 0.386973* sample[,2])^2* 2.11029)
colnames(sample)=c("t1","t2","y","x")


p=ggplot(sample,aes(x=x,y=y,color=y))+
  geom_point()+
  theme_bw(base_size =12) + labs(x=expression(g(theta[1]+theta[2])), y=expression(f(theta[1]+theta[2])))+
  scale_color_gradientn(guide =FALSE,colors=c("dodgerblue1","lightgreen","gold","firebrick1"),values=c(0,0.05,0.25,1),limits=c(1,max))+
  theme(axis.title.y=element_text(vjust=1.5))+
  theme(legend.key = element_blank(),legend.background=element_blank())+
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color = "black",size=0.4),
        axis.line.y = element_line(color = "black",size=0.4))+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))

p
ggsave(file="Exp2d.pdf", width = 3.25, height = 3.25)


# 
# 
# 
# 
# 
# 
# #################
# 
# data=matrix(NA,100*100*100,5)
# data=as.data.frame(data)
# colnames(data)=c("x","y","z1","z2")
# # data$a=rep(seq(0,1,length.out=100),each=100*100)
# # data$b=rep(rep(seq(0,1,length.out=100),each=100),100)
# # data$c=rep(seq(0,1,length.out=100),100*100)
# data$x=rep(seq(0,1,length.out=100),each=200)
# data$y=rep(seq(0,2,length.out=200),100)
# data$z1=data$x+log(1+data$y)
# data$z2=data$x-data$y
# #max=max(data$z1)
# 
# #a=0.04, b=0.73, c=0.11
# chi = as.matrix(c(1,log(1+0.31),0.04/(1+0.31)))
# Fmat= matrix(chi%o%chi,nrow=3)
# eigendecomp=eigen(Fmat)
# eigenvalues=eigendecomp$values
# eigenvectors=eigendecomp$vectors
# 
# vec1=c(1,1)/2
# vec2=c(1,-1)/2
# vec3=c(1,1)/2
# vec4=c(1,-1)/2
# 
# p=ggplot(data,aes(x=x,y=y,fill=z1))+
#   coord_fixed()+
#   geom_raster()+
#   geom_contour(breaks=exp(0.3+0.924057),aes(z=z1),color="black",size=2)+
#   scale_fill_gradientn(guide =guide_colorbar(title=expression(f(theta[1],theta[2]))),colors=c("dodgerblue1","lightgreen","gold","firebrick1"),values=c(0,0.05,0.25,1))+
#   annotate("point",x=0.30,y=0.924057,size=5)+
#   annotate("segment", x = 0.30, xend = (0.3+vec1[1]), y = 0.924057, yend = (0.924057+vec1[2]), colour="black", size=2, arrow=arrow())+
#   annotate("segment", x = 0.30, xend = (0.3+vec2[1]), y = 0.924057, yend = (0.924057+vec2[2]), colour="black", size=2, arrow=arrow())+
#   theme_bw(base_size = 12) + labs(x=expression(theta[1]), y=expression(theta[2]))+
#   theme(axis.title.y=element_text(vjust=1.5))+
#   theme(legend.key = element_blank(),legend.background=element_blank())+
#   theme(panel.border = element_blank(),
#         panel.background = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.line.x = element_line(color = "black",size=0.4),
#         axis.line.y = element_line(color = "black",size=0.4))+
#   scale_x_continuous(limits=c(0,1),expand=(c(0,0)))+
#   scale_y_continuous(limits=c(0,2),expand=c(0,0))
# 
# 
# p
# 
# p=ggplot(data,aes(x=x,y=y,fill=z2))+
#   coord_fixed()+
#   geom_raster()+
#   scale_fill_gradientn(guide =guide_colorbar(title=expression(g(theta[1],theta[2]))),colors=c("dodgerblue1","lightgreen","gold","firebrick1"),values=c(0,0.05,0.25,1))+
#   annotate("point",x=0.30,y=0.924057,size=5)+
#   annotate("segment", x = 0.30, xend = (0.3+vec3[1]), y = 0.924057, yend = (0.924057+vec3[2]), colour="black", size=2, arrow=arrow())+
#   annotate("segment", x = 0.30, xend = (0.3+vec4[1]), y = 0.924057, yend = (0.924057+vec4[2]), colour="black", size=2, arrow=arrow())+
#   theme_bw(base_size = 12) + labs(x=expression(theta[1]), y=expression(theta[2]))+
#   theme(axis.title.y=element_text(vjust=1.5))+
#   theme(legend.key = element_blank(),legend.background=element_blank())+
#   theme(panel.border = element_blank(),
#         panel.background = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.line.x = element_line(color = "black",size=0.4),
#         axis.line.y = element_line(color = "black",size=0.4))+
#   scale_x_continuous(limits=c(0,1),expand=(c(0,0)))+
#   scale_y_continuous(limits=c(0,2),expand=c(0,0))
# 
# p
# 
