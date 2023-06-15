library("rstan")
library("bayesplot")
library("ggplot2")
library("gridExtra")
###########Priori and Posterior plot##########
setwd("E:/zeroCOVID_NPI/Version0504")
set.seed(123)
####################alpha##################
u<-as.data.frame(rgamma(5000,0.1667,1)-log(1.05)/6.0)
colnames(u)<-"Priori"
plot_title <- expression(paste("Priori distributions of ", alpha, " ~ gamma(0.1667,1) - log(1.05)/6"))

N = length(u$Priori) 
mean = mean(u$Priori)
sd = sd(u$Priori)
se = sd/sqrt(N) 
dense = data.frame(density(u$Priori)[c('x','y')])
p1<-ggplot()+
  geom_area(data=dense,aes(x,y),
            colour="#606882",fill="#c8d5e5",size=0.5)+
  #mcmc_areas(u,#par=c("Priori","Posterior"),par=c("Priori"),prob = 0.95) + 
  labs(title = plot_title,x="",y = "Density")+
  geom_vline(aes(xintercept=0),
             color="#606882", linetype="dashed", size=0.5)+
  theme_classic()+
  theme(legend.position = "none",
        legend.background = element_rect(fill=NA),
        legend.key.size = unit(9, "pt"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title=element_blank(),
        plot.title = element_text(color="black",hjust = 0,vjust=0,size = unit(7, "pt")),
        axis.text.x= element_text(color="black",size = unit(9, "pt")),
        axis.title.x=element_blank(),
        plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"),
        panel.background=element_rect(fill = "transparent",colour = NA))
###############################beta############################
u<-as.data.frame(rnorm(5000,0,0.5))
colnames(u)<-"Priori"
plot_title <- expression(paste("Priori distributions of ", beta, " ~ normal(0,0.5) "))

N = length(u$Priori) 
mean = mean(u$Priori)
sd = sd(u$Priori)
se = sd/sqrt(N) 
dense = data.frame(density(u$Priori)[c('x','y')])

p2<-ggplot()+
  geom_area(data=dense,aes(x,y),
            colour="#606882",fill="#c8d5e5",size=0.5)+
  #mcmc_areas(u,#par=c("Priori","Posterior"),par=c("Priori"),prob = 0.95) + 
  labs(title = plot_title,x="",y = "Density")+
  geom_vline(aes(xintercept=0),
             color="#606882", linetype="dashed", size=0.5)+
  theme_classic()+
  theme(legend.position = "none",
        legend.background = element_rect(fill=NA),
        legend.key.size = unit(9, "pt"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title=element_blank(),
        axis.title.x=element_blank(),
        plot.title = element_text(color="black",hjust = 0,vjust=0,size = unit(7, "pt")),
        axis.text.x= element_text(color="black",size = unit(9, "pt")),
        plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"),
        panel.background=element_rect(fill = "transparent",colour = NA))

plot<-grid.arrange(arrangeGrob(grobs =list(p1,p2),ncol =2))
outpath = "Figures_Reviewed"
ggsave(paste0(outpath,"/Priori_NPIeffect.pdf"),plot,units="mm",width=150,height=60,device = cairo_pdf)
ggsave(paste0(outpath,"/Priori_NPIeffect.jpg"),plot,units="mm",width=150,height=60)


##########################Serival interval#####################
u<-as.data.frame(rgamma(5000, 5.80, 3.20))
colnames(u)<-"Priori"
plot_title <- expression(paste("Original SARS-CoV-2: Gamma(5.80,3.20)"))

N = length(u$Priori) 
mean = mean(u$Priori)
sd = sd(u$Priori)
se = sd/sqrt(N) 
dense = data.frame(density(u$Priori)[c('x','y')])

p1<-ggplot()+
  geom_area(data=dense,aes(x,y),
            colour="#606882",fill="#c8d5e5",size=0.5)+
  #mcmc_areas(u,#par=c("Priori","Posterior"),par=c("Priori"),prob = 0.95) + 
  labs(title = plot_title,x="",y = "Density")+
  geom_vline(aes(xintercept=mean),
             color="#606882", linetype="dashed", size=0.5)+
  theme_classic()+
  theme(legend.position = "none",
        legend.background = element_rect(fill=NA),
        legend.key.size = unit(9, "pt"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title=element_blank(),
        axis.title.x=element_blank(),
        plot.title = element_text(color="black",hjust = 0,vjust=0,size = unit(9, "pt")),
        axis.text.x= element_text(color="black",size = unit(9, "pt")),
        plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"),
        panel.background=element_rect(fill = "transparent",colour = NA))
##alpha##
u<-as.data.frame(rgamma(5000, 3.18, 4.36))
colnames(u)<-"Priori"
plot_title <- expression(paste("Alpha: Gamma(3.18,4.36)"))

N = length(u$Priori) 
mean = mean(u$Priori)
sd = sd(u$Priori)
se = sd/sqrt(N) 
dense = data.frame(density(u$Priori)[c('x','y')])

p2<-ggplot()+
  geom_area(data=dense,aes(x,y),
            colour="#606882",fill="#c8d5e5",size=0.5)+
  #mcmc_areas(u,#par=c("Priori","Posterior"),par=c("Priori"),prob = 0.95) + 
  labs(title = plot_title,x="",y = "Density")+
  geom_vline(aes(xintercept=mean),
             color="#606882", linetype="dashed", size=0.5)+
  theme_classic()+
  theme(legend.position = "none",
        legend.background = element_rect(fill=NA),
        legend.key.size = unit(9, "pt"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title=element_blank(),
        axis.title.x=element_blank(),
        plot.title = element_text(color="black",hjust = 0,vjust=0,size = unit(9, "pt")),
        axis.text.x= element_text(color="black",size = unit(9, "pt")),
        plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"),
        panel.background=element_rect(fill = "transparent",colour = NA))
##Delta##
u<-as.data.frame(rgamma(5000, 2.3, 3.4))
colnames(u)<-"Priori"
plot_title <- expression(paste("Delta: Gamma(2.30,3.40) "))

N = length(u$Priori) 
mean = mean(u$Priori)
sd = sd(u$Priori)
se = sd/sqrt(N) 
dense = data.frame(density(u$Priori)[c('x','y')])

p3<-ggplot()+
  geom_area(data=dense,aes(x,y),
            colour="#606882",fill="#c8d5e5",size=0.5)+
  #mcmc_areas(u,#par=c("Priori","Posterior"),par=c("Priori"),prob = 0.95) + 
  labs(title = plot_title,x="",y = "Density")+
  geom_vline(aes(xintercept=mean),
             color="#606882", linetype="dashed", size=0.5)+
  theme_classic()+
  theme(legend.position = "none",
        legend.background = element_rect(fill=NA),
        legend.key.size = unit(9, "pt"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title=element_blank(),
        axis.title.x=element_blank(),
        plot.title = element_text(color="black",hjust = 0,vjust=0,size = unit(9, "pt")),
        axis.text.x= element_text(color="black",size = unit(9, "pt")),
        plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"),
        panel.background=element_rect(fill = "transparent",colour = NA))
##Omicron##
u<-as.data.frame(rgamma(5000, 2.9, 1.6))
colnames(u)<-"Priori"
plot_title <- expression(paste("Omicron: Gamma(2.90,1.60) "))

N = length(u$Priori) 
mean = mean(u$Priori)
sd = sd(u$Priori)
se = sd/sqrt(N) 
dense = data.frame(density(u$Priori)[c('x','y')])

p4<-ggplot()+
  geom_area(data=dense,aes(x,y),
            colour="#606882",fill="#c8d5e5",size=0.5)+
  #mcmc_areas(u,#par=c("Priori","Posterior"),par=c("Priori"),prob = 0.95) + 
  labs(title = plot_title,x="",y = "Density")+
  geom_vline(aes(xintercept=mean),
             color="#606882", linetype="dashed", size=0.5)+
  theme_classic()+
  theme(legend.position = "none",
        legend.background = element_rect(fill=NA),
        legend.key.size = unit(9, "pt"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title=element_blank(),
        axis.title.x=element_blank(),
        plot.title = element_text(color="black",hjust = 0,vjust=0,size = unit(9, "pt")),
        axis.text.x= element_text(color="black",size = unit(9, "pt")),
        plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"),
        panel.background=element_rect(fill = "transparent",colour = NA))

plot<-grid.arrange(arrangeGrob(grobs =list(p1,p2,p3,p4),ncol =2))
outpath = "Figures_Reviewed"
ggsave(paste0(outpath,"/Priori_serial_interval.pdf"),plot,units="mm",width=150,height=120,device = cairo_pdf)
ggsave(paste0(outpath,"/Priori_serial_interval.jpg"),plot,units="mm",width=150,height=120)


######incubation period####
##########################Serival interval#####################
u<-as.data.frame(rlnorm(5000, meanlog = 1.78, sdlog = 0.52))
colnames(u)<-"Priori"
plot_title <- expression(paste("Incubation period of original SARS-CoV-2: lognormal(1.78,0.52) "))

N = length(u$Priori) 
mean = mean(u$Priori)
sd = sd(u$Priori)
se = sd/sqrt(N) 
dense = data.frame(density(u$Priori)[c('x','y')])

p1<-ggplot()+
  geom_area(data=dense,aes(x,y),
            colour="#606882",fill="#c8d5e5",size=0.5)+
  #mcmc_areas(u,#par=c("Priori","Posterior"),par=c("Priori"),prob = 0.95) + 
  labs(title = plot_title,x="",y = "Density")+
  geom_vline(aes(xintercept=mean),
             color="#606882", linetype="dashed", size=0.5)+
  theme_classic()+
  theme(legend.position = "none",
        legend.background = element_rect(fill=NA),
        legend.key.size = unit(9, "pt"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title=element_blank(),
        axis.title.x=element_blank(),
        plot.title = element_text(color="black",hjust = 0,vjust=0,size = unit(7, "pt")),
        axis.text.x= element_text(color="black",size = unit(9, "pt")),
        plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"),
        panel.background=element_rect(fill = "transparent",colour = NA))
##alpha##
u<-as.data.frame(rlnorm(5000, meanlog = 1.50, sdlog = 0.46))
colnames(u)<-"Priori"
plot_title <- expression(paste("Incubation period of variant Alpha: lognormal(1.50,0.46) "))

N = length(u$Priori) 
mean = mean(u$Priori)
sd = sd(u$Priori)
se = sd/sqrt(N) 
dense = data.frame(density(u$Priori)[c('x','y')])

p2<-ggplot()+
  geom_area(data=dense,aes(x,y),
            colour="#606882",fill="#c8d5e5",size=0.5)+
  #mcmc_areas(u,#par=c("Priori","Posterior"),par=c("Priori"),prob = 0.95) + 
  labs(title = plot_title,x="",y = "Density")+
  geom_vline(aes(xintercept=mean),
             color="#606882", linetype="dashed", size=0.5)+
  theme_classic()+
  theme(legend.position = "none",
        legend.background = element_rect(fill=NA),
        legend.key.size = unit(9, "pt"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title=element_blank(),
        axis.title.x=element_blank(),
        plot.title = element_text(color="black",hjust = 0,vjust=0,size = unit(7, "pt")),
        axis.text.x= element_text(color="black",size = unit(9, "pt")),
        plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"),
        panel.background=element_rect(fill = "transparent",colour = NA))
##Delta##
u<-as.data.frame(rlnorm(5000, meanlog = 1.25, sdlog = 0.34))
colnames(u)<-"Priori"
plot_title <- expression(paste("Incubation period of variant Delta: lognormal(1.25,0.34) "))

N = length(u$Priori) 
mean = mean(u$Priori)
sd = sd(u$Priori)
se = sd/sqrt(N) 
dense = data.frame(density(u$Priori)[c('x','y')])

p3<-ggplot()+
  geom_area(data=dense,aes(x,y),
            colour="#606882",fill="#c8d5e5",size=0.5)+
  #mcmc_areas(u,#par=c("Priori","Posterior"),par=c("Priori"),prob = 0.95) + 
  labs(title = plot_title,x="",y = "Density")+
  geom_vline(aes(xintercept=mean),
             color="#606882", linetype="dashed", size=0.5)+
  theme_classic()+
  theme(legend.position = "none",
        legend.background = element_rect(fill=NA),
        legend.key.size = unit(9, "pt"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title=element_blank(),
        axis.title.x=element_blank(),
        plot.title = element_text(color="black",hjust = 0,vjust=0,size = unit(7, "pt")),
        axis.text.x= element_text(color="black",size = unit(9, "pt")),
        plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"),
        panel.background=element_rect(fill = "transparent",colour = NA))
##Omicron##
u<-as.data.frame(rlnorm(5000, meanlog = 1.02, sdlog = 0.45))
colnames(u)<-"Priori"
plot_title <- expression(paste("Incubation period of variant Omicron: lognormal(1.02,0.45)"))

N = length(u$Priori) 
mean = mean(u$Priori)
sd = sd(u$Priori)
se = sd/sqrt(N) 
dense = data.frame(density(u$Priori)[c('x','y')])

p4<-ggplot()+
  geom_area(data=dense,aes(x,y),
            colour="#606882",fill="#c8d5e5",size=0.5)+
  #mcmc_areas(u,#par=c("Priori","Posterior"),par=c("Priori"),prob = 0.95) + 
  labs(title = plot_title,x="",y = "Density")+
  geom_vline(aes(xintercept=mean),
             color="#606882", linetype="dashed", size=0.5)+
  theme_classic()+
  theme(legend.position = "none",
        legend.background = element_rect(fill=NA),
        legend.key.size = unit(9, "pt"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title=element_blank(),
        axis.title.x=element_blank(),
        plot.title = element_text(color="black",hjust = 0,vjust=0,size = unit(7, "pt")),
        axis.text.x= element_text(color="black",size = unit(9, "pt")),
        plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"),
        panel.background=element_rect(fill = "transparent",colour = NA))

plot<-grid.arrange(arrangeGrob(grobs =list(p1,p2,p3,p4),ncol =2))
outpath = "Figures_Reviewed"
ggsave(paste0(outpath,"/Priori_Incubation period.pdf"),plot,units="mm",width=150,height=120,device = cairo_pdf)
ggsave(paste0(outpath,"/Priori_Incubation period.jpg"),plot,units="mm",width=150,height=120)


##report delay##
u<-as.data.frame(rlnorm(5000, meanlog = 0.82, sdlog = 0.84))
colnames(u)<-"Priori"
plot_title <- expression(paste("lognormal(0.82,0.84)"))

N = length(u$Priori) 
mean = mean(u$Priori)
sd = sd(u$Priori)
se = sd/sqrt(N) 
dense = data.frame(density(u$Priori)[c('x','y')])

p1<-ggplot()+
  geom_area(data=dense,aes(x,y),
            colour="#606882",fill="#c8d5e5",size=0.5)+
  #mcmc_areas(u,#par=c("Priori","Posterior"),par=c("Priori"),prob = 0.95) + 
  labs(title = plot_title,x="",y = "Density")+
  geom_vline(aes(xintercept=mean),
             color="#606882", linetype="dashed", size=0.5)+
  theme_classic()+
  theme(legend.position = "none",
        legend.background = element_rect(fill=NA),
        legend.key.size = unit(9, "pt"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title=element_blank(),
        axis.title.x=element_blank(),
        plot.title = element_text(color="black",hjust = 0,vjust=0,size = unit(7, "pt")),
        axis.text.x= element_text(color="black",size = unit(9, "pt")),
        plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"),
        panel.background=element_rect(fill = "transparent",colour = NA))
###onset-to-report
u<-as.data.frame(rbinom(5000, meanlog = 1, sdlog = 0.3))
colnames(u)<-"Priori"
plot_title <- expression(paste("binomial(1,0.3)"))

N = length(u$Priori) 
mean = mean(u$Priori)
sd = sd(u$Priori)
se = sd/sqrt(N) 
dense = data.frame(density(u$Priori)[c('x','y')])

p2<-ggplot()+
  geom_area(data=dense,aes(x,y),
            colour="#606882",fill="#c8d5e5",size=0.5)+
  #mcmc_areas(u,#par=c("Priori","Posterior"),par=c("Priori"),prob = 0.95) + 
  labs(title = plot_title,x="",y = "Density")+
  geom_vline(aes(xintercept=mean),
             color="#606882", linetype="dashed", size=0.5)+
  theme_classic()+
  theme(legend.position = "none",
        legend.background = element_rect(fill=NA),
        legend.key.size = unit(9, "pt"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title=element_blank(),
        axis.title.x=element_blank(),
        plot.title = element_text(color="black",hjust = 0,vjust=0,size = unit(7, "pt")),
        axis.text.x= element_text(color="black",size = unit(9, "pt")),
        plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"),
        panel.background=element_rect(fill = "transparent",colour = NA))

plot<-grid.arrange(arrangeGrob(grobs =list(p1,p2),ncol =2))
outpath = "Figures_Reviewed"
ggsave(paste0(outpath,"/Priori_timelag_report.pdf"),plot,units="mm",width=150,height=60,device = cairo_pdf)
ggsave(paste0(outpath,"/Priori_timelag_report.jpg"),plot,units="mm",width=150,height=60)









