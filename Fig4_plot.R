library(ggplot2)
library(openxlsx)
library(reshape)
library(ggalt)
library(ggtext)
library(scales)
library(dplyr)
element_textbox_highlight <- function(..., hi.labels = NULL, hi.fill = NULL,
                                      hi.col = NULL, hi.box.col = NULL, hi.family = NULL) {
  structure(
    c(element_textbox(...),
      list(hi.labels = hi.labels, hi.fill = hi.fill, hi.col = hi.col, hi.box.col = hi.box.col, hi.family = hi.family)
    ),
    class = c("element_textbox_highlight", "element_textbox", "element_text", "element")
  )
}
detach_package <- function(pkg, character.only = FALSE)
{
  if(!character.only)
  {
    pkg <- deparse(substitute(pkg))
  }
  search_item <- paste("package", pkg, sep = ":")
  while(search_item %in% search())
  {
    detach(search_item, unload = TRUE, character.only = TRUE)
  }
}
setwd("E:/zeroCOVID_NPI/Rt_0709/codeV0803")
dir = "SEIR_simulation_P4_1129/"
indir = "1025dataset"
dat = read.csv(paste0(indir,"/Rt&smooth_NPI.csv"),stringsAsFactors = F)
ctlist = list(c("210800","220200","450600","210600","340400"),
              c("140100","361100","131000","220100","210100"),
              c("370200","441900","110000","310000","320500"))
cn = c("SC","MC","LC")
A = do.call(rbind,lapply(seq(1,length(ctlist)),function(i){
  ct = ctlist[[i]]
  call = do.call(rbind,lapply(ct,function(ci){
    c = list.files(dir,pattern=ci)
    f = read.xlsx(paste0(dir,c), sheet = 1)
    f = f[,c(1:12)]
    f = do.call(rbind,lapply(split(f,f$R0),function(r){
      r = do.call(rbind,lapply(split(r,r$Latent),function(l){
        l = do.call(rbind,lapply(split(l,l$C),function(c){
          c = do.call(rbind,lapply(split(c,c$BI),function(b){
            b = do.call(rbind,lapply(split(b,b$Start.NPI),function(s){
              s$Mean = round(s$Mean,2)
              s$CI05 = round(s$CI05,2)
              s$CI95 = round(s$CI95,2)
              s$Days = seq(1,nrow(s))
              s$Date = as.Date(as.character(s$original_start),format="%Y%m%d")+s$Days
              if(length(which(s$Mean>=1))>1){
                s = s[min(which(s$Mean>=1)):max(which(s$Mean>=1)),]
                dur = nrow(s)
              }else{
                dur = 0
              }
              s1 = data.frame(cum_infection = sum(s$Mean),cum_infection_CI05 = sum(s$CI05),
                              cum_infection_CI95 = sum(s$CI95),Duration =  dur,
                              max_infection = max(s$Mean),max_infection_CI05 = max(s$CI05),
                              max_infection_CI95 = max(s$CI95),
                              Start.NPI = unique(s$Start.NPI),
                              SD = unique(s$BI),
                              CC = unique(s$C),
                              Latent = unique(s$Latent),
                              R0 = unique(s$R0))
              return(s1)
            }))
            return(b)
          }))
          return(c)
        }))
        l = subset(l, l$Duration<90&l$Duration>7&l$cum_infection<=10000)
        if(!is.null(l)){
          return(l)
        }else(return(NULL))
      }))
      return(r)
    }))
    if(nrow(f)>0){
      f$citycode = ci
      f$pop = unique(dat$pop[which(dat$citycode==ci)])
    }
    
    if(!is.null(f)){
      return(f)
    }else(return(NULL))
  }))
  
  call = do.call(rbind,lapply(split(call,call$R0),function(cs){
    cs = do.call(rbind,lapply(split(cs,cs$Latent),function(csc){
      csc = do.call(rbind,lapply(split(csc,csc$Start.NPI),function(cr){
        cr = do.call(rbind,lapply(split(cr,cr$CC),function(crl){
          sdk = sort(unique(crl$SD))
          crl$SDD = crl$SD
          for(k in seq(1,length(sdk))){
            if(length(which(crl$SDD==sdk[k]))>=3){
              crl=crl[which(crl$SD==sdk[k]),]
            }else{
              crl$SDD[which(crl$SDD==sdk[k])]=sdk[k+1]
            }
          }
          
          if(nrow(crl)>=3){
            
            c = do.call(rbind,lapply(split(crl,crl$SDD),function(crll){
              data.frame(cum_infection = mean(crll$cum_infection),cum_infection_CI05=mean(crll$cum_infection_CI05),
                         cum_infection_CI95=mean(crll$cum_infection_CI95),Duration=mean(crll$Duration),
                         max_infection = mean(crll$max_infection),max_infection_CI05 = mean(crll$max_infection_CI05),
                         max_infection_CI95 = mean(crll$max_infection_CI95),Start.NPI=unique(crll$Start.NPI),
                         SD = unique(crll$SDD),CC = unique(crll$CC),Latent = unique(crll$Latent),R0=unique(crll$R0),
                         pop=mean(crll$pop))
            } ))
              
            return(c)
          }else(return(NULL))
        }))
        
        return(cr)
      }))
      return(csc)
    }))
    return(cs)
  }))
  #call = subset(call, call$Duration<75&call$Duration>7&call$cum_infection<=50000)
  call$citytype = cn[i]
  return(call)
  print(paste(cn[i],"has been processed"))
}))
require(reshape2)
A$infection_pre = A$cum_infection#/(A$pop)
A$FM= 0.5
df <- A[,c("infection_pre","Duration","Start.NPI",
           "SD","CC","R0","Latent","citytype","pop")] %>% 
  group_by(infection_pre, Duration) %>%
  melt(id.vars=c("R0","infection_pre","Duration","Latent","citytype","pop"))%>%
  ungroup()

detach_package(reshape2)
require(dplyr)
colormap = list(
  "LC" = c("Start.NPI"="grey30","CC"="#8CBCBD",
           "SD" = "#00808c"),
  "MC" = c("Start.NPI"="grey30","CC"="#a52a2a",
           "SD" = "#6b2b2c"),
  "SC" =c("Start.NPI"="grey30","CC"="#a1a7c2",
          "SD" = "#455787")
)
my_plot_fun <- function(data){
  city = unique(data$citytype)
  data$value = as.numeric(data$value)
  data$value[which(data$variable=="Start.NPI")] = data$value[which(data$variable=="Start.NPI")]/21 
  ggplot(data, aes(variable, value, fill = variable)) + 
    geom_bar(stat = "identity", alpha = 0.85,width = 0.85,color= "transparent") + 
    scale_y_continuous(limits= c(0,1),expand = c(0,0))+
    coord_polar(theta = "y") + 
    scale_fill_manual(values = colormap[[city]])+
    scale_color_manual(values = colormap[[city]])+
    theme_void()+ guides(fill = "none",color = "none")
}
annotation_fun1 <- function(data, infection_pre,Duration,pop, plot_fun) {
  subplot = plot_fun(data)
  sub_grob <- annotation_custom(ggplotGrob(subplot), 
                                x = sqrt(infection_pre-sqrt(as.numeric(pop))*1), 
                                y = Duration-sqrt(as.numeric(pop))-1, 
                                xmax =  sqrt(infection_pre+sqrt(as.numeric(pop))*1),
                                ymax = Duration+sqrt(as.numeric(pop))+1)
}
annotation_fun2 <- function(data, infection_pre,Duration,pop, plot_fun) {
  subplot = plot_fun(data)
  sub_grob <- annotation_custom(ggplotGrob(subplot), 
                                x =  sqrt(infection_pre-sqrt(as.numeric(pop))*2.5), 
                                y = Duration-sqrt(as.numeric(pop))-2, 
                                xmax =   sqrt(infection_pre+sqrt(as.numeric(pop))*2.5),
                                ymax = Duration+sqrt(as.numeric(pop))+2)
}
annotation_fun3 <- function(data, infection_pre,Duration,pop, plot_fun) {
  subplot = plot_fun(data)
  sub_grob <- annotation_custom(ggplotGrob(subplot), 
                                x = sqrt(infection_pre-sqrt(as.numeric(pop))*3), 
                                y = Duration-sqrt(as.numeric(pop))-2, 
                                xmax =  sqrt(infection_pre+sqrt(as.numeric(pop))*3),
                                ymax = Duration+sqrt(as.numeric(pop))+2)
}
library(maps)
library(tidyverse)
outdir = "Figures"

for(r in split(df,df$R0)){
  
  for(l in split(r,r$Latent)){
    annotation_fun = ifelse(unique(l$R0=="3"),annotation_fun1,ifelse(unique(l$R0=="8"),annotation_fun2,annotation_fun3))
    k1 = split(l,l$citytype)[[1]]
    subgrobs1 <-k1 %>% 
      nest(-Duration,-infection_pre,-pop)  %>%
      pmap(annotation_fun,plot_fun = my_plot_fun)
    l1 = unique(l[,c("infection_pre","Duration")])
    g  =  ggplot(data =l1, aes(x=infection_pre,y=Duration))+
      geom_point(color = "transparent")+
      subgrobs1+
      theme_bw()+
      scale_y_continuous(limits = c(10,100))+
      geom_hline(aes(yintercept = 30),lty="dashed",color= "grey30",size = 0.5)+
      geom_hline(aes(yintercept = 60),lty="dashed",color= "grey50",size = 0.3)+
      geom_hline(aes(yintercept = 90),lty="dashed",color= "grey70",size = 0.2)+
      theme(legend.position ="top",
            legend.background = element_rect(
              fill = "transparent", # Ìî³äÉ«
              colour = "transparent", # ¿òÏßÉ«
              size = 0.5),
            legend.title=element_blank(),
            panel.grid = element_blank(),
            panel.border = element_rect(linetype = "dashed", color="grey70"),
            #axis.line=element_line(color="black",size=0.5),
            plot.title = element_text(color="black",hjust = 0,vjust=0,size =unit(8,"pt")),
            axis.title.x= element_text(color="black",size = unit(10,"pt")),
            axis.text.x = element_text(color="black",size = unit(8,"pt")),
            axis.title.y= element_text(color="black",size = unit(10,"pt")),
            axis.text.y = element_text(color="black",size = unit(8,"pt")),
            plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"))
    
    g = g + scale_x_sqrt(limits = c(round(min(r$infection_pre),0)-100,round(max(r$infection_pre),0)+100))
    
    if(length(unique(l$citytype))==2){
      k2 = split(l,l$citytype)[[2]]
      subgrobs2 <-k2 %>% 
        nest(-Duration,-infection_pre,-pop)  %>%
        pmap(annotation_fun,plot_fun = my_plot_fun)
      g = g+subgrobs2
    }else if(length(unique(l$citytype))==3){
      k2 = split(l,l$citytype)[[2]]
      subgrobs2 <-k2 %>% 
        nest(-Duration,-infection_pre,-pop)  %>%
        pmap(annotation_fun,plot_fun = my_plot_fun)
      k3 = split(l,l$citytype)[[3]]
      subgrobs3 <-k3 %>% 
        nest(-Duration,-infection_pre,-pop)  %>%
        pmap(annotation_fun,plot_fun = my_plot_fun)
      g = g+subgrobs2+subgrobs3
    }
    g = g + labs(x="Total number of the infections",y = "Duration(days)")
    ggsave(paste0(outdir,"/Fig4_R0",unique(l$R0),"_Latent",unique(l$Latent),".pdf"),g, width=7, height=6, units="cm", scale=1)
    
  }
}

write.csv(A,paste0(outdir,"/Fig4_data.csv"),row.names = F)

