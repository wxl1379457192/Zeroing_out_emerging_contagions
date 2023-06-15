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
setwd("E:/zeroCOVID_NPI/Version0504")
dir = "SEIR_simulation_P3/"
ctlist = list(c("210800","220200","450600","210600","340400"),
              c("140100","361100","131000","220100","210100"),
           c("370200","441900","110000","310000","320500"))
cn = c("SC","MC","LC")
A = do.call(rbind,lapply(seq(1,length(ctlist)),function(i){
  ct = ctlist[[i]]
  call = do.call(rbind,lapply(ct,function(ci){
    c = list.files(dir,pattern=ci)
    f = read.xlsx(paste0(dir,c), sheet = 1)
    f = do.call(rbind,lapply(split(f,f$R0),function(r){
      r = do.call(rbind,lapply(split(r,r$Latent),function(l){
        
        l = l[which(l$Start.NPI>0),]
        l$dailymax = l$max/l$duration
        l = do.call(rbind,lapply(seq(1,nrow(l)/11),function(j){
          l1 = l[c(((j-1)*11+1):(j*11)),]
          if(max(l1$BI)>0){
            for(k in c(1,2,3,7,8,9,10,11,12,22)){
              l1[,k]= round(l1[,k],4)
              if(l1[which(l1$BI==0),k]>=1&!is.na(l1[which(l1$BI==0),k])){
                l1[,k] = (1-l1[,k]/l1[which(l1$BI==0),k])
              }else{
                l1[,k] = 1
              }
            }
            l1$var = "BI"
            l1$intensity = l1$BI
          }else if(max(l1$FI)>0){
            for(k in c(1,2,3,7,8,9,10,11,12,22)){
              l1[,k]= round(l1[,k],4)
              if(l1[which(l1$FI==0),k]>=1&!is.na(l1[which(l1$FI==0),k])){
                l1[,k] = (1-l1[,k]/l1[which(l1$FI==0),k])
              }else{
                l1[,k] = 1
              }
            }
            l1$var = "FI"
            l1$intensity = l1$FI
          }else if(max(l1$SI)>0){
            for(k in c(1,2,3,7,8,9,10,11,12,22)){
              l1[,k]= round(l1[,k],4)
              if(l1[which(l1$SI==0),k]>=1&!is.na(l1[which(l1$SI==0),k])){
                l1[,k] = (1-l1[,k]/l1[which(l1$SI==0),k])
              }else{
                l1[,k] = 1
              }
            }
            l1$var = "SI"
            l1$intensity = l1$SI
          }else{
            for(k in c(1,2,3,7,8,9,10,11,12,22)){
              l1[,k]= round(l1[,k],4)
              if(l1[which(l1$C==0),k]>=1&!is.na(l1[which(l1$C==0),k])){
                l1[,k] = (1-l1[,k]/l1[which(l1$C==0),k])
              }else{
                l1[,k] = 1
              }
            }
            l1$var="C"
            l1$intensity = l1$C
          }
          return(l1)
        }))
        return(l)
      }))
      return(r)
    }))
    return(f)
  }))
  for(k in c(1,2,3,7,8,9,10,11,12,22)){
    call[,k]= round(call[,k],5)
    call[which(call[,k]<0),k] = 0
  }
  call = call%>%group_by(R0,Start.NPI,Latent,var,intensity)%>%
    summarise(Mean = mean(Mean[!is.na(Mean)]),CI05=mean(CI05[!is.na(CI05)]),
              CI95=mean(CI95[!is.na(CI95)]),Sum =  mean(Sum[!is.na(Sum)]),
              Sum_CI05=mean(Sum_CI05[!is.na(Sum_CI05)]),
              Sum_CI95 =mean(Sum_CI95[!is.na(Sum_CI95)]),
              max =  mean(max[!is.na(max)]),max_CI05=mean(max_CI05[!is.na(max_CI05)]),
              max_CI95 =mean(max_CI95[!is.na(max_CI95)]),
              dailymax = mean(dailymax[!is.na(dailymax)]))
  call$citytype = cn[i]
  return(call)
}))
A1 = subset(A,A$var=="SI")
A1 = subset(A1,A1$R0=="3")
#indir = "1025dataset"
#dat = read.csv(paste0(indir,"/Rt&smooth_NPI.csv"),stringsAsFactors = F)

#f1 = subset(f1,f1$Latent==2)
f2 = A[,c("Mean","citytype","R0","Start.NPI","var","intensity","Latent")]
f2$Start.NPI=as.numeric(f2$Start.NPI)
#f1 = subset(f2,f2$var=="BI")
f2$intensity=as.numeric(f2$intensity)
#f2$dailymax[which(f2$dailymax<0)] = 0
f2$citytype = factor(f2$citytype,ordered=TRUE,levels = c("LC","MC","SC"))

f2$R0 = factor(f2$R0,ordered=TRUE,levels = c("3","8","13"))
f2$R0=as.numeric(f2$R0)
f2$citytype=as.numeric(f2$citytype)
#f2$R = paste("R0=",f2$R0,"Latent",f2$Latent)
#f2$dailymax[is.na(f2$dailymax)] = 0
g = ggplot(f2)+geom_raster(aes(x=Start.NPI,y = intensity,fill = Mean))+
  facet_grid(paste(R0,Latent)~paste(citytype,var))+
  scale_x_continuous(breaks = seq(5, 25, by = 10), expand = c(0,0))+
  scale_fill_gradientn(colours = c('#365a87','#13a69a',"grey80",
                                   '#fbd830','#c21a13'),
                       limits = c(0, 1))+
  scale_y_continuous(breaks = seq(0, 1, by = 0.3), expand = c(0,0))+
  labs(y="Intensity of NPI",x="Start day of NPI implementation")+
  theme(legend.position = "top",
        legend.background = element_rect(
          fill = "transparent", # Ìî³äÉ«
          colour = "transparent", # ¿òÏßÉ«
          size = 1.5),
        strip.text = element_textbox_highlight(
          size = 10, 
          fill = "transparent", box.color = "transparent", color = "gray10",
          halign = .5, linetype = 1, r = unit(0, "pt"), width = unit(1, "npc"),
          padding = margin(2, 0, 1, 0), margin = margin(3, 3, 3, 3)),
        strip.background.x = element_rect(fill = "white", colour = "transparent"),
        legend.title=element_blank(),
        panel.grid = element_blank(),
        axis.line=element_line(color="black",size=0.5),
        plot.title = element_text(color="black",hjust = 0,vjust=0,size =unit(9,"pt")),
        axis.title.x= element_text(color="black",hjust = 0.5,vjust=0,size =unit(9,"pt")),
        axis.text.x = element_text(color="black",size = unit(9,"pt")),
        axis.title.y= element_text(color="black",hjust = 0.5,vjust=0,size =unit(9,"pt")),
        axis.text.y = element_text(color="black",size = unit(9,"pt")),
        plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"))


outdir = "Figures_Reviewed"
ggsave(paste0(outdir,"/fig3.pdf"),g, width=18, height=17, units="cm", scale=1)

