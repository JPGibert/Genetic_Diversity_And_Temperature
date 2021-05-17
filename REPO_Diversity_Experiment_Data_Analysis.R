		####################################
		##### Singleton & Gibert, 2021 #####
		####################################
		
# install and load packages for data handling 
library("reshape2")
library("dplyr")
library("rstatix")
library("broom")
library("emmeans")
library("viridis")

##------------------------------------------------------------------------------------------
## 1) Data handling
# set working directory to the folder where the data file is 
setwd("~/Desktop/JP/Papers_in_progress/Alex_tetra/Data/Alex's_code")

# load data and visualize using function head 
tdata <- read.csv("Tetra_data_clones.csv")
head(tdata)

# Converting data from short format to long format
tdata_new <- melt(tdata,
id.vars=c("Temperature","Clone_number","Combination","Clone_A","Clone_B","Clone_C","Clone_D","Clone_E"))

# Calculating per capita growth rate (using as Ni=50in/ml and time=1 day)
tdata_new$growth <- log(tdata_new$value)-log(50)

# Redefine variables as needed (eg. Clone presence/absence as factor)
tdata_new <- tdata_new %>%
mutate(lgrowth=log10(growth+3),lclone=log10(Clone_number),Clone_A=as.factor(Clone_A),Clone_B=as.factor(Clone_B),Clone_C=as.factor(Clone_C),Clone_D=as.factor(Clone_D),Clone_E=as.factor(Clone_E))



##------------------------------------------------------------------------------------------
## 2) FIGURE 2a
# 1) Non-transformed data: Preliminary data analysis
# data analysis
mod <- lm(growth~Clone_number*Temperature,data=tdata_new)
summary(mod)

df_19 <- data.frame(Clone_number=seq(1,5,length.out=100),Temperature=rep(19,100))
df_22 <- data.frame(Clone_number=seq(1,5,length.out=100),Temperature=rep(22,100))
df_25 <- data.frame(Clone_number=seq(1,5,length.out=100),Temperature=rep(25,100))
df_28 <- data.frame(Clone_number=seq(1,5,length.out=100),Temperature=rep(28,100))
pred_19 <- t(t(predict(mod,df_19)))
pred_22 <- t(t(predict(mod,df_22)))
pred_25 <- t(t(predict(mod,df_25)))
pred_28 <- t(t(predict(mod,df_28)))

plot(growth~jitter(Clone_number),data=tdata_new,pch=16,col=c("blue","sky blue","orange","darkorange2")[as.factor(tdata_new$Temperature)])
lines(pred_19~df_19$Clone_number,col="blue", lwd=2)
lines(pred_22~df_22$Clone_number,col="sky blue", lwd=2)
lines(pred_25~df_25$Clone_number,col="orange", lwd=2)
lines(pred_28~df_28$Clone_number,col="darkorange2", lwd=2)


# 2) Log-transformed data: Preliminary data analysis
#plot(lgrowth~jitter(lclone),data=tdata_new,pch=16,col=c("blue","sky blue","orange","red")[as.factor(tdata_new$Temperature)],axes=FALSE)
#box(lwd=2)
#axis(1,at=)

# data analysis
modl <- lm(lgrowth~lclone*Temperature,data=tdata_new)
summary(modl)

df_19 <- data.frame(lclone=seq(0,0.7,length.out=100),Temperature=rep(19,100))
df_22 <- data.frame(lclone=seq(0,0.7,length.out=100),Temperature=rep(22,100))
df_25 <- data.frame(lclone=seq(0,0.7,length.out=100),Temperature=rep(25,100))
df_28 <- data.frame(lclone=seq(0,0.7,length.out=100),Temperature=rep(28,100))
pred_19 <- t(t(predict(modl,df_19)))
pred_22 <- t(t(predict(modl,df_22)))
pred_25 <- t(t(predict(modl,df_25)))
pred_28 <- t(t(predict(modl,df_28)))

# Figure
#pdf("~/Documents/Protist_Lab/fig_2.pdf",height=6,width=6)
plot(lgrowth~jitter(lclone),data=tdata_new,pch=16,col=c("blue","sky blue","orange","darkorange2")[as.factor(tdata_new$Temperature)],axes=FALSE,xlab="",ylab="",xlim=c(log10(0.8),log10(5.5)),ylim=c(log10(1.6),log10(6.5)))
box(lwd=2,bty="l")
axis(1,at=log10(c(1,2,3,4,5)),labels=c("1","2","3","4","5"),tck=0.02,padj=-1.25)
mtext("Genetic Diversity (number of clones)",1,line=1.8,cex=1.1)
axis(2,at=log10(c(2,3,4,5,6)),labels=c("-1","0","1","2","3"),tck=0.02,las=TRUE,hadj=-0.5)
mtext("Intrinsic Growth Rate (r)",2,line=1.7,cex=1.1)
lines(pred_19~df_19$lclone,col="blue", lwd=2)
lines(pred_22~df_22$lclone,col="sky blue", lwd=2)
lines(pred_25~df_25$lclone,col="orange", lwd=2)
lines(pred_28~df_28$lclone,col="darkorange2", lwd=2)
#dev.off()


##------------------------------------------------------------------------------------------
## 2) TABLE 1

## We select the monoclonal populations to understand the mechanisms through which these effects occur. Notice we're not using log values here.
tdata_single_cl <- tdata_new %>%
	group_by(Clone_number) %>%
	filter(Clone_number==1) %>%
			mutate(Combination=as.factor(Combination), Temperature=as.factor(Temperature))

# GxE
tdata_single_cl %>%
  anova_test(
    growth ~ Combination*Temperature
  )

# E
tdata_single_cl %>%
  group_by(Combination) %>%
  anova_test(growth ~ Temperature)

# G
tdata_single_cl %>%
  group_by(Temperature) %>%
  anova_test(growth ~ Combination)


## Post-hoc differences mong clones across temperatures.
pwc <- tdata_single_cl %>% 
  group_by(Temperature) %>%
  emmeans_test(
    growth ~ Combination,
    p.adjust.method = "bonferroni"
    )
pwc %>% print(n=120)


##------------------------------------------------------------------------------------------
## 2) FIGURE 2b

## Prep data
tdata_new_2 <- tdata_new %>%
	group_by(Temperature,Combination) %>%
	summarize(mean_growth=mean(growth, na.rm=TRUE),se=sd(growth, na.rm=TRUE)/sqrt(4))

## Offset the means in the horizontal axis for better visualization
Aplot <- data.frame(mean_growth=tdata_new_2[which(tdata_new_2$Combination=="A"),]$mean_growth,Temperature=tdata_new_2[which(tdata_new_2$Combination=="A"),]$Temperature-0.5,Combination=rep("A",4),se=tdata_new_2[which(tdata_new_2$Combination=="A"),]$se)
Bplot <- data.frame(mean_growth=tdata_new_2[which(tdata_new_2$Combination=="B"),]$mean_growth,Temperature=tdata_new_2[which(tdata_new_2$Combination=="B"),]$Temperature-0.25,Combination=rep("B",4),se=tdata_new_2[which(tdata_new_2$Combination=="B"),]$se)
Cplot <- data.frame(mean_growth=tdata_new_2[which(tdata_new_2$Combination=="C"),]$mean_growth,Temperature=tdata_new_2[which(tdata_new_2$Combination=="C"),]$Temperature,Combination=rep("C",4),se=tdata_new_2[which(tdata_new_2$Combination=="C"),]$se)
Dplot <- data.frame(mean_growth=tdata_new_2[which(tdata_new_2$Combination=="D"),]$mean_growth,Temperature=tdata_new_2[which(tdata_new_2$Combination=="D"),]$Temperature+0.25,Combination=rep("D",4),se=tdata_new_2[which(tdata_new_2$Combination=="D"),]$se)
Eplot <- data.frame(mean_growth=tdata_new_2[which(tdata_new_2$Combination=="E"),]$mean_growth,Temperature=tdata_new_2[which(tdata_new_2$Combination=="E"),]$Temperature+0.5,Combination=rep("E",4),se=tdata_new_2[which(tdata_new_2$Combination=="E"),]$se)

finalplot <- rbind(Aplot,Bplot,Cplot,Dplot,Eplot)

## Figure
plot(mean_growth~Temperature,data=finalplot[which(finalplot$Combination=="A"),],pch=16,type="b",xlim=c(18,29),ylim=c(-1,3),col=viridis(72)[71],axes=FALSE,xlab="",ylab="", lwd=2)
box(lwd=1.5,bty="l")
points(mean_growth~Temperature,data=finalplot[which(finalplot$Combination=="B"),],pch=16,type="b",col=viridis(72)[52], lwd=1.5)
points(mean_growth~Temperature,data=finalplot[which(finalplot$Combination=="C"),],pch=16,type="b",col=viridis(72)[32], lwd=1.5)
points(mean_growth~Temperature,data=finalplot[which(finalplot$Combination=="D"),],pch=16,type="b",col=viridis(72)[15], lwd=1.5)
points(mean_growth~Temperature,data=finalplot[which(finalplot$Combination=="E"),],pch=16,type="b",col=viridis(72)[64], lwd=1.5)

axis(1,(c(19,22,25,28)),labels=c("19","22","25","28"),tck=0.02,padj=-1)
mtext("Temperature",1,line=1.7,cex=1.1)
mtext("Intrinsic Growth Rate (r)",2,line=2,cex=1.1)
axis(2,c(-2,-1,0,1,2,3),tck=0.02,hadj=0.5,las=TRUE)

arrows(x0=finalplot$Temperature, 
y0=finalplot$mean_growth-finalplot$se, 
x1=finalplot$Temperature, 
y1=finalplot$mean_growth+finalplot$se, 
code=3, angle=90, length=0.025, col=c(viridis(72)[71],viridis(72)[52],viridis(72)[32],viridis(72)[15],viridis(72)[64])[as.factor(finalplot$Combination)], lwd=1)



##THE END