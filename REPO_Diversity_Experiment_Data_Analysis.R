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
library("nls.multstart")

##------------------------------------------------------------------------------------------
## 1) Data handling
# set working directory to the folder where the data file is 
setwd("~/Desktop/JP/Papers_in_review_submitted/Alex_tetra/Data/Alex's_code")

# load data and visualize using function head 
tdata <- read.csv("Tetra_data_clones.csv")
tpcdata <- read.csv("Tetra_TPC.csv")
head(tdata)
head(tpcdata)
# Converting data from short format to long format
tdata_new <- melt(tdata,
id.vars=c("Temperature","Clone_number","Combination","Clone_A","Clone_B","Clone_C","Clone_D","Clone_E"))

# Calculating per capita growth rate (using as Ni=50ind/ml and time=1 day)
tdata_new$growth <- log(tdata_new$value)-log(50)

# Redefine variables as needed (eg. Clone presence/absence as factor)
tdata_new <- tdata_new %>%
mutate(lgrowth=log10(growth+3),lclone=log10(Clone_number),Clone_A=as.factor(Clone_A),Clone_B=as.factor(Clone_B),Clone_C=as.factor(Clone_C),Clone_D=as.factor(Clone_D),Clone_E=as.factor(Clone_E))

# Calculating r for TPC data
tpcdata <- tpcdata %>%	
	mutate(initial= rep(15,dim(tpcdata)[1]),r=log10((Ind.ml+3)/(initial+3)))


##------------------------------------------------------------------------------------------
## 1) Fig 1c 
# define the Sharpe-Schoolfield equation
schoolfield_high <- function(lnc, E, Eh, Th, temp, Tc) {
Tc <- 273.15 + Tc
k <- 8.62e-5
boltzmann.term <- lnc + log(exp(E/k*(1/Tc - 1/temp)))
inactivation.term <- log(1/(1 + exp(Eh/k*(1/Th - 1/temp))))
return(boltzmann.term + inactivation.term)
}
# Fit
fits <- nls_multstart(r ~ schoolfield_high(lnc, E, Eh, Th, temp = Temp+273, Tc = 20),
data = tpcdata,
iter = 500,
start_lower = c(lnc=-10, E=0.1, Eh=0.5, Th=285),
start_upper = c(lnc=10, E=2, Eh=5, Th=330),
lower = c(lnc=-10, E=0, Eh=0, Th=0),
supp_errors = 'Y')

## Plot
pdf("~/Desktop/JP/Papers_in_review_submitted/Alex_tetra/Manuscript/Figures/Fig_1/fig_1c.pdf",height=5.5,width=12)
## Figure
par(mfrow=c(1,2))
plot(schoolfield_high(coef(fits)[1], coef(fits)[2], coef(fits)[3], coef(fits)[4], temp = Temp+273, Tc = 20) ~ Temp, data=tpcdata,type='b', ylim=c(-1,2.5),xlim=c(10,40), col="white", axes=FALSE, ylab="", xlab="")
schoolfield_high(coef(fits)[1], coef(fits)[2], coef(fits)[3], coef(fits)[4], temp = c(19,22,25,28)+273, Tc = 20)

points(tpcdata$r ~ tpcdata$Temp, col="gray",cex=1)
TdegK <- seq(5,50,0.1)+273
TdegC <- seq(5,50,0.1)
lines(schoolfield_high(coef(fits)[1], coef(fits)[2], coef(fits)[3], coef(fits)[4], temp = TdegK, Tc = 20) ~ TdegC, lwd=2)
box(lwd=2,bty="l")
axis(1,at=c(19,22,25,28),tck=0.01,padj=-1.25)
#mtext("Genetic Diversity (number of clones)",1,line=1.8,cex=1.1)
axis(2,at=c(-1,0,1,2,3),tck=0.01,las=TRUE,hadj=-0.5)
abline(a=0,b=0,lty=2)

dev.off()

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
			mutate(Combination=as.factor(Combination), Temperature=as.factor(Temperature))%>%
			#mutate(Star=as.factor(as.numeric(Combination=="B" | Combination=="D")))
			mutate(Background=ifelse(Combination=="B" | Combination=="D","A",ifelse(Combination=="A","B","C")))

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

## Post-hoc differences among clones across temperatures.
#################################################################################################
## Author's note: in some computers we have noticed an eammeans_test() error "Nonconforming number of contrast coefficients", which results in this analysis aborting before completing. We have not been able to understand why this sometimes works and why it doesn't, whether it's a problem with the package version or some other type of version conflict. However, we have been able to reproduce these results in two different computers, and report the R output of this analysis in SI.
#################################################################################################  
pwc <- tdata_single_cl %>% 
  group_by(Temperature) %>%
  emmeans_test(
    growth ~ Combination,
    p.adjust.method = "bonferroni"
    )
pwc %>% print(n=120)


## Effect of star genomes.
## Author's note: this analysis works well.
tdata_single_cl %>%
  anova_test(
    growth ~ Background*Temperature
  )
  summary(lm(growth ~ Background*Temperature,data=tdata_single_cl))

pwc <- tdata_single_cl %>% 
  group_by(Temperature) %>%
  emmeans_test(
    growth ~ Background,
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


##------------------------------------------------------------------------------------------
## 2) FIGURE 2c

## Prep data
tdata_new_star <- tdata_new %>%
	group_by(Clone_number) %>%
	filter(Clone_number==1) %>%
	mutate(Background=ifelse(Combination=="B" | Combination=="D","A",ifelse(Combination=="A","B","C"))) %>%
	group_by(Temperature,Background) %>%
	summarize(mean_growth=mean(growth, na.rm=TRUE),se_A=sd(growth, na.rm=TRUE)/sqrt(8),se_B=sd(growth, na.rm=TRUE)/sqrt(4),se_C=sd(growth, na.rm=TRUE)/sqrt(8)) 

## Offset the means in the horizontal axis for better visualization
A_plot <- data.frame(mean_growth=tdata_new_star[which(tdata_new_star$Background=="A"),]$mean_growth,Temperature=tdata_new_star[which(tdata_new_star$Background=="A"),]$Temperature-0.25,Background=rep("A",4),se=tdata_new_star[which(tdata_new_star$Background=="A"),]$se_A)

B_plot <- data.frame(mean_growth=tdata_new_star[which(tdata_new_star$Background=="B"),]$mean_growth,Temperature=tdata_new_star[which(tdata_new_star$Background=="B"),]$Temperature,Background=rep("B",4),se=tdata_new_star[which(tdata_new_star$Background=="B"),]$se_B)

C_plot <- data.frame(mean_growth=tdata_new_star[which(tdata_new_star$Background=="C"),]$mean_growth,Temperature=tdata_new_star[which(tdata_new_star$Background=="C"),]$Temperature+0.25,Background=rep("C",4),se=tdata_new_star[which(tdata_new_star$Background=="C"),]$se_C)

finalplot <- rbind(A_plot,B_plot,C_plot)

pdf("~/Desktop/JP/Papers_in_review_submitted/Alex_tetra/Manuscript/Figures/Fig_2/fig_2c.pdf",height=5.5,width=12)
## Figure
par(mfrow=c(1,2))
plot(mean_growth~Temperature,data=finalplot[which(finalplot$Background=="A"),],pch=16,type="b",xlim=c(18,29),ylim=c(-1,3),col=viridis(72)[71],axes=FALSE,xlab="",ylab="", lwd=2)
box(lwd=1.5,bty="l")
points(mean_growth~Temperature,data=finalplot[which(finalplot$Background=="B"),],pch=16,type="b",col=viridis(72)[52], lwd=2)
points(mean_growth~Temperature,data=finalplot[which(finalplot$Background=="C"),],pch=16,type="b",col=viridis(72)[15], lwd=2)


axis(1,(c(19,22,25,28)),labels=c("19","22","25","28"),tck=0.01,padj=-1)
mtext("Temperature",1,line=1.7,cex=1.1)
mtext("Intrinsic Growth Rate (r)",2,line=2,cex=1.1)
axis(2,c(-2,-1,0,1,2,3),tck=0.01,hadj=0.5,las=TRUE)

arrows(x0=finalplot$Temperature, 
y0=finalplot$mean_growth-finalplot$se, 
x1=finalplot$Temperature, 
y1=finalplot$mean_growth+finalplot$se, 
code=3, angle=90, length=0.025, col=c(viridis(72)[72],viridis(72)[52],viridis(72)[15])[as.factor(finalplot$Background)], lwd=1)
dev.off()


##THE END