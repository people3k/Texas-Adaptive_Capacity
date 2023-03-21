library(splines)
library(reshape2)
library(ggplot2)
library(ggfortify)
library(nlme)
library(rcarbon)
library(cowplot)

require(remotes)
remotes::install_version("sf", version = "1.0-9")
###SPDs of population components

#Load the data and subset for Cental Texas
box<- read.csv("raw/FinalRCDTexas2.csv")
box2<- subset(box, Region=="CTx")

###Calibrate the radiocarbon ages
cptcal <- calibrate(x = box2$DatesUse,  errors = box2$Error, calCurves = "intcal20",  normalised = FALSE)
boxbins <- binPrep(sites = box2$Trinomial, ages = box2$DatesUse, h = 100)

###Bin the ages by site trinomial. Correcting for site ``over sampling." Note that h = 100,
#See rcarbon description for more details. 
boxbins <- binPrep(sites = box2$Trinomial, ages = box2$DatesUse, h = 100)

##Run SPD for central Texas to fit to each repective component logistic from 12500 to 8400
spd.CTx <- spd(cptcal, bins=boxbins, runm=100, timeRange=c(12500,8400))
plot(spd.CTx, runm=100, xlim=c(12500,8400), type="simple")

Den<-spd.CTx$grid$PrDens
Time<-spd.CTx$grid$calBP

#Run and Plot null model logistic
logFit <- nls(PrDens~SSlogis(calBP, Asym, xmid, scale), data=spd.CTx$grid,control=nls.control(maxiter=300),start=list(Asym=0.1, xmid=10700, scale=-100))
# Generate a data frame containing the fitted values
logFitDens=data.frame(calBP=spd.CTx$grid$calBP,PrDens=SSlogis(input=spd.CTx$grid$calBP,Asym=coefficients(logFit)[1],xmid=coefficients(logFit)[2],scal=coefficients(logFit)[3]))
# Use the modelTest function (returning the raw simulation output - see below)
LogNull <- modelTest(cptcal, errors=box$Error,nsim=1000,
                     timeRange=c(12500, 8400), model="custom",predgrid=logFitDens, bins=boxbins, runm=100, raw=TRUE)

###Return fit of logistic model and estimated parameters
summary(logFit)
cor(spd.CTx$grid$PrDens,predict(logFit))


##Global significance test of the fit of the three parameter logistic model
round(LogNull$pval,4) #p-value

###Data.Frame of simulation
Con<-cbind(LogNull$result,predict(logFit))
write.table(Con, file = "CTxPhase1Confit.csv", sep = ",")

###Plot SPD for Component #1
Con1<- read.csv("CTxPhase1Confit.csv")

p <- ggplot(Con1,aes(x=(calBP), y=(PrDens))) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_line(aes(color=PhaseID)) + 
  geom_line(aes(y=LogisticFit), color="blue", size=1.1) +
  theme_bw() +
  scale_x_reverse(breaks=c(12000, 11000, 10000, 9000, 8000), limits=c(12500,8400))+
  #scale_y_continuous(limits=c(0,0.42))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text( 
          size=28), plot.title = element_text(size=20, face = "bold"))+
  labs(x = "Years Cal BP", y="SPD of radiocarbon ages", title = "A. Central Tx. Component I")+
  geom_vline(xintercept = 10023)+
  #geom_vline(xintercept = 1300)+
  geom_vline(xintercept = 8650)
p


####Run analysis for component 2 logistic 8400 to 3350=======================================================
#============================================================================================================
spd.CTx <- spd(cptcal, bins=boxbins, runm=100, timeRange=c(8400,3350))
plot(spd.CTx, runm=100, xlim=c(8400,3350), type="simple")

Den<-spd.CTx$grid$PrDens
Time<-spd.CTx$grid$calBP

######Four parameter logistic model
#Plot null model logistic
logFit2 <- nls(PrDens~SSfpl(calBP, A, B, xmid, scale), data=spd.CTx$grid,control=nls.control(maxiter=300),start=list(A=0.03, B=0.15, xmid=5900, scale=-100))
# Generate a data frame containing the fitted values
logFitDens2=data.frame(calBP=spd.CTx$grid$calBP,PrDens=SSfpl(input=spd.CTx$grid$calBP,A=coefficients(logFit2)[1],B=coefficients(logFit2)[2],xmid=coefficients(logFit2)[3],scal=coefficients(logFit2)[4]))
# Use the modelTest function (returning the raw simulation output - see below)
LogNull2 <- modelTest(cptcal, errors=box2$Error,nsim=1000,
                      timeRange=c(8400, 3350), model="custom",predgrid=logFitDens2, bins=boxbins, runm=100, raw=TRUE)

###Plot model over data and return logfit2 model summary for four parameterlogistic 
plot(LogNull2, xlim = c(8400,3350))
lines(Time,predict(logFit2),col="red",lty=2,lwd=3)
summary(logFit2)

plot(Den~Time)


##Global significance test of the fit of the 4 parameter logistic
round(LogNull2$pval,4) #p-value

###Data.Frame of simulation
Con<-cbind(LogNull2$result,predict(logFit2))
write.table(Con, file = "Collapse/CTxPhase2ConfitLogfit.csv", sep = ",")

###GGplot of component 2
Con2<- read.csv("CTxPhase2ConfitLogfit.csv")

p2 <- ggplot(Con2,aes(x=(calBP), y=(PrDens))) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_line(aes(color=PhaseID)) + 
  scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  geom_line(aes(y=logFit2), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse(breaks=c(8500, 7500, 6500, 5500, 4500, 3500), limits=c(8400,3350))+
  # scale_y_continuous(limits=c(0,0.42))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text( 
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="SPD of radiocarbon ages", title = "D. Population Growth Component II")+
  geom_vline(xintercept = 5789)
#geom_vline(xintercept = 1300)+
#geom_vline(xintercept = 600)
p2


####Run analysis for component 3 logistic 3400 to 150
spd.CTx <- spd(cptcal, bins=boxbins, runm=100, timeRange=c(3400,200))
plot(spd.CTx, runm=100, xlim=c(3400,200), type="simple")

Den<-spd.CTx$grid$PrDens
Time<-spd.CTx$grid$calBP

######Four parameter logistic model
#Plot null model logistic
logFit3 <- nls(PrDens~SSfpl(calBP, A, B, xmid, scale), data=spd.CTx$grid,control=nls.control(maxiter=700),start=list(A=0.12, B=0.35, xmid=1711, scale=-100))
# Generate a data frame containing the fitted values
logFitDens3=data.frame(calBP=spd.CTx$grid$calBP,PrDens=SSfpl(input=spd.CTx$grid$calBP,A=coefficients(logFit3)[1],B=coefficients(logFit3)[2],xmid=coefficients(logFit3)[3],scal=coefficients(logFit3)[4]))
# Use the modelTest function (returning the raw simulation output - see below)
LogNull3 <- modelTest(cptcal, errors=box2$Error,nsim=1000,
                      timeRange=c(3400, 200), model="custom",predgrid=logFitDens3, bins=boxbins, runm=100, raw=TRUE)

##Plot 4 parameter logistic over data and return logfit3 model summary
plot(LogNull3, xlim = c(3400,200))
lines(Time,predict(logFit3),col="red",lty=2,lwd=3)
summary(logFit3)

##Global significance test of the fit of the 4 parameter logistic model
round(LogNull3$pval,4) #p-value

Pre1<-predict(logFit3)

##Data.Frame of simulation
Con<-cbind(LogNull3$result,predict(logFit3))
write.table(Con, file = "Collapse/CTxPhase3ConfitLogfit.csv", sep = ",")

###GGplot of component #3

Con3<- read.csv("CTxPhase3ConfitLogfit.csv")

p3 <- ggplot(Con3,aes(x=(calBP), y=(PrDens))) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_line(aes(color=PhaseID)) + 
  scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3500,200))+
  # scale_y_continuous(limits=c(0,0.42))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text( 
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="SPD of radiocarbon ages", title = "H. Population Growth Component III")+
  geom_vline(xintercept = 1800)+
  #geom_vline(xintercept = 1300)+
  geom_vline(xintercept = 650)
p3

pct <- ggplot(Con3,aes(x=(calBP), y=(PrDens))) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(color=DeltaN), size=2.5) +
  scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3500,250))+
  # scale_y_continuous(limits=c(0,0.42))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text( 
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="SPD of radiocarbon ages", title = "B. Central Texas Population and Delta 15N Collagen")+
  geom_vline(xintercept = 1800)
  #geom_vline(xintercept = 1300)+
  #geom_vline(xintercept = 650)
pct

#####All three phases combined

Con3<- read.csv("CTxThreePhaseCompined.csv")

p3 <- ggplot(Con3,aes(x=(calBP), y=(PrDens))) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(color=PhaseID)) + 
  scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  geom_line(aes(y=logFit2), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse(limits=c(12500,200))+
  # scale_y_continuous(limits=c(0,0.42))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text( 
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="SPD of radiocarbon ages", title = "A. Central Texas Three Logistic Components")
  #geom_vline(xintercept = 1800)+
  #geom_vline(xintercept = 1300)+
 # geom_vline(xintercept = 650)
p3

p4 <- ggplot(Con3,aes(x=(calBP), y=(PrDens))) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(color=Violence)) + 
  scale_color_gradient(low ="#619CFF", high = "#F8766D") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  geom_line(aes(y=logFit2), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse(limits=c(12500,200))+
  # scale_y_continuous(limits=c(0,0.42))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text( 
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="SPD of radiocarbon ages", title = "Central Texas Population Growth Components")
#geom_vline(xintercept = 1800)+
#geom_vline(xintercept = 1300)+
# geom_vline(xintercept = 650)
p4

p5 <- ggplot(Con3,aes(x=(calBP), y=(PrDens))) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(color=MeanFCRArea)) + 
  scale_color_gradient(low ="#619CFF", high = "#F8766D") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  geom_line(aes(y=logFit2), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse(limits=c(12500,200))+
  # scale_y_continuous(limits=c(0,0.42))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text( 
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="SPD of radiocarbon ages", title = "Central Texas Population Growth Components")
#geom_vline(xintercept = 1800)+
#geom_vline(xintercept = 1300)+
# geom_vline(xintercept = 650)
p5

p5 <- ggplot(Con3,aes(x=(calBP), y=(PrDens))) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(color=BRMspd)) + 
  scale_color_gradient(low ="#619CFF", high = "#F8766D") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  geom_line(aes(y=logFit2), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse(limits=c(12500,200))+
  # scale_y_continuous(limits=c(0,0.42))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text( 
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="SPD of radiocarbon ages", title = "Central Texas Population Growth Components")
#geom_vline(xintercept = 1800)+
#geom_vline(xintercept = 1300)+
# geom_vline(xintercept = 650)
p5

library(cowplot)
Figp<-plot_grid(p3, p, p2, ncol=3, align="h", axis = "rl")
Figp

pdf("Collapse/LogisticComparison.pdf", width=26, height=13.55)
Figp
dev.off()

#####Next, Need to Sum-up the SPDs and Check Spearmans correlations and mutual information with Paleo-climate NPP and just climate paramters
###Assess mutual information of climate-PrDens series by growth components and by phases
#Sum-up the SPDs fit separately to each component of growth.


##Basic Graphs of NPP vs. Time and Density
keep<-read.csv(file="SPD50.csv", header=T)
keep2<-subset(keep, CompID1=="Component I")
keep3<-subset(keep, CompID1=="Component II")
keep4<-subset(keep, CompID1=="Component III")

p5 <- ggplot(keep,aes(x=(calBP), y=(PrDens))) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(color=CRR), size=5) + 
  scale_color_gradient(low ="#619CFF", high = "#F8766D") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  geom_line(aes(y=logFit), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse(limits=c(12500,200))+
  # scale_y_continuous(limits=c(0,0.42))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text( 
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="SPD of radiocarbon ages", title = "Central Texas Population Growth Components")
#geom_vline(xintercept = 1800)+
#geom_vline(xintercept = 1300)+
# geom_vline(xintercept = 650)
p5


###Spearman correlation by component and phase. First subset each component dataset and then run correlation
Phase<-subset(keep4, PhaseID=="Recession")
cor.test(Phase$NPP, Phase$PrDens, method= "spearman")

##Scatter plots of NPP and SPD values by component of growth

nppAs2 <- ggplot(keep3,aes(x=(NPP), y=(PrDens))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(color=PhaseID)) + 
  #geom_line(aes(y=logFit2, color=factor(CompID)), size=.5) +
  theme_bw() +
  #scale_y_continuous(breaks=c(), limits=c(12500,8400))+
  scale_x_continuous(breaks=c(650, 700, 750, 800, 850), limits=c(625,875))+
  facet_wrap(factor(CompID1)~.)+
  scale_y_continuous(breaks=c(0,5,10,15,20), limits=c(0,20))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text( 
          size=28), plot.title = element_text(size=20, face = "bold"))+
  labs(x = "Modeled NPP", y="SPD of radiocarbon ages", title = "B. Component II NPP and SPD")+
  #geom_vline(xintercept = 10023)+
  geom_smooth(se=FALSE)
nppAs2

nppAs3 <- ggplot(keep4,aes(x=(NPP), y=(PrDens))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(color=PhaseID)) + 
  #geom_line(aes(y=logFit2, color=factor(CompID)), size=.5) +
  theme_bw() +
  scale_y_continuous(breaks=c(0,5,10,15,20), limits=c(0,20))+
  scale_x_continuous(breaks=c(650, 700, 750, 800, 850), limits=c(625,875))+
  facet_wrap(factor(CompID1)~.)+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text( 
          size=28), plot.title = element_text(size=20, face = "bold"))+
  labs(x = "Modeled NPP", y="SPD of radiocarbon ages", title = "C. Component III NPP and SPD")+
  #geom_vline(xintercept = 10023)+
  geom_smooth(se=FALSE)
nppAs3

nppAs1 <- ggplot(keep2,aes(x=(NPP), y=(PrDens))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(color=PhaseID)) + 
  #geom_line(aes(y=logFit2, color=factor(CompID)), size=.5) +
  theme_bw() +
  scale_y_continuous(breaks=c(0,5,10,15,20), limits=c(0,20))+
  scale_x_continuous(breaks=c(650, 700, 750, 800, 850), limits=c(625,875))+
  facet_wrap(factor(CompID1)~.)+
  # scale_y_continuous(limits=c(0,0.42))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text( 
          size=28), plot.title = element_text(size=20, face = "bold"))+
  labs(x = "Modeled NPP", y="SPD of radiocarbon ages", title = "A. Component I NPP and SPD")+
  #geom_vline(xintercept = 10023)+
  geom_smooth(se=FALSE)
nppAs1

library(cowplot)
Figp<-plot_grid(nppAs1, nppAs2, nppAs3, ncol=3, align="vh", axis = "rl")
Figp

####Basic plots of Estimated NPP and SPDs over time


p1npp <- ggplot(keep2,aes(x=(calBP), y=(NPP))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_line(aes(color=PhaseID)) + 
  #geom_line(aes(y=logFit2, color=factor(CompID)), size=.5) +
  theme_bw() +
  scale_x_reverse(breaks=c(12000, 11000, 10000, 9000, 8000), limits=c(12500,8400))+
  scale_y_continuous(breaks=c(650, 700, 750, 800, 850), limits=c(625,875))+
  facet_wrap(factor(CompID1)~.)+
  # scale_y_continuous(limits=c(0,0.42))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text( 
          size=28), plot.title = element_text(size=20, face = "bold"))+
  labs(x = "Years Cal BP", y="Modeled NPP", title = "A. Component I Modeled NPP")+
  geom_vline(xintercept = 10023)+
  geom_smooth(se=FALSE)
p1npp

p1 <- ggplot(keep2,aes(x=(calBP), y=(PrDens))) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_line(aes(color=PhaseID)) + 
  geom_line(aes(y=logFit2), color="blue", size=.5) +
  theme_bw() +
  scale_x_reverse(breaks=c(12000, 11000, 10000, 9000, 8000), limits=c(12500,8400))+
  facet_wrap(factor(CompID1)~.)+
  # scale_y_continuous(limits=c(0,0.42))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text( 
          size=28), plot.title = element_text(size=20, face = "bold"))+
  labs(x = "Years Cal BP", y="SPD of radiocarbon ages", title = "B. Component I Population Growth")+
  geom_vline(xintercept = 10023)
#geom_smooth(se=FALSE)
p1


p2npp <- ggplot(keep3,aes(x=(calBP), y=(NPP))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_line(aes(color=PhaseID)) + 
  #geom_line(aes(y=logFit2, color=factor(CompID)), size=.5) +
  theme_bw() +
  scale_x_reverse(breaks=c(8500, 7500, 6500, 5500, 4500, 3500), limits=c(8400,3350))+
  scale_y_continuous(breaks=c(650, 700, 750, 800, 850), limits=c(625,875))+
  facet_wrap(factor(CompID1)~.)+
  # scale_y_continuous(limits=c(0,0.42))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text( 
          size=28), plot.title = element_text(size=20, face = "bold"))+
  
  labs(x = "Years Cal BP", y="Modeled NPP", title = "C. Component II Modeled NPP")+
  geom_vline(xintercept = 5789)+
  geom_smooth(se=FALSE)
p2npp

p2 <- ggplot(keep3,aes(x=(calBP), y=(PrDens))) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_line(aes(color=PhaseID)) + 
  geom_line(aes(y=logFit2), color="blue", size=.5) +
  theme_bw() +
  scale_x_reverse(breaks=c(8500, 7500, 6500, 5500, 4500, 3500), limits=c(8400,3350))+
  facet_wrap(factor(CompID1)~.)+
  # scale_y_continuous(limits=c(0,0.42))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text( 
          size=28), plot.title = element_text(size=20, face = "bold"))+
  labs(x = "Years Cal BP", y="SPD of radiocarbon ages", title = "D. Component II Population Growth")+
  geom_vline(xintercept = 5789)
#geom_smooth(se=FALSE)
p2


p3npp <- ggplot(keep4,aes(x=(calBP), y=(NPP))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_line(aes(color=PhaseID)) + 
  #geom_line(aes(y=logFit2, color=factor(CompID)), size=.5) +
  theme_bw() +
  scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3500,200))+
  scale_y_continuous(breaks=c(650, 700, 750, 800, 850), limits=c(625,875))+
  facet_wrap(factor(CompID1)~.)+
  # scale_y_continuous(limits=c(0,0.42))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text( 
          size=28), plot.title = element_text(size=20, face = "bold"))+
  
  labs(x = "Years Cal BP", y="Modeled NPP", title = "E. Component III Modeled NPP")+
  geom_vline(xintercept = 1800)+
  geom_smooth(se=FALSE)
p3npp

p3 <- ggplot(keep4,aes(x=(calBP), y=(PrDens))) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_line(aes(color=PhaseID)) + 
  geom_line(aes(y=logFit2), color="blue", size=.5) +
  theme_bw() +
  scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3500,200))+
  facet_wrap(factor(CompID1)~.)+
  # scale_y_continuous(limits=c(0,0.42))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text( 
          size=28), plot.title = element_text(size=20, face = "bold"))+
  labs(x = "Years Cal BP", y="SPD of radiocarbon ages", title = "F. Component III Population Growth")+
  geom_vline(xintercept = 1800)
#geom_smooth(se=FALSE)
p3

library(cowplot)
Figp<-plot_grid(p1npp, p2npp, p3npp, p1, p2,p3, ncol=3, align="vh", axis = "rl")
Figp

###Calculate mutual information for 50 year sums

####Mutual Information
install.packages("remotes")
remotes::install_github("mdscheuerell/muti")

library(muti)
###Run MI for each data frame--keep2, keep3, and keep4, Components 1,2,and3 repectively
muti(keep2$NPP, keep2$PrDens, sym = TRUE, n_bins = NULL, normal = TRUE, lags = seq(-6, 6),
     mc = 100, alpha = 0.1)

####Isotope data analysis

######Box Plots and Violin Plots for Isotopes by Copial I, II and Malthusian Collapse Phases.
d2 <- read.csv("CTexIsotopeSubmit.csv")
dCII<- subset(d2,CompID1=="Component II")
dCIII<-subset(d2,CompID1=="Component III")



pcar1 <- ggplot(dCII, aes(factor(CompID2), (DeltaC13Car), fill=CompID2))+
  geom_violin()+
  geom_boxplot(notch = FALSE, width=0.2)+
  stat_boxplot(geom ='errorbar')+
  scale_fill_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  facet_wrap( ~ factor(CompID1))+
  scale_y_continuous(limits=c(-15.5, -4.5))+
  #scale_fill_manual(values=c("#FC4E07", "#00A4CCFF"))+
  theme_bw() +
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text( 
          size=28),  plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Predicted phase of growth", y="Delta 13C apatite", title = "A. Component II Delta 13C Apatite by Predicted Growth Phase ")
pcar1

pcol1 <- ggplot(dCII, aes(factor(CompID2), (DeltaC13Col), fill=CompID2))+
  geom_violin()+
  geom_boxplot(notch = FALSE, width=0.2)+
  stat_boxplot(geom ='errorbar')+
  scale_fill_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  scale_y_continuous(limits=c(-21, -9))+
  #facet_wrap( ~ factor(CompID1))+
  #scale_fill_manual(values=c("#FC4E07", "#00A4CCFF"))+
  theme_bw() +
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text( 
          size=28),  plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Predicted phase of growth", y="Delta 13C collagen", title = "B. Component II Delta 13C Collagen by Predicted Phase of Growth")
pcol1

pn1 <- ggplot(dCII, aes(factor(CompID2), (DeltaN), fill=CompID2))+
  geom_violin()+
  geom_boxplot(notch = FALSE, width=0.2)+
  scale_fill_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  stat_boxplot(geom ='errorbar')+
  #facet_wrap( ~ factor(CompID1))+
  scale_y_continuous(limits=c(6, 13.5))+
  #scale_fill_manual(values=c("#FC4E07", "#00A4CCFF"))+
  theme_bw() +
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text( 
          size=28),  plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Predicted phase of growth", y="Delta 15N collagen", title = "C. Component II Delta 15N by Predicted Phase of Growth")
pn1


pcar2 <- ggplot(dCIII, aes(factor(CompID2), (DeltaC13Car), fill=CompID2))+
  geom_violin()+
  geom_boxplot(notch = FALSE, width=0.2)+
  stat_boxplot(geom ='errorbar')+
  scale_fill_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  facet_wrap( ~ factor(CompID1))+
  scale_y_continuous(limits=c(-15.5, -4.5))+
  #scale_fill_manual(values=c("#FC4E07", "#00A4CCFF"))+
  theme_bw() +
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text( 
          size=28),  plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Predicted phase of growth", y="Delta 13C apatite", title = "E. Component III Delta 13C Apatite by Predicted Growth Phase ")
pcar2

pcol2 <- ggplot(dCIII, aes(factor(CompID2), (DeltaC13Col), fill=CompID2))+
  geom_violin()+
  geom_boxplot(notch = FALSE, width=0.2)+
  stat_boxplot(geom ='errorbar')+
  scale_fill_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  scale_y_continuous(limits=c(-21, -9))+
  #facet_wrap( ~ factor(CompID1))+
  #scale_fill_manual(values=c("#FC4E07", "#00A4CCFF"))+
  theme_bw() +
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text( 
          size=28),  plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Predicted phase of growth", y="Delta 13C collagen", title = "F. Component III Delta 13C Collagen by Predicted Phase of Growth")
pcol2

pn2 <- ggplot(dCIII, aes(factor(CompID2), (DeltaN), fill=CompID2))+
  geom_violin()+
  geom_boxplot(notch = FALSE, width=0.2)+
  scale_fill_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  stat_boxplot(geom ='errorbar')+
  #facet_wrap( ~ factor(CompID1))+
  scale_y_continuous(limits=c(6, 13.5))+
  #scale_fill_manual(values=c("#FC4E07", "#00A4CCFF"))+
  theme_bw() +
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text( 
          size=28),  plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Predicted phase of growth", y="Delta 15N collagen", title = "G. Component III Delta 15N by Predicted Phase of Growth")
pn2

library(cowplot)

Fig1Rev<-plot_grid(pcar1,pcar2, pcol1, pcol2,pn1, pn2, p2, p3, ncol=2, align="hv", axis = "rl")
Fig1Rev

pdf("Collapse/BoxPlots3.pdf", width=21, height=17.55)
Fig1Rev
dev.off()


####Pairwise Mann-Whittney U tests
##Component II
dCIICopial<- subset(d2, CompID1=="Component II" & WilcoxID<3)
dCIIMalthus<-subset(d2,CompID1=="Component II" & WilcoxID>1)

###13C Apatite
wilcox.test(DeltaC13Car ~ PhaseComp, data=dCIICopial) 
wilcox.test(DeltaC13Car ~ PhaseComp, data=dCIIMalthus) 

##13C Collagen
wilcox.test(DeltaC13Col ~ PhaseComp, data=dCIICopial) 
wilcox.test(DeltaC13Col ~ PhaseComp, data=dCIIMalthus) 

###15N Collage
wilcox.test(DeltaN ~ PhaseComp, data=dCIICopial) 
wilcox.test(DeltaN ~ PhaseComp, data=dCIIMalthus) 

###Component III===============================================================
dCIIICopial2<- subset(d2, CompID1=="Component III" & WilcoxID<3)
dCIIIMalthus2<-subset(d2,CompID1=="Component III" & WilcoxID>1)

###13C Apatite
wilcox.test(DeltaC13Car ~ PhaseComp, data=dCIIICopial2) 
wilcox.test(DeltaC13Car ~ PhaseComp2, data=dCIIIMalthus2) 

##13C Collagen
wilcox.test(DeltaC13Col ~ PhaseComp, data=dCIIICopial2) 
wilcox.test(DeltaC13Col ~ PhaseComp2, data=dCIIIMalthus) 

###15N Collage
wilcox.test(DeltaN ~ PhaseComp, data=dCIIICopial2) 
wilcox.test(DeltaN ~ PhaseComp2, data=dCIIIMalthus2) 

###Cooking Feature SPD and Median Age Calibration
SurArea<- read.csv("raw/SurfaceAreaRaw.csv")

library(dplyr)
#####Average raw radiocarbon ages by site and feature=======================================
MeanSD<-SurArea %>%
  group_by(Site) %>%
  mutate(MTime = mean(RCYBP),
         SDTime = mean(Error))

by_vs_am <- SurArea %>% group_by(Site, Feature)
by_vs <- by_vs_am %>% summarise(n = n(), RCYBP = mean(RCYBP),Error=mean(Error), SurfaceArea=mean(SurfaceArea), RockNumber=mean(RockNumber), RockWeight=mean(RockWeight))
write.table(by_vs, file = "raw/SurfaceMeans2.csv", sep = ",")

#1 Calculate SPD for Dated Earth Oven and Burned Rock Midden Features in Central Texas

##Load rcarbon
library(rcarbon)

#Load the data and subset for Cental Texas
box<- read.csv("raw/SurfaceAreaRaw.csv")
box1 <- subset(box, EarthOven=="Yes")

###Calibrate the radiocarbon ages
cptcal <- calibrate(x = box$RCYB,  errors = box$Error, calCurves = "intcal20",  normalised = FALSE)

####Calculate median ages of all dated features
MedCal<-medCal(cptcal)
write.table(MedCal, file = "raw/MedCalSurface.csv", sep = ",")

#boxbins <- binPrep(sites = box2$Trinomial, ages = box2$DatesUse, h = 100)

##Run SPD for central Texas 
spd.CTx <- spd(cptcal, bins=NA, runm=300, timeRange=c(12500,200))
plot(spd.CTx, runm=300, xlim=c(12500,200), type="simple")
write.table(spd.CTx, file = "NewEarthOvenSPD.csv", sep = ",")
 
####Correlation between surface area and FCR
SurArea2<- read.csv("raw/SurfaceMeans2.csv")

cor.test(SurArea2$RockNumber, SurArea2$SurfaceArea, method= "spearman")

####Graphs of Population dynamics and Earth Oven Dynamics
library(RColorBrewer)   

Con1<- read.csv("CTxPhase1Confit.csv")

p <- ggplot(Con1,aes(x=(calBP), y=(PrDens))) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_line(aes(color=PhaseID),size=2.5) + 
  geom_line(aes(y=LogisticFit), color="blue", size=1.1) +
  theme_bw() +
  scale_color_manual(values=c("#619CFF","#00BA38","#F8766D"))+
  scale_x_reverse(breaks=c(12000, 11000, 10000, 9000, 8000), limits=c(12500,8400))+
  #scale_y_continuous(limits=c(0,0.42))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text( 
          size=28), plot.title = element_text(size=20, face = "bold"))+
  labs(x = "Years Cal BP", y="SPD of radiocarbon ages", title = "A. Component I Population Growth")+
  geom_vline(xintercept = 10023)
#geom_vline(xintercept = 1300)+
#geom_vline(xintercept = 8650)
p


br <- ggplot(Con1,aes(x=(calBP), y=(BRMspd))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_line(aes(color=(MeanFCRArea)), size=2.5) + 
  # scale_color_viridis(direction=-1)+
  scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  theme_bw() +
  scale_x_reverse(breaks=c(12000, 11000, 10000, 9000, 8000), limits=c(12500,8400))+
  #scale_y_continuous(limits=c(0,0.42))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text( 
          size=28), plot.title = element_text(size=20, face = "bold"))+
  labs(x = "Years Cal BP", y="SPD of earth oven radiocarbon ages", title = "B. Component I Earth Oven SPD and Surface Area")+
  geom_vline(xintercept = 10023)
#geom_vline(xintercept = 1300)+
# geom_vline(xintercept = 8650)
br



##Component II
Con2<- read.csv("CTxPhase2ConfitLogfit.csv")

p2 <- ggplot(Con2,aes(x=(calBP), y=(PrDens))) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_line(aes(color=PhaseID), size=2.5) + 
  geom_line(aes(y=logFit2), color="blue", size=1.1) +
  theme_bw() +
  scale_color_manual(values=c("#619CFF","#00BA38","#F8766D"))+
  scale_x_reverse(breaks=c(8500, 7500, 6500, 5500, 4500, 3500), limits=c(8400,3350))+
  # scale_y_continuous(limits=c(0,0.42))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text( 
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="SPD of radiocarbon ages", title = "C. Component II Population Growth")+
  geom_vline(xintercept = 5789)
#geom_vline(xintercept = 4200)
#geom_vline(xintercept = 600)
p2

br2 <- ggplot(Con2,aes(x=(calBP), y=(BRMspd))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_line(aes(color=(MeanFCRArea)), size=2.5) + 
  #scale_color_viridis()+
  scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  theme_bw() +
  #scale_color_manual(values=c("#00BA38","#619CFF","#F8766D"))+
  scale_x_reverse(breaks=c(8500, 7500, 6500, 5500, 4500, 3500), limits=c(8400,3350))+
  #scale_y_continuous(limits=c(0,0.42))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text( 
          size=28), plot.title = element_text(size=20, face = "bold"))+
  labs(x = "Years Cal BP", y="SPD of earth oven radiocarbon ages", title = "D. Component II Earth Oven SPD and Surface Area")+
  geom_vline(xintercept = 5789)
#geom_vline(xintercept = 4200)
br2

###Component III

Con3<- read.csv("CTxPhase3ConfitLogfit.csv")

p3 <- ggplot(Con3,aes(x=(calBP), y=(PrDens))) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_line(aes(color=PhaseID), size=2.5) + 
  geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3500,200))+
  scale_color_manual(values=c("#619CFF","#00BA38","#F8766D"))+
  # scale_y_continuous(limits=c(0,0.42))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text( 
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="SPD of radiocarbon ages", title = "E. Component III Population Growth")+
  geom_vline(xintercept = 1800)
#geom_vline(xintercept = 1300)+
#geom_vline(xintercept = 650)
p3

br3 <- ggplot(Con3,aes(x=(calBP), y=(BRMspd))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_line(aes(color=(MeanFCRArea)), size=2.5) + 
  #  scale_color_viridis()+
  scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  theme_bw() +
  #scale_color_manual(values=c("#00BA38","#619CFF","#F8766D"))+
  scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3500,200))+
  #scale_y_continuous(limits=c(0,0.42))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text( 
          size=28), plot.title = element_text(size=20, face = "bold"))+
  labs(x = "Years Cal BP", y="SPD of earth oven radiocarbon ages", title = "F. Component III Earth Oven SPD and Surface Area")+
  geom_vline(xintercept = 1800)
#geom_vline(xintercept = 650)
br3

library(cowplot)
Figp<-plot_grid(p, p2, p3, br, br2, br3, ncol=3, align="h", axis = "rl")
Figp

pdf("PopBRMComponents.pdf", width=25, height=19.55)
Figp
dev.off()
