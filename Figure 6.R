### This code produces Figure 6 and Table S1 in the manuscript: 
### Burgess SC, Johnston EC, Speare KE, McLachlan RH, 
# Adam TC, Vega Thurber R, Burkepile DE. Differential effects of 
# nutrients and fish consumer pressure on sympatric 
# cryptic coral species (Pocillopora spp.).
# Code finalized July 2024
# Code written by Scott Burgess
# Any comments or error reporting, please contact Scott Burgess: sburgess@bio.fsu.edu

library('dplyr') 
library('tidyr')
library('glmmTMB')
library('emmeans')
library('AICcmodavg')
library('car')
library('DHARMa')


# sessionInfo()
# R version 4.4.1 (2024-04-24)
# Platform: x86_64-apple-darwin20
# Running under: macOS Sonoma 14.5
# other attached packages:
# [1] patchwork_1.2.0  ggplot2_3.5.1   
# [3] emmeans_1.10.2   tidyr_1.3.1     
# [5] car_3.1-2        carData_3.0-5   
# [7] dplyr_1.1.4      AICcmodavg_2.3-3
# [9] glmmTMB_1.1.9    DHARMa_0.4.6  



# set the colors for plotting
cols <- cbind.data.frame(Species = c("P. tuahiniensis",
                                     "P. meandrina",
                                     "P. grandis",
                                     "P. verrucosa",
                                     "P. cf. effusa",
                                     "Unknown"),
                         Color = c("#D55E00",
                                   "#0072B2",
                                   "#56B4E9",
                                   "#E69F00",
                                   "#009E73",
                                   "grey"))




################### Import data ################### 
dat <- read.csv("Predation and Polyp data.csv")
SpeciesID <- read.csv("Species IDs.csv")


################### Prepare data ################### 
# Add species ID's to dat
unique(dat$Un_ID)
unique(SpeciesID$Un_ID)
dat <- cbind.data.frame(dat,
                        SpeciesID[match(dat$Un_ID,SpeciesID$Un_ID),c(2:4,7:12)] # check column indices are correct
)

# these should all be removed (were in spreadsheet but no sample provided to us, or bad sample)
dat %>% filter(flag_missing_sample==1)
dat <- dat %>% filter(flag_missing_sample!=1)

# could not get sequence from these corals, the only sample without a species ID
dat %>% filter(Species=="")
# Call the blank species names "Unknown"
dat$Species <- ifelse(dat$Species=="", "Unknown",dat$Species)
dat %>% filter(Species=="Unknown")


# Rename the herbivory treatments
dat$Herbivory_trt <- ifelse(dat$Herbivory_trt=="1X1","Very low",
                            ifelse(dat$Herbivory_trt=="2X2","Low",
                                   ifelse(dat$Herbivory_trt=="3X3","Medium","High")))
dat$Herbivory_trt <- factor(dat$Herbivory_trt, 
                            levels=c("Very low","Low","Medium","High"))

# How many Pocillopora in total?
nrow(dat) # n=406

# Convert data from wide to long format, to analyze as repeated measure
d <- dat %>% dplyr::select(Un_ID,
                           Plot,
                           Herbivory_trt,
                           Enrichment_trt,
                           Species,
                           percent_polyps_extended_Nov21_T10,
                           Bites_present_Nov21_T10,
                           percent_polyps_extended_Apr22_T11,
                           Bites_present_Apr22_T11,
                           percent_polyps_extended_Jul22_T12,
                           Bites_present_Jul22_T12)
d_long <- tidyr::gather(d, Bite.time.point, Bites_present, Bites_present_Nov21_T10, Bites_present_Apr22_T11, Bites_present_Jul22_T12)
foo <- tidyr::gather(d, polyps.time.point, polyps_extended, percent_polyps_extended_Nov21_T10,percent_polyps_extended_Apr22_T11,percent_polyps_extended_Jul22_T12)
d_long <- cbind.data.frame(d_long[,-8:-8],foo[,9:10])
# names(d_long)

# Change yes/no to 1/0 in Bites_present
d_long$Bites_present <- ifelse(d_long$Bites_present=="yes",1,
                               ifelse(d_long$Bites_present=="no",0,NA))

# Make a block variable for each spatial block
d_long$Block <- substr(d_long$Plot,1,1)


# Function to compare models using AICc 
AnalyzeAICc <- function(dat,yname,family){
  dat$yvar <- dat[,names(dat)==yname,]
  
  m1 <- glmmTMB(yvar ~ Herbivory_trt * Enrichment_trt * Species + (1|Un_ID), family=family, data=dat)
  m2 <- update(m1,.~.-Herbivory_trt:Enrichment_trt:Species)
  m3 <- update(m2,.~.-Herbivory_trt:Enrichment_trt)
  m4 <- update(m2,.~.-Enrichment_trt:Species)
  m5 <- update(m2,.~.-Herbivory_trt:Species)
  m6 <- glmmTMB(yvar ~ Herbivory_trt + Enrichment_trt + Species + Herbivory_trt:Enrichment_trt + (1|Un_ID), family=family, data=dat)
  m7 <- glmmTMB(yvar ~ Herbivory_trt + Enrichment_trt + Species + Enrichment_trt:Species + (1|Un_ID), family=family, data=dat)
  m8 <- glmmTMB(yvar ~ Herbivory_trt + Enrichment_trt + Species + Herbivory_trt:Species + (1|Un_ID), family=family, data=dat)
  m9 <- glmmTMB(yvar ~ Herbivory_trt + Enrichment_trt + Species + (1|Un_ID), family=family, data=dat)
  m10 <- glmmTMB(yvar ~ Herbivory_trt + Enrichment_trt + Herbivory_trt:Enrichment_trt + (1|Un_ID), family=family, data=dat)
  m11 <- glmmTMB(yvar ~ Herbivory_trt + Enrichment_trt + (1|Un_ID), family=family, data=dat)
  aic_table <- aictab(list(m1=m1,m2=m2,m3=m3,m4=m4,m5=m5,m6=m6,m7=m7,m8=m8,m9=m9,m10=m10,m11=m11))
  Anova_table <- Anova(get(aic_table$Modnames[1]),type="III")
  
  pred <- expand.grid(Herbivory_trt=unique(dat$Herbivory_trt),
                      Enrichment_trt=unique(dat$Enrichment_trt),
                      Species=unique(dat$Species),
                      Un_ID=NA)
  p <- predict(get(aic_table$Modnames[1]),newdat=pred,se.fit=T)

if(family=="binomial"){
  pred$fit <- plogis(p$fit)
  pred$lwr95 <- plogis(p$fit - qnorm(1-(0.05/2))*p$se.fit)
  pred$upr95 <- plogis(p$fit + qnorm(1-(0.05/2))*p$se.fit)
  pred$lwr85 <- plogis(p$fit - qnorm(1-(0.15/2))*p$se.fit)
  pred$upr85 <- plogis(p$fit + qnorm(1-(0.15/2))*p$se.fit)}
  
  if(family=="guassian"){
    pred$fit <- (p$fit)
    pred$lwr95 <- (p$fit - qnorm(1-(0.05/2))*p$se.fit)
    pred$upr95 <- (p$fit + qnorm(1-(0.05/2))*p$se.fit)
    pred$lwr85 <- (p$fit - qnorm(1-(0.15/2))*p$se.fit)
    pred$upr85 <- (p$fit + qnorm(1-(0.15/2))*p$se.fit)}
  
  DHARMaOutput <- simulateResiduals(fittedModel = get(aic_table$Modnames[1]), plot = F)
  
  list(pred=pred,
       model=get(aic_table$Modnames[1]),
       DHARMaOutput=DHARMaOutput,
       Anova_table=Anova_table,
       aic_table=aic_table)
}


# Analyze bites data and save results
d <- d_long %>% filter(Species %in% c("P. meandrina", "P. tuahiniensis"))
Results_bites <- AnalyzeAICc(dat=d, yname="Bites_present",family="binomial")
# plot(Results_bites$DHARMaOutput)

# Analyze polyp extension data and save results
d <- d_long %>% filter(Species %in% c("P. meandrina", "P. tuahiniensis"))
d$polyps_extended_01 <- ifelse(d$polyps_extended==0,0,1)
Results_polyps <- AnalyzeAICc(dat=d, yname="polyps_extended_01",family="binomial")
# plot(Results_polyps$DHARMaOutput) # looks ok

############################################################################




################### Results ###################
### Bites
Results_bites$aic_table
Results_bites$Anova_table #: χ^2=5.29, df = 1, p = 0.021

Results_bites$model
EMM <- emmeans(Results_bites$model,~Enrichment_trt:Species)
### Simple pairwise comparisons
# pairs(EMM, simple = "Species")     
predict(Results_bites$model,
        newdat=expand.grid(Herbivory_trt=c("Very low"),
                           Enrichment_trt=c("Ambient","Enriched"),
                           Species=c("P. meandrina","P. tuahiniensis"),
                           Un_ID=NA),type='link')
# Ambient vs Enriched for P. tuahiniensis
# -4.757795 - -3.573985 = -1.18381
# exp(-4.757795); exp(-3.573985)
# [1] 0.008584517 # odds Ambient
# [1] 0.02804388 # odds Enriched
pairs(EMM, simple = "Enrichment_trt")
(exp(-1.184)-1)*100 # odds were 69% (28 - 87, 95% CI) less in ambient compared to nutrients, P. tua
(exp(-1.184 + c(-2,2)*0.425)-1)*100

# Ambient vs Enriched for P. meandrina
# -2.854173 -  -3.096636 = 0.242463
# exp(-2.854173); exp(-3.096636)
# [1] 0.05760344 # odds Ambient
# [1] 0.045201 # odds Enriched
(exp(0.242)-1)*100 
(exp(0.242 + c(-2,2)*0.450)-1)*100

Results_bites$pred %>% filter(Herbivory_trt=="Low", Enrichment_trt=="Ambient")
Results_bites$pred %>% filter(Herbivory_trt=="Low", Enrichment_trt=="Enriched")



### Polyp extension
Results_polyps$aic_table
Results_polyps$Anova_table # : χ^2=9.65, df = 3, p = 0.022


EMM <- emmeans(Results_polyps$model,~Herbivory_trt:Species)
predict(Results_polyps$model,
        newdat=expand.grid(Herbivory_trt=c("Medium"),
                           Enrichment_trt=c("Ambient"),
                           Species=c("P. meandrina","P. tuahiniensis"),
                           Un_ID=NA),type='link')
# -0.1432152 - 2.1937419 = -2.336957
exp(-0.1432152);exp(2.1937419)
# [1] 0.8665676 # Odds P.m
# [1] 8.96871 # Odds Pt
pairs(EMM, simple = "Species")
(exp(-2.337)-1)*100 # 90% (73 - 97, 95% CI)
(exp(-2.337 + c(-2,2)*0.514)-1)*100

pairs(EMM, simple = "Herbivory_trt")

EMM <- emmeans(Results_polyps$model,~Enrichment_trt)
pairs(EMM, simple = "Enrichment_trt")
(exp(-0.537)-1)*100 # odds ? by 42% (-67 - 3, 95% CI) under nutrients
(exp(-0.537 + c(-2,2)*0.281)-1)*100
############################################################################





################### FIGURE 6 ###################
quartz(width=5,height=5)
par(mfrow=c(2,2),mar=c(4,2,1,2),oma=c(1,4,1,1))
ylims <- c(0,1) #30
xlims <- c(0.8,4.2)
xlabs <- unique(Results_bites$pred$Herbivory_trt)
x_Pm <- as.numeric(xlabs)-0.1
x_H10 <- as.numeric(xlabs)+0.1


##### BITES
# a) Ambient
plot(xlims,ylims, type="n", xaxt="n", yaxt="n",bty="l",ylab="",xlab="")

y <- Results_bites$pred %>% filter(Enrichment_trt=="Ambient" & Species=="P. meandrina")
with(y, points(x_Pm, fit,pch=19,cex=1.5,
               col=adjustcolor(cols[2,2],alpha.f=0.8)))	
with(y, segments(x_Pm, lwr95, x_Pm, upr95,
                 col=adjustcolor(cols[2,2],alpha.f=0.8),lwd=2,lend=2))
with(y, segments(x_Pm, lwr85, x_Pm, upr85,
                 col=adjustcolor(cols[2,2],alpha.f=0.8),lwd=4,lend=2))

# P. tuahiniensis
y <- Results_bites$pred %>% filter(Enrichment_trt=="Ambient" & Species=="P. tuahiniensis")
with(y, points(x_H10, fit,pch=19,cex=1.5,
               col=adjustcolor(cols[1,2],alpha.f=0.8)))	
with(y, segments(x_H10, lwr95, x_H10, upr95,
                 col=adjustcolor(cols[1,2],alpha.f=0.8),lwd=2,lend=2))
with(y, segments(x_H10, lwr85, x_H10, upr85,
                 col=adjustcolor(cols[1,2],alpha.f=0.8),lwd=4,lend=2))
axis(side=1,at=1:4,c("","","",""))
# text(1:4,rep(-0.15,4),xlabs,xpd=T,srt=45,adj=1)
axis(side=2,at=seq(0,1,0.2),las=1,cex.axis=1.2)
mtext("a) Ambient",side=3,line=0,adj=0)

mtext("Proportion of colonies\nwith bites", side=2, line=3,at=0.5,adj=0.5,cex=1.2)


# b) Enrichment
plot(xlims,ylims, type="n", xaxt="n", yaxt="n",bty="l",ylab="",xlab="")

# P. meandrina
y <- Results_bites$pred %>% filter(Enrichment_trt=="Enriched" & Species=="P. meandrina")
with(y, points(x_Pm, fit,pch=19,cex=1.5,
               col=adjustcolor(cols[2,2],alpha.f=0.8)))	
with(y, segments(x_Pm, lwr95, x_Pm, upr95,
                 col=adjustcolor(cols[2,2],alpha.f=0.8),lwd=2,lend=2))
with(y, segments(x_Pm, lwr85, x_Pm, upr85,
                 col=adjustcolor(cols[2,2],alpha.f=0.8),lwd=4,lend=2))


# P. tuahiniensis
y <- Results_bites$pred %>% filter(Enrichment_trt=="Enriched" & Species=="P. tuahiniensis")
with(y, points(x_H10, fit,pch=19,cex=1.5,
               col=adjustcolor(cols[1,2],alpha.f=0.8)))	
with(y, segments(x_H10, lwr95, x_H10, upr95,
                 col=adjustcolor(cols[1,2],alpha.f=0.8),lwd=2,lend=2))
with(y, segments(x_H10, lwr85, x_H10, upr85,
                 col=adjustcolor(cols[1,2],alpha.f=0.8),lwd=4,lend=2))
axis(side=1,at=1:4,c("","","",""))
# text(1:4,rep(-0.15,4),xlabs,xpd=T,srt=45,adj=1)
axis(side=2,at=seq(0,1,0.2),las=1,cex.axis=1.2)
mtext("b) Nutrient Enrichment",side=3,line=0,adj=0)

legend(2.5,0.3,legend=c(
  expression(paste(italic("P. meandrina"))),
  expression(paste(italic("P. tuahiniensis")))),
  pch=19,
  pt.cex=1,
  cex=0.9,
  bty="n",
  xpd=T,
  col=c(adjustcolor(cols[2,2],alpha.f=0.8),
        adjustcolor(cols[1,2],alpha.f=0.8)))


##### POLYPS

ylims <- c(0,1) #30
xlims <- c(0.8,4.2)
xlabs <- unique(Results_polyps$pred$Herbivory_trt)
x_Pm <- as.numeric(xlabs)-0.1
x_H10 <- as.numeric(xlabs)+0.1

# a) Ambient
plot(xlims,ylims, type="n", xaxt="n", yaxt="n",bty="l",ylab="",xlab="")

y <- Results_polyps$pred %>% filter(Enrichment_trt=="Ambient" & Species=="P. meandrina")
with(y, points(x_Pm, fit,pch=19,cex=1.5,
               col=adjustcolor(cols[2,2],alpha.f=0.8)))	
with(y, segments(x_Pm, lwr95, x_Pm, upr95,
                 col=adjustcolor(cols[2,2],alpha.f=0.8),lwd=2,lend=2))
with(y, segments(x_Pm, lwr85, x_Pm, upr85,
                 col=adjustcolor(cols[2,2],alpha.f=0.8),lwd=4,lend=2))

# P. tuahiniensis
y <- Results_polyps$pred %>% filter(Enrichment_trt=="Ambient" & Species=="P. tuahiniensis")
with(y, points(x_H10, fit,pch=19,cex=1.5,
               col=adjustcolor(cols[1,2],alpha.f=0.8)))	
with(y, segments(x_H10, lwr95, x_H10, upr95,
                 col=adjustcolor(cols[1,2],alpha.f=0.8),lwd=2,lend=2))
with(y, segments(x_H10, lwr85, x_H10, upr85,
                 col=adjustcolor(cols[1,2],alpha.f=0.8),lwd=4,lend=2))
axis(side=1,at=1:4,c("","","",""))
text(1:4,rep(-0.15,4),xlabs,xpd=T,srt=45,adj=1)
axis(side=2,at=seq(0,1,0.2),las=1,cex.axis=1.2)
mtext("c) Ambient",side=3,line=0,adj=0)

mtext("Proportion of colonies\nwith polyps extended", side=2, line=3,at=0.5,adj=0.5,cex=1.2)


# b) Enrichment
plot(xlims,ylims, type="n", xaxt="n", yaxt="n",bty="l",ylab="",xlab="")

# P. meandrina
y <- Results_polyps$pred %>% filter(Enrichment_trt=="Enriched" & Species=="P. meandrina")
with(y, points(x_Pm, fit,pch=19,cex=1.5,
               col=adjustcolor(cols[2,2],alpha.f=0.8)))	
with(y, segments(x_Pm, lwr95, x_Pm, upr95,
                 col=adjustcolor(cols[2,2],alpha.f=0.8),lwd=2,lend=2))
with(y, segments(x_Pm, lwr85, x_Pm, upr85,
                 col=adjustcolor(cols[2,2],alpha.f=0.8),lwd=4,lend=2))

# P. tuahiniensis
y <- Results_polyps$pred %>% filter(Enrichment_trt=="Enriched" & Species=="P. tuahiniensis")
with(y, points(x_H10, fit,pch=19,cex=1.5,
               col=adjustcolor(cols[1,2],alpha.f=0.8)))	
with(y, segments(x_H10, lwr95, x_H10, upr95,
                 col=adjustcolor(cols[1,2],alpha.f=0.8),lwd=2,lend=2))
with(y, segments(x_H10, lwr85, x_H10, upr85,
                 col=adjustcolor(cols[1,2],alpha.f=0.8),lwd=4,lend=2))
axis(side=1,at=1:4,c("","","",""))
text(1:4,rep(-0.15,4),xlabs,xpd=T,srt=45,adj=1)
axis(side=2,at=seq(0,1,0.2),las=1,cex.axis=1.2)
mtext("d) Nutrient Enrichment",side=3,line=0,adj=0)

mtext("Consumer pressure treatment", side=1, line=0, cex=1.2, outer=T)
########################################################################
