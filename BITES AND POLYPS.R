# Code to analyse results from the 'RECHARGE' plots
# Code written by Scott Burgess July 2023

# R version 4.2.3 (2021-03-31) -- "Shortstop Beagle"
library('dplyr') # v1.2.0
library('tidyr')
library('glmmTMB') # v1.1.3
library('emmeans')
library('AICcmodavg')
library('car')
library('DHARMa')

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




# Import and prepare data
dat <- read.csv("Pred&Polyp_from_Deron 2023_02_24.csv")
SpeciesID <- read.csv("Species IDs.csv")


# Add species ID's to dat
unique(dat$Un_ID)
unique(SpeciesID$Un_ID)
dat <- cbind.data.frame(dat,
                        SpeciesID[match(dat$Un_ID,SpeciesID$Un_ID),c(2:4,11:16)] # check column indices are correct
)

# these should all be removed (were in spreadsheet but no sample provided to us, or bad sample, see notes)
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

### Look at proportion of Pocillopora with bites from each species in each plot
d <- d_long %>% filter(Species %in% c("P. meandrina", "P. tuahiniensis"))

Results <- AnalyzeAICc(dat=d, yname="Bites_present",family="binomial")

# plot(Results$DHARMaOutput) # looks ok

Results$aic_table
Results$Anova_table

Results$model
EMM <- emmeans(Results$model,~Enrichment_trt:Species)
### Simple pairwise comparisons
# pairs(EMM, simple = "Species")     
pairs(EMM, simple = "Enrichment_trt")
(exp(0.242)-1)*100 # odds ? by 27% (-48 - 213, 95% CI), P. mea
(exp(0.242 + c(-2,2)*0.450)-1)*100

(exp(-1.184)-1)*100 # odds increased by 69% (28 - 87, 95% CI) in nutrients, P. tua
(exp(-1.184 + c(-2,2)*0.425)-1)*100
((1/exp(-1.184))-1)*100 # odds increased by 69% (28 - 87, 95% CI), P. tua
(exp(-1.184 + c(-2,2)*0.425)-1)*100

Results$pred %>% filter(Herbivory_trt=="Low", Enrichment_trt=="Ambient")
Results$pred %>% filter(Herbivory_trt=="Low", Enrichment_trt=="Enriched")


##### FIGURE 6 top
quartz(width=5,height=5)
par(mfrow=c(2,2),mar=c(4,2,1,2),oma=c(1,4,1,1))
ylims <- c(0,1) #30
xlims <- c(0.8,4.2)
xlabs <- unique(Results$pred$Herbivory_trt)
x_Pm <- as.numeric(xlabs)-0.1
x_H10 <- as.numeric(xlabs)+0.1

# a) Ambient
plot(xlims,ylims, type="n", xaxt="n", yaxt="n",bty="l",ylab="",xlab="")

y <- Results$pred %>% filter(Enrichment_trt=="Ambient" & Species=="P. meandrina")
with(y, points(x_Pm, fit,pch=19,cex=1.5,
               col=adjustcolor(cols[2,2],alpha.f=0.8)))	
with(y, segments(x_Pm, lwr95, x_Pm, upr95,
                 col=adjustcolor(cols[2,2],alpha.f=0.8),lwd=2,lend=2))
with(y, segments(x_Pm, lwr85, x_Pm, upr85,
                 col=adjustcolor(cols[2,2],alpha.f=0.8),lwd=4,lend=2))

# P. tuahiniensis
y <- Results$pred %>% filter(Enrichment_trt=="Ambient" & Species=="P. tuahiniensis")
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
y <- Results$pred %>% filter(Enrichment_trt=="Enriched" & Species=="P. meandrina")
with(y, points(x_Pm, fit,pch=19,cex=1.5,
               col=adjustcolor(cols[2,2],alpha.f=0.8)))	
with(y, segments(x_Pm, lwr95, x_Pm, upr95,
                 col=adjustcolor(cols[2,2],alpha.f=0.8),lwd=2,lend=2))
with(y, segments(x_Pm, lwr85, x_Pm, upr85,
                 col=adjustcolor(cols[2,2],alpha.f=0.8),lwd=4,lend=2))


# P. tuahiniensis
y <- Results$pred %>% filter(Enrichment_trt=="Enriched" & Species=="P. tuahiniensis")
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

mtext("Fish consumer pressure treatment", side=1, line=0, cex=1.2, outer=T)
########################################################################







### Polyp extension
d <- d_long %>% filter(Species %in% c("P. meandrina", "P. tuahiniensis"))
d$polyps_extended_01 <- ifelse(d$polyps_extended==0,0,1)

Results01 <- AnalyzeAICc(dat=d, yname="polyps_extended_01",family="binomial")
# Results <- AnalyzeAICc(dat=d, yname="polyps_extended",family="gaussian")

# plot(Results01$DHARMaOutput) # looks ok
# plot(Results$DHARMaOutput) # looks bad? 

Results01$aic_table
Results01$Anova_table

Results01$model
EMM <- emmeans(Results01$model,~Herbivory_trt:Species)
### Simple pairwise comparisons
# pairs(EMM, simple = "Species")     
pairs(EMM, simple = "Species")
(exp(-0.948)-1)*100 # odds ? by 27% (-48 - 213, 95% CI), P. mea
(exp(-0.948 + c(-2,2)*0.532)-1)*100

(exp(-2.337))*100 # 10% (3 - 27, 95% CI), P. mea
(exp(-2.337 + c(-2,2)*0.514))*100

pairs(EMM, simple = "Herbivory_trt")

EMM <- emmeans(Results01$model,~Enrichment_trt)
pairs(EMM, simple = "Enrichment_trt")
(exp(-0.537)-1)*100 # odds ? by 42% (-67 - 3, 95% CI) under nutrients
(exp(-0.537 + c(-2,2)*0.281)-1)*100


##### FIGURE bottom

ylims <- c(0,1) #30
xlims <- c(0.8,4.2)
xlabs <- unique(Results01$pred$Herbivory_trt)
x_Pm <- as.numeric(xlabs)-0.1
x_H10 <- as.numeric(xlabs)+0.1

# a) Ambient
plot(xlims,ylims, type="n", xaxt="n", yaxt="n",bty="l",ylab="",xlab="")

y <- Results01$pred %>% filter(Enrichment_trt=="Ambient" & Species=="P. meandrina")
with(y, points(x_Pm, fit,pch=19,cex=1.5,
               col=adjustcolor(cols[2,2],alpha.f=0.8)))	
with(y, segments(x_Pm, lwr95, x_Pm, upr95,
                 col=adjustcolor(cols[2,2],alpha.f=0.8),lwd=2,lend=2))
with(y, segments(x_Pm, lwr85, x_Pm, upr85,
                 col=adjustcolor(cols[2,2],alpha.f=0.8),lwd=4,lend=2))

# P. tuahiniensis
y <- Results01$pred %>% filter(Enrichment_trt=="Ambient" & Species=="P. tuahiniensis")
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
y <- Results01$pred %>% filter(Enrichment_trt=="Enriched" & Species=="P. meandrina")
with(y, points(x_Pm, fit,pch=19,cex=1.5,
               col=adjustcolor(cols[2,2],alpha.f=0.8)))	
with(y, segments(x_Pm, lwr95, x_Pm, upr95,
                 col=adjustcolor(cols[2,2],alpha.f=0.8),lwd=2,lend=2))
with(y, segments(x_Pm, lwr85, x_Pm, upr85,
                 col=adjustcolor(cols[2,2],alpha.f=0.8),lwd=4,lend=2))

# P. tuahiniensis
y <- Results01$pred %>% filter(Enrichment_trt=="Enriched" & Species=="P. tuahiniensis")
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
########################################################################
