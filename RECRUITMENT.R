# Code to analyse results from the 'RECHARGE' plots
# Code written by Scott Burgess July 2023

# R version 4.2.3 (2021-03-31) -- "Shortstop Beagle"
library('dplyr') # v1.2.0
library('tidyr') # v1.3.0
library('glmmTMB') # v1.1.3
library('AICcmodavg')
library('car')
library('emmeans')
# library('bbmle') # 1.0.25
library('DHARMa') # 0.4.6
library('ggplot2')
library('patchwork') # 1.1.2
# library(broom.mixed)




# set the colors for plotting
cols <- cbind.data.frame(Species = c("P. tuahiniensis",
                                     "P. meandrina",
                                     "P. grandis",
                                     "P. verrucosa",
                                     "P. cf. effusa"),
                         Color = c("#D55E00",
                                   "#0072B2",
                                   "#56B4E9",
                                   "#E69F00",
                                   "#009E73"))




# Import and prepare data
algae <- read.csv('plot_macroalgae_summary_to_T10.csv')

dat <- read.csv("sampled_corals_for_haplo_v3.csv")
SpeciesID <- read.csv("Species IDs.csv")


# Rename the herbivory treatments
dat$Herbivory_trt <- ifelse(dat$Herbivory_trt=="1X1","Very low",
                            ifelse(dat$Herbivory_trt=="2X2","Low",
                                   ifelse(dat$Herbivory_trt=="3X3","Medium","High")))
dat$Herbivory_trt <- factor(dat$Herbivory_trt, 
                            levels=c("Very low","Low","Medium","High"))

algae$Herbivory_trt <- ifelse(algae$Herbivory_trt=="1X1","Very low",
                              ifelse(algae$Herbivory_trt=="2X2","Low",
                                     ifelse(algae$Herbivory_trt=="3X3","Medium","High")))
algae$Herbivory_trt <- factor(algae$Herbivory_trt, 
                              levels=c("Very low","Low","Medium","High"))

# Change Plot to A, B, C, or D
dat$Block <- substr(dat$Plot,1,1)
algae$Block <- substr(algae$Block_Plot,1,1)

# Add species ID's to dat
unique(dat$Un_ID)
unique(SpeciesID$Un_ID)
dat <- cbind.data.frame(dat,
                 SpeciesID[match(dat$Un_ID,SpeciesID$Un_ID),11:16] # check column indices are correct
)

# these should all be removed (were in spreadsheet but no sample provided to us, or bad sample, see notes)
dat %>% filter(flag_missing_sample==1) %>% 
  select(Un_ID,Plot,Herbivory_trt,Enrichment_trt,Coral_ID_past, Coral_ID_field, Coral_ID_corrected, Coral_ID, Size_T10, Recruitment_timepoint,Important_notes)
dat <- dat %>% filter(flag_missing_sample!=1)

# could not get sequence from these corals, the only sample without a species ID
dat %>% filter(Species=="")
# Call the blank species names "Unknown"
dat$Species <- ifelse(dat$Species=="", "Unknown",dat$Species)
dat %>% filter(Species=="Unknown")

# only 5 corals flagged as partial mortality
dat %>% filter(flag_partial_mortality==1) 




# How many Pocillopora in total?
nrow(dat) # n=406
nrow(dat) / 32 / 1.25 # 10.15m2 sampling density; 32 cages, each 1.25m2 according to Tom

# How many were not identified to species? n=1 
dat %>% filter(Species=="Unknown") %>% select(Species)

# What proportion were each species - RESULTS
dat %>% group_by(Species) %>% summarise(count=n()) %>% mutate(freq=(count/sum(count)*100))


# Calculate the number of colonies of each species
# First, remove the individuals with uncertain recruit timing
tmp <- dat %>% filter(flag_uncertain_recruit_time==0) 
# number of colonies of each species at each time point
counts_by_timepoint <- aggregate(Un_ID ~ Species + Recruitment_timepoint,
                       FUN=length,
                       data=tmp,
                       drop=F)
names(counts_by_timepoint)[which(names(counts_by_timepoint)=="Un_ID")] <- "Freq"
counts_by_timepoint$Freq[is.na(counts_by_timepoint$Freq)] <- 0

# number of colonies of each species overall
counts_by_all <- aggregate(Un_ID ~ Species,
                       FUN=length,
                       data=dat,
                       drop=F)
names(counts_by_all)[which(names(counts_by_all)=="Un_ID")] <- "Freq"
counts_by_all$Freq[is.na(counts_by_all$Freq)] <- 0


# control the order of Species
counts_by_timepoint$Species <- factor(counts_by_timepoint$Species, levels=c("P. tuahiniensis",
                                                        "P. meandrina",
                                                        "P. grandis",
                                                        "P. verrucosa",
                                                        "P. cf. effusa",
                                                        "Unknown"))
counts_by_all$Species <- factor(counts_by_all$Species, levels=levels(counts_by_timepoint$Species))

# Remove the one colony that could not be genetically identified
counts_by_timepoint <- counts_by_timepoint %>% filter(Species!="Unknown") 
counts_by_all <- counts_by_all %>% filter(Species!="Unknown") 

######## FIGURE 2
quartz(width=4,height=7)
# par(mfrow=c(5,1),mar=c(1,0,1.5,1),oma=c(5,4,1,0))
# 2018
y <- counts_by_timepoint %>% filter(Recruitment_timepoint=="T0")
y <- y[match(cols$Species,y$Species),]
p1 <- ggplot(y, aes(fill=Species, y=Freq, x=Species)) + 
  ylim(0,90) +
  geom_bar(position='stack', stat='identity') +
  theme_minimal() + 
  labs(x='', y='', title='a) Recruited in 2018') +
  theme(plot.title = element_text(hjust=0, size=10, face='bold'),
        axis.text.x = element_blank(),
        legend.position = "none") +
  scale_fill_manual('Species', values=cols$Color)

# 2019
y <- counts_by_timepoint %>% filter(Recruitment_timepoint=="T3")
y <- y[match(cols$Species,y$Species),]
p2 <- ggplot(y, aes(fill=Species, y=Freq, x=Species)) + 
  ylim(0,90) +
  geom_bar(position='stack', stat='identity') +
  theme_minimal() + 
  labs(x='', y='', title='b) Recruited in 2019') +
  theme(plot.title = element_text(hjust=0, size=10, face='bold'),
        axis.text.x = element_blank(),
        legend.position = "none") +
  scale_fill_manual('Species', values=cols$Color)

y_test <- aggregate(Freq ~ Species, FUN=sum, data=y)
y_test <- y_test %>% filter(Species %in% c("P. tuahiniensis","P. meandrina"))
CT_2019 <- chisq.test(as.table(y_test$Freq))

# 2020
y <- counts_by_timepoint %>% filter(Recruitment_timepoint=="T6")
y <- y[match(cols$Species,y$Species),]
p3 <- ggplot(y, aes(fill=Species, y=Freq, x=Species)) + 
  ylim(0,90) +
  geom_bar(position='stack', stat='identity') +
  theme_minimal() + 
  labs(x='', y='Number of colonies in 2021', title='c) Recruited in 2020') +
  theme(plot.title = element_text(hjust=0, size=10, face='bold'),
        axis.text.x = element_blank(),
        legend.position = "none") +
  scale_fill_manual('Species', values=cols$Color)

y_test <- aggregate(Freq ~ Species, FUN=sum, data=y)
y_test <- y_test %>% filter(Species %in% c("P. tuahiniensis","P. meandrina"))
CT_2020 <- chisq.test(as.table(y_test$Freq))

# 2021
y <- counts_by_timepoint %>% filter(Recruitment_timepoint=="T9")
y <- y[match(cols$Species,y$Species),]
p4 <- ggplot(y, aes(fill=Species, y=Freq, x=Species)) + 
  ylim(0,90) +
  geom_bar(position='stack', stat='identity') +
  theme_minimal() + 
  labs(x='', y='', title='d) Recruited in 2021') +
  theme(plot.title = element_text(hjust=0, size=10, face='bold'),
        axis.text.x = element_blank(),
        legend.position = "none") +
  scale_fill_manual('Species', values=cols$Color)

y_test <- aggregate(Freq ~ Species, FUN=sum, data=y)
y_test <- y_test %>% filter(Species %in% c("P. tuahiniensis","P. meandrina"))
CT_2021 <- chisq.test(as.table(y_test$Freq))

# Overall
y <- counts_by_all
y <- y[match(cols$Species,y$Species),]
p5 <- ggplot(y, aes(fill=Species, y=Freq, x=Species)) + 
  ylim(0,200) +
  geom_bar(position='stack', stat='identity') +
  theme_minimal() + 
  labs(x='', y='', title='e) Overall') +
  theme(plot.title = element_text(hjust=0, size=10, face='bold'),
        axis.text.x = element_text(angle = 45, hjust=1.1,vjust=1.1),
        legend.position = "none") +
  scale_fill_manual('Species', values=cols$Color) + geom_blank()

p1/p2/p3/p4/p5  

y_test <- aggregate(Freq ~ Species, FUN=sum, data=y)
y_test <- y_test %>% filter(Species %in% c("P. tuahiniensis","P. meandrina"))
CT_Overall <- chisq.test(as.table(y_test$Freq))

CT_2019
CT_2020
CT_2021
CT_Overall



# Calculate the number of colonies of each species in each plot and recruit time point
counts_by <- aggregate(Un_ID ~ Species + Herbivory_trt + Enrichment_trt + Plot,
                       FUN=length,
                       data=dat,
                       drop=F)
# Drop the false unwanted treatment combinations created by previous operation
# with(dat,table(Plot,Enrichment_trt))
counts_by <- counts_by %>% filter((Plot=="A3" & Enrichment_trt=="Enriched") |
                       (Plot=="A4" & Enrichment_trt=="Ambient") |
                       (Plot=="B1" & Enrichment_trt=="Enriched") |
                       (Plot=="B2" & Enrichment_trt=="Ambient") |
                       (Plot=="C2" & Enrichment_trt=="Ambient") |
                       (Plot=="C4" & Enrichment_trt=="Enriched") |
                       (Plot=="D1" & Enrichment_trt=="Ambient") |
                       (Plot=="D3" & Enrichment_trt=="Enriched"))
names(counts_by)[which(names(counts_by)=="Un_ID")] <- "Freq"
counts_by$Freq[is.na(counts_by$Freq)] <- 0
# View(counts_by %>% filter(Plot=="A3"))

# Add the 'Block' variable bank 
counts_by$Block <- substr(counts_by$Plot,1,1)


# Function to compare models using AICc 
AnalyzeAICc <- function(dat,yname){
  dat$yvar <- dat[,names(dat)==yname,]
  
  m1 <- glmmTMB(yvar ~ Herbivory_trt * Enrichment_trt * Species + (1|Block), family="poisson", data=dat)
  m2 <- update(m1,.~.-Herbivory_trt:Enrichment_trt:Species)
  m3 <- update(m2,.~.-Herbivory_trt:Enrichment_trt)
  m4 <- update(m2,.~.-Enrichment_trt:Species)
  m5 <- update(m2,.~.-Herbivory_trt:Species)
  m6 <- glmmTMB(yvar ~ Herbivory_trt + Enrichment_trt + Species + Herbivory_trt:Enrichment_trt + (1|Block), family="poisson", data=dat)
  m7 <- glmmTMB(yvar ~ Herbivory_trt + Enrichment_trt + Species + Enrichment_trt:Species + (1|Block), family="poisson", data=dat)
  m8 <- glmmTMB(yvar ~ Herbivory_trt + Enrichment_trt + Species + Herbivory_trt:Species + (1|Block), family="poisson", data=dat)
  m9 <- glmmTMB(yvar ~ Herbivory_trt + Enrichment_trt + Species + (1|Block), family="poisson", data=dat)
  m10 <- glmmTMB(yvar ~ Herbivory_trt + Enrichment_trt + Herbivory_trt:Enrichment_trt + (1|Block), family="poisson", data=dat)
  m11 <- glmmTMB(yvar ~ Herbivory_trt + Enrichment_trt + (1|Block), family="poisson", data=dat)
  candidate.models <- list(m1=m1,m2=m2,m3=m3,m4=m4,m5=m5,m6=m6,m7=m7,m8=m8,m9=m9,m10=m10,m11=m11)
  aic_table <- aictab(candidate.models)
  Anova_table <- Anova(get(aic_table$Modnames[1]),type="III")
  
  pred <- expand.grid(Herbivory_trt=unique(dat$Herbivory_trt),
                      Enrichment_trt=unique(dat$Enrichment_trt),
                      Species=unique(dat$Species),
                      Block=NA)
  p <- predict(get(aic_table$Modnames[1]),newdat=pred,se.fit=T)

  crit95 <- qt((1-0.95)/2, df = df.residual(get(aic_table$Modnames[1])), lower.tail = FALSE)
  crit85 <- qt((1-0.85)/2, df = df.residual(get(aic_table$Modnames[1])), lower.tail = FALSE)
  
  pred$fit <- exp(p$fit)
  pred$lwr95 <- exp(p$fit - crit95*p$se.fit)
  pred$upr95 <- exp(p$fit + crit95*p$se.fit)
  pred$lwr85 <- exp(p$fit - crit85*p$se.fit)
  pred$upr85 <- exp(p$fit + crit85*p$se.fit)
  
  DHARMaOutput <- simulateResiduals(fittedModel = get(aic_table$Modnames[1]), plot = F)
  
  list(pred=pred,
       model=get(aic_table$Modnames[1]),
       DHARMaOutput=DHARMaOutput,
       Anova_table=Anova_table,
       aic_table=aic_table)
}

# Select the two most common species to compare
y <- counts_by %>% filter(Species %in% c("P. meandrina","P. tuahiniensis"))
# rough expectations:
# aggregate(Freq ~ Herbivory_trt * Enrichment_trt * Species,FUN=mean,data=y,drop=F) 

Results <- AnalyzeAICc(dat=y,yname="Freq")

# plot(Results$DHARMaOutput) # looks ok

Results$aic_table
Results$Anova_table

# Results$model
EMM <- emmeans(Results$model,~Herbivory_trt:Species)
### Simple pairwise comparisons
# pairs(EMM, simple = "Species")     
pairs(EMM, simple = "Herbivory_trt")
# test(pairs(EMM, by = "Species"), by = NULL, adjust = "mvt")
(exp(0.5480)-1)*100 # increase of 73% (14 - 162, 95% CI), P. tua
(exp(0.5480 + c(-2,2)*0.207)-1)*100
#vs.
(exp(-0.2877)-1)*100 # reduction of 25% (-52 - 16,  95% CI), P. mean
(exp(-0.2877 + c(-2,2)*0.220)-1)*100

(exp(-0.6763)-1)*100 # reduction of 49% (20 - 68, 95% CI), P. mea
(exp(-0.6763 + c(-2,2)*0.224)-1)*100

EMM <- emmeans(Results$model,~Enrichment_trt)
pairs(EMM, simple = "Enrichment_trt")
((1/exp(-0.268))-1)*100 # nutrients increases recruitment by 30% (6 - 62%, 95% CI)
(1/exp(-0.268 + c(-2,2)*0.106)-1)*100





# Function to compare models using AICc for the **ALGAE data
AnalyzeAICc_algae <- function(dat,yname){
  
  dat$yvar <- dat[,names(dat)==yname,]
  
  m1 <- glmmTMB(yvar ~ Herbivory_trt * Enrichment_trt + (1|Block), family=beta_family(link = "logit"), data=dat)
  m2 <- glmmTMB(yvar ~ Herbivory_trt + Enrichment_trt + (1|Block), family=beta_family(link = "logit"), data=dat)
  candidate.models <- list(m1=m1,m2=m2)
  aic_table <- aictab(candidate.models)
  Anova_table <- Anova(get(aic_table$Modnames[1]),type="III")
  
  pred <- expand.grid(Herbivory_trt=unique(dat$Herbivory_trt),
                      Enrichment_trt=unique(dat$Enrichment_trt),
                      Block=NA)
  p <- predict(get(aic_table$Modnames[1]),newdat=pred,se.fit=T)
  
  crit95 <- qt((1-0.95)/2, df = df.residual(get(aic_table$Modnames[1])), lower.tail = FALSE)
  crit85 <- qt((1-0.85)/2, df = df.residual(get(aic_table$Modnames[1])), lower.tail = FALSE)
  
  pred$fit <- plogis(p$fit)
  pred$lwr95 <- plogis(p$fit - crit95*p$se.fit)
  pred$upr95 <- plogis(p$fit + crit95*p$se.fit)
  pred$lwr85 <- plogis(p$fit - crit85*p$se.fit)
  pred$upr85 <- plogis(p$fit + crit85*p$se.fit)
  
  DHARMaOutput <- simulateResiduals(fittedModel = get(aic_table$Modnames[1]), plot = F)
  
  list(pred=pred,
       model=get(aic_table$Modnames[1]),
       DHARMaOutput=DHARMaOutput,
       Anova_table=Anova_table,
       aic_table=aic_table)
}

# convert the percentage data to proportions to fit Gamma GLM.
algae$mean_macroalgae_prop <- algae$mean_macroalgae/100
algae$median_macroalgae_prop <- algae$median_macroalgae/100

Results_algae_mean <- AnalyzeAICc_algae(dat=algae, yname='mean_macroalgae_prop')
Results_algae_median <- AnalyzeAICc_algae(dat=algae, yname='median_macroalgae_prop')
# plot(Results_algae_mean$DHARMaOutput) # looks ok
# plot(Results_algae_median$DHARMaOutput) # looks ok

Results_algae_mean$aic_table
Results_algae_median$aic_table

Results_algae_mean$Anova_table
Results_algae_median$Anova_table


# Results$model
EMM <- emmeans(Results_algae_mean$model,~Herbivory_trt)
pairs(EMM, simple = "Herbivory_trt")

EMM <- emmeans(Results_algae_mean$model,~Enrichment_trt)
pairs(EMM, simple = "Enrichment_trt")






################################ FIGURE 3
quartz(width=6,height=5)
par(mfrow=c(2,2),mar=c(5,2,0,1),oma=c(2,4,2,0))
ylims <- c(0,17) #30
ylims_algae <- c(0,0.4) 
xlims <- c(0.8,4.2)
xlabs <- unique(Results$pred$Herbivory_trt)
x_Pm <- as.numeric(xlabs)-0.1
x_ma <- as.numeric(xlabs)
x_H10 <- as.numeric(xlabs)+0.1

# a) Ambient
plot(xlims,ylims, type="n", xaxt="n", yaxt="n",bty="l",ylab="",xlab="")

legend('topleft',legend=c(
  expression(paste(italic("P. meandrina"))),
  expression(paste(italic("P. tuahiniensis")))),
       pch=19,
       pt.cex=1.2,
       bty="n",
       xpd=T,
       col=c(adjustcolor(cols[2,2],alpha.f=0.8),
             adjustcolor(cols[1,2],alpha.f=0.8)))

y <- counts_by %>% filter(Species %in% c("P. meandrina"),
                          Enrichment_trt %in% c("Ambient"))
set.seed(99); x <- as.numeric(factor(y$Herbivory_trt))-0.1
with(y, points(jitter(x,0.2), 
               Freq, 
               pch=1,
               cex=1,
               col=adjustcolor(cols[2,2],alpha.f=0.4)))	

y <- Results$pred %>% filter(Enrichment_trt=="Ambient" & Species=="P. meandrina")
with(y, points(x_Pm, fit,pch=19,cex=1.5,
               col=adjustcolor(cols[2,2],alpha.f=0.8)))	
with(y, segments(x_Pm, lwr95, x_Pm, upr95,
                 col=adjustcolor(cols[2,2],alpha.f=0.8),lwd=2,lend=2))
with(y, segments(x_Pm, lwr85, x_Pm, upr85,
                 col=adjustcolor(cols[2,2],alpha.f=0.8),lwd=4,lend=2))
# Haplotype 10
y <- counts_by %>% filter(Species %in% c("P. tuahiniensis"),
                          Enrichment_trt %in% c("Ambient"),
                          !is.na(Freq))
set.seed(99); x <- as.numeric(factor(y$Herbivory_trt))+0.1
with(y, points(jitter(x,0.2), 
               Freq, 
               pch=1,
               cex=1,
               col=adjustcolor(cols[1,2],alpha.f=0.4)))	

y <- Results$pred %>% filter(Enrichment_trt=="Ambient" & Species=="P. tuahiniensis")
with(y, points(x_H10, fit,pch=19,cex=1.5,
               col=adjustcolor(cols[1,2],alpha.f=0.8)))	
with(y, segments(x_H10, lwr95, x_H10, upr95,
                 col=adjustcolor(cols[1,2],alpha.f=0.8),lwd=2,lend=2))
with(y, segments(x_H10, lwr85, x_H10, upr85,
                 col=adjustcolor(cols[1,2],alpha.f=0.8),lwd=4,lend=2))
axis(side=1,at=1:4,c("","","",""))
text(1:4,rep(-3,4),xlabs,xpd=T,srt=45,adj=1)
axis(side=2,at=seq(0,30,2),las=1,cex.axis=1.2)
mtext("a) Ambient",side=3,line=0,adj=0)

mtext("Number of colonies", side=2, line=3.7, adj=0.5,cex=1.2)
mtext("per plot", side=2, line=2.5, adj=0.5,cex=1.2)


# b) Enrichment
plot(xlims,ylims, type="n", xaxt="n", yaxt="n",bty="l",ylab="",xlab="")
# P.meandrina
y <- counts_by %>% filter(Species %in% c("P. meandrina"),
                          Enrichment_trt %in% c("Enriched"),
                          !is.na(Freq))
set.seed(99); x <- as.numeric(factor(y$Herbivory_trt))-0.1
with(y, points(jitter(x,0.2), 
               Freq, 
               pch=1,
               cex=1,
               col=adjustcolor(cols[2,2],alpha.f=0.4)))	

y <- Results$pred %>% filter(Enrichment_trt=="Enriched" & Species=="P. meandrina")
with(y, points(x_Pm, fit,pch=19,cex=1.5,
               col=adjustcolor(cols[2,2],alpha.f=0.8)))	
with(y, segments(x_Pm, lwr95, x_Pm, upr95,
                 col=adjustcolor(cols[2,2],alpha.f=0.8),lwd=2,lend=2))
with(y, segments(x_Pm, lwr85, x_Pm, upr85,
                 col=adjustcolor(cols[2,2],alpha.f=0.8),lwd=4,lend=2))

# P. tuahiniensis
y <- counts_by %>% filter(Species %in% c("P. tuahiniensis"),
                          Enrichment_trt %in% c("Enriched"),
                          !is.na(Freq))
set.seed(99); x <- as.numeric(factor(y$Herbivory_trt))+0.1
with(y, points(jitter(x,0.2), 
               Freq, 
               pch=1,
               cex=1,
               col=adjustcolor(cols[1,2],alpha.f=0.4)))	

y <- Results$pred %>% filter(Enrichment_trt=="Enriched" & Species=="P. tuahiniensis")
with(y, points(x_H10, fit,pch=19,cex=1.5,
               col=adjustcolor(cols[1,2],alpha.f=0.8)))	
with(y, segments(x_H10, lwr95, x_H10, upr95,
                 col=adjustcolor(cols[1,2],alpha.f=0.8),lwd=2,lend=2))
with(y, segments(x_H10, lwr85, x_H10, upr85,
                 col=adjustcolor(cols[1,2],alpha.f=0.8),lwd=4,lend=2))
axis(side=1,at=1:4,c("","","",""))
text(1:4,rep(-3,4),xlabs,xpd=T,srt=45,adj=1)
axis(side=2,at=seq(0,30,2),las=1,cex.axis=1.2)
mtext("b) Nutrient Enrichment",side=3,line=0,adj=0)



# c) Ambient - ALGAE
plot(xlims,ylims_algae, type="n", xaxt="n", yaxt="n",bty="l",ylab="",xlab="")

y <- algae %>% filter(Enrichment_trt %in% c("Ambient"))
set.seed(99); x <- as.numeric(factor(y$Herbivory_trt))
with(y, points(jitter(x,0.2), 
               mean_macroalgae_prop, 
               pch=1,
               cex=1,
               col=adjustcolor('black',alpha.f=0.4)))	

y <- Results_algae_mean$pred %>% filter(Enrichment_trt=="Ambient")
with(y, points(x_ma, fit,pch=19,cex=1.5,
               col=adjustcolor('black',alpha.f=0.8)))	
with(y, segments(x_ma, lwr95, x_ma, upr95,
                 col=adjustcolor('black',alpha.f=0.8),lwd=2,lend=2))
with(y, segments(x_ma, lwr85, x_ma, upr85,
                 col=adjustcolor('black',alpha.f=0.8),lwd=4,lend=2))

axis(side=1,at=1:4,c("","","",""))
text(1:4,rep(-0.075,4),xlabs,xpd=T,srt=45,adj=1)
axis(side=2,at=seq(0,1,0.1),las=1,cex.axis=1.2)
mtext("c) Ambient",side=3,line=0,adj=0)

mtext(expression(paste("Proportion cover of")), side=2, line=4, adj=0.5, cex=1.2)
mtext(expression(paste("macroalgae per plot")), side=2, line=2.5, adj=0.5, cex=1.2)


# d) Enrichment
plot(xlims,ylims_algae, type="n", xaxt="n", yaxt="n",bty="l",ylab="",xlab="")


y <- algae %>% filter(Enrichment_trt %in% c("Enriched"))
set.seed(99); x <- as.numeric(factor(y$Herbivory_trt))
with(y, points(jitter(x,0.2), 
               mean_macroalgae_prop, 
               pch=1,
               cex=1,
               col=adjustcolor('black',alpha.f=0.4)))	

y <- Results_algae_mean$pred %>% filter(Enrichment_trt=="Enriched")
with(y, points(x_ma, fit,pch=19,cex=1.5,
               col=adjustcolor('black',alpha.f=0.8)))	
with(y, segments(x_ma, lwr95, x_ma, upr95,
                 col=adjustcolor('black',alpha.f=0.8),lwd=2,lend=2))
with(y, segments(x_ma, lwr85, x_ma, upr85,
                 col=adjustcolor('black',alpha.f=0.8),lwd=4,lend=2))

axis(side=1,at=1:4,c("","","",""))
text(1:4,rep(-0.075,4),xlabs,xpd=T,srt=45,adj=1)
axis(side=2,at=seq(0,1,0.1),las=1,cex.axis=1.2)
mtext("d) Nutrient Enrichment",side=3,line=0,adj=0)

mtext("Fish consumer pressure treatment", side=1, line=0, cex=1.2, outer=T)
########################################################################


# y <- counts_by %>% filter(Species %in% c("P. meandrina"))
# with(y, interaction.plot(Enrichment_trt,Herbivory_trt,mean_macroalgae,col=Herbivory_trt,lwd=3,lty=1))
# with(y, interaction.plot(Enrichment_trt,Herbivory_trt,median_macroalgae,col=Herbivory_trt,lwd=3,lty=1))




