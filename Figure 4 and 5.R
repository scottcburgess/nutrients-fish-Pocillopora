### This code produces Figure 4 and Figure 5 and Table S1 in the manuscript: 
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
library('AICcmodavg')
library('car')
library('DHARMa') 
library('emmeans')


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


# T0 = July/Aug 2018
# T3 = July/Aug 2019
# T6 = July/Aug 2020
# T9 = July/Aug 2021
# T10 = November 2021
# Colony size over time.csv has data on the size of corals at each timepoint 
# for the subset of corals that were sampled at T10 for genetics. 
# Used to calculate growth rate

Herb.labels <- c("Very low","Low","Medium","High")

################### Import data ################### 
coral_sizes <- read.csv('Colony size over time.csv')
length(unique(coral_sizes$Un_ID)) # 404 corals
coral_sizesT10 <- read.csv("Recruitment timepoint and size T10.csv")
SpeciesID <- read.csv("Species IDs.csv")

# Coral sizes were measured using ImageJ at timepoints T0, T3, T6, T9 in haplo_coral_sizes_v2.csv
# Coral sizes were measured in situ to the nearest half centimeter at timepoint T10
# in Recruitment timepoint and size T10.csv (November 2021) 
# when corals were sampled for genetic analysis.  




################### Prepare data ################### 
# Add species names 
coral_sizes$Species <- SpeciesID[match(coral_sizes$Un_ID, SpeciesID$Un_ID),which(names(SpeciesID)=="Species")]
coral_sizes$Species <- ifelse(coral_sizes$Species=="", "Unknown",coral_sizes$Species)

coral_sizesT10$Species <- SpeciesID[match(coral_sizesT10$Un_ID, SpeciesID$Un_ID),which(names(SpeciesID)=="Species")]
coral_sizesT10$Species <- ifelse(coral_sizesT10$Species=="", "Unknown",coral_sizesT10$Species)


# Add unique code for each spatial block
coral_sizes$Block <- substr(coral_sizes$Block_plot,1,1)
coral_sizesT10$Block <- substr(coral_sizesT10$Plot,1,1)


# Rename and re-order the herbivory treatments
coral_sizes$Herbivory_trt <- ifelse(coral_sizes$Herbivory_trt=="1X1","Very low",
                                     ifelse(coral_sizes$Herbivory_trt=="2X2","Low",
                                            ifelse(coral_sizes$Herbivory_trt=="3X3","Medium","High")))
coral_sizes$Herbivory_trt <- factor(coral_sizes$Herbivory_trt, 
                                     levels=c("Very low","Low","Medium","High"))

coral_sizesT10$Herbivory_trt <- ifelse(coral_sizesT10$Herbivory_trt=="1X1","Very low",
                                    ifelse(coral_sizesT10$Herbivory_trt=="2X2","Low",
                                           ifelse(coral_sizesT10$Herbivory_trt=="3X3","Medium","High")))
coral_sizesT10$Herbivory_trt <- factor(coral_sizesT10$Herbivory_trt, 
                                    levels=c("Very low","Low","Medium","High"))




# Rearrange to wide format to calculate growth
idcols <- c("Un_ID","Block_plot","Herbivory_trt","Enrichment_trt","Recruitment_timepoint","Species","Block")
coral_sizes_wide <- coral_sizes %>% tidyr::pivot_wider(id_cols=idcols,
                                        names_from=Timepoint,
                                        values_from=Size)

# Get the corals that recruited in T3 = July/Aug 2019 or T6 = July/Aug 2020
coral_sizes_wide <- coral_sizes_wide %>% filter(Recruitment_timepoint %in% c("T3","T6"))

# Make column for 'start' and 'end' size based on recruit size and size after first year
coral_sizes_wide$start_size <- ifelse(coral_sizes_wide$Recruitment_timepoint=="T3",coral_sizes_wide$T3,coral_sizes_wide$T6)
coral_sizes_wide$end_size <- ifelse(coral_sizes_wide$Recruitment_timepoint=="T3",coral_sizes_wide$T6,coral_sizes_wide$T9)

# Calculate log of size
coral_sizes_wide$log_start_size <- log(coral_sizes_wide$start_size)
coral_sizes_wide$log_end_size <- log(coral_sizes_wide$end_size)

# Based on Kelly re-checking the images (Large_recruits_to_check_KES_notes.docx),
# these corals can be removed because it was unclear when they recruited
remove <- c("C2_3X3_P29","C4_1X1_P19","C4_2X2_P27")
coral_sizes_wide <- coral_sizes_wide %>% filter(!(Un_ID %in% remove))
# B2_open_P11 is good to use

# Calculate Relative Growth Rate = RGR
# Size is longest diameter in cm
# RGR is in units of cm cm^-1 yr^-1
coral_sizes_wide$RGR <- (coral_sizes_wide$log_end_size - coral_sizes_wide$log_start_size) / 1
# View(coral_sizes_wide)


# Function to compare models of coral size using AICc 
AnalyzeAICc <- function(dat,yname,family){
  dat <- data.frame(dat)
  dat$yvar <- dat[,names(dat)==yname,]

  m1 <- glmmTMB(yvar ~ Herbivory_trt * Enrichment_trt * Species + (1|Block), family=family, data=dat)
  m2 <- update(m1,.~.-Herbivory_trt:Enrichment_trt:Species)
  m3 <- update(m2,.~.-Herbivory_trt:Enrichment_trt)
  m4 <- update(m2,.~.-Enrichment_trt:Species)
  m5 <- update(m2,.~.-Herbivory_trt:Species)
  m6 <- glmmTMB(yvar ~ Herbivory_trt + Enrichment_trt + Species + Herbivory_trt:Enrichment_trt + (1|Block), family=family, data=dat)
  m7 <- glmmTMB(yvar ~ Herbivory_trt + Enrichment_trt + Species + Enrichment_trt:Species + (1|Block), family=family, data=dat)
  m8 <- glmmTMB(yvar ~ Herbivory_trt + Enrichment_trt + Species + Herbivory_trt:Species + (1|Block), family=family, data=dat)
  m9 <- glmmTMB(yvar ~ Herbivory_trt + Enrichment_trt + Species + (1|Block), family=family, data=dat)
  candidate.models <- list(m1=m1,m2=m2,m3=m3,m4=m4,m5=m5,m6=m6,m7=m7,m8=m8,m9=m9)
  aic_table <- aictab(candidate.models)
  Anova_table <- Anova(get(aic_table$Modnames[1]),type="III")
  
  pred <- expand.grid(Herbivory_trt=unique(dat$Herbivory_trt),
                      Enrichment_trt=unique(dat$Enrichment_trt),
                      Species=unique(dat$Species),
                      Block=NA)
  p <- predict(get(aic_table$Modnames[1]),newdat=pred,se.fit=T)

  crit95 <- qt((1-0.95)/2, df = df.residual(get(aic_table$Modnames[1])), lower.tail = FALSE)
  crit85 <- qt((1-0.85)/2, df = df.residual(get(aic_table$Modnames[1])), lower.tail = FALSE)
  
  if(family$family=='gaussian'){
  pred$fit <- p$fit
  pred$lwr95 <- p$fit - crit95*p$se.fit
  pred$upr95 <- p$fit + crit95*p$se.fit
  pred$lwr85 <- p$fit - crit85*p$se.fit
  pred$upr85 <- p$fit + crit85*p$se.fit}
  
  if(family$family=="Gamma"){
    pred$fit <- exp(p$fit)
    pred$lwr95 <- exp(p$fit - crit95*p$se.fit)
    pred$upr95 <- exp(p$fit + crit95*p$se.fit)
    pred$lwr85 <- exp(p$fit - crit85*p$se.fit)
    pred$upr85 <- exp(p$fit + crit85*p$se.fit)}
  
  DHARMaOutput <- simulateResiduals(fittedModel = get(aic_table$Modnames[1]), plot = F)
  
  list(pred=pred,
       model=get(aic_table$Modnames[1]),
       DHARMaOutput=DHARMaOutput,
       Anova_table=Anova_table,
       aic_table=aic_table)
}


### Analyze coral size data and save results
# T3 
y <- coral_sizes_wide %>% filter(Species %in% c("P. meandrina","P. tuahiniensis"),
                                 Recruitment_timepoint=="T3")
Results_start_size_2019 <- AnalyzeAICc(dat=y,yname="start_size", family=Gamma(link='log'))
# plot(Results_start_size_2019$DHARMaOutput) # looks ok
Results_RGR_2019 <- AnalyzeAICc(dat=y,yname="RGR", family=gaussian())
# plot(Results_RGR_2019$DHARMaOutput) # looks ok

# T6
y <- coral_sizes_wide %>% filter(Species %in% c("P. meandrina","P. tuahiniensis"),
                                 Recruitment_timepoint=="T6")
Results_start_size_2020 <- AnalyzeAICc(dat=y,yname="start_size", family=Gamma(link='log'))
# plot(Results_start_size_2020$DHARMaOutput) # looks ok
Results_RGR_2020 <- AnalyzeAICc(dat=y,yname="RGR", family=gaussian())
# plot(Results_RGR_2020$DHARMaOutput) # looks ok

y <- coral_sizesT10 %>% filter(Species %in% c("P. meandrina","P. tuahiniensis"))
Results_final_size_Nov2021 <- AnalyzeAICc(dat=y,yname="Size_T10", family=Gamma(link='log'))
# plot(Results_final_size_Nov2021$DHARMaOutput) # looks ok



# Function to compare growth using AICc 
AnalyzeAICc <- function(dat){
  dat <- data.frame(dat)
  
  m1 <- suppressWarnings(glmmTMB(log_end_size ~ log_start_size * Herbivory_trt * Enrichment_trt * Species + (1|Block), family=gaussian, data=dat))
  m2 <- update(m1,.~.-Herbivory_trt:Enrichment_trt:Species)
  m3 <- update(m2,.~.-Herbivory_trt:Enrichment_trt)
  m4 <- update(m2,.~.-Enrichment_trt:Species)
  m5 <- update(m2,.~.-Herbivory_trt:Species)
  m6 <- glmmTMB(log_end_size ~   log_start_size + Herbivory_trt + Enrichment_trt + Species +
                  log_start_size:Herbivory_trt + log_start_size:Enrichment_trt + log_start_size:Species +
                  log_start_size:Herbivory_trt:Enrichment_trt + log_start_size:Herbivory_trt:Species + log_start_size:Enrichment_trt:Species + 
                  log_start_size:Herbivory_trt:Enrichment_trt:Species +
                  (1|Block), family=gaussian, data=dat)
  m7 <- update(m6,.~.-log_start_size:Herbivory_trt:Enrichment_trt:Species)
  m8 <- update(m7,.~.-log_start_size:Enrichment_trt:Species)
  m9 <- update(m7,.~.-log_start_size:Herbivory_trt:Species)
  m10 <- update(m7,.~.-log_start_size:Herbivory_trt:Enrichment_trt)
  m11 <- glmmTMB(log_end_size ~   log_start_size + Herbivory_trt + Enrichment_trt + Species + 
                   log_start_size:Herbivory_trt + log_start_size:Enrichment_trt + log_start_size:Species +
                   (1|Block), family=gaussian, data=dat)
  m12 <- update(m11,.~.-log_start_size:Species)
  m13 <- update(m11,.~.-log_start_size:Enrichment_trt)
  m14 <- update(m11,.~.-log_start_size:Herbivory_trt)
  m15 <- glmmTMB(log_end_size ~ log_start_size + Herbivory_trt + Enrichment_trt + Species + (1|Block), family=gaussian, data=dat)
  
  # because there is not enough data to estimate this model for T6
  if(unique(dat$Recruitment_timepoint)[1]=="T3"){
    candidate.models <- list(m1=m1,m2=m2,m3=m3,m4=m4,m5=m5,m6=m6,m7=m7,m8=m8,m9=m9,m10=m10,m11=m11,m12=m12,m13=m13,m14=m14,m15=m15)}

  if(unique(dat$Recruitment_timepoint)[1]=="T6"){
    candidate.models <- list(m2=m2,m3=m3,m4=m4,m5=m5,m6=m6,m7=m7,m8=m8,m9=m9,m10=m10,m11=m11,m12=m12,m13=m13,m14=m14,m15=m15)}
  
  aic_table <- aictab(candidate.models)
  Anova_table <- Anova(get(aic_table$Modnames[1]),type="III")
  
  pred <- expand.grid(log_start_size=seq(min(dat$log_start_size,na.rm=T),
                                         max(dat$log_start_size,na.rm=T),
                                         length.out=100),
                      Herbivory_trt=unique(dat$Herbivory_trt),
                      Enrichment_trt=unique(dat$Enrichment_trt),
                      Species=unique(dat$Species),
                      Block=NA)

  p <- predict(get(aic_table$Modnames[1]),newdat=pred,se.fit=T)
  
  crit95 <- qt((1-0.95)/2, df = df.residual(get(aic_table$Modnames[1])), lower.tail = FALSE)
  crit85 <- qt((1-0.85)/2, df = df.residual(get(aic_table$Modnames[1])), lower.tail = FALSE)
  
  pred$fit <- p$fit
  pred$lwr95 <- p$fit - crit95*p$se.fit
  pred$upr95 <- p$fit + crit95*p$se.fit
  pred$lwr85 <- p$fit - crit85*p$se.fit
  pred$upr85 <- p$fit + crit85*p$se.fit
  
  DHARMaOutput <- simulateResiduals(fittedModel = get(aic_table$Modnames[1]), plot = F)
  
  list(pred=pred,
       model=get(aic_table$Modnames[1]),
       DHARMaOutput=DHARMaOutput,
       Anova_table=Anova_table,
       aic_table=aic_table)
}


### Analyze coral growth data and save results
y <- coral_sizes_wide %>% filter(Species %in% c("P. meandrina","P. tuahiniensis"),
                                 Recruitment_timepoint=="T3")
Results_growth_2019 <- AnalyzeAICc(dat=y)
# plot(Results_growth_2019$DHARMaOutput)

y <- coral_sizes_wide %>% filter(Species %in% c("P. meandrina","P. tuahiniensis"),
                                 Recruitment_timepoint=="T6")
Results_growth_2020 <- AnalyzeAICc(dat=y)
# plot(Results_growth_2020$DHARMaOutput)  

y <- coral_sizes_wide %>% filter(Species %in% c("P. meandrina","P. tuahiniensis"))
Results_growth <- AnalyzeAICc(dat=y)
# plot(Results_growth$DHARMaOutput) # looks ok





################### Results ###################
## Growth
# Results_growth_2019$Anova_table
# Results_growth_2020$Anova_table
Results_growth$Anova_table 
# (χ^2=26.75, df = 1, p < 0.001)
# (χ^2=5.26, df = 3, p = 0.15) 
# (χ^2=5.77, df = 3, p = 0.02). 

# Results_growth_2019$aic_table
# Results_growth_2020$aic_table
Results_growth$aic_table

# emmeans(Results_growth$model,~log_start_size + Species)
EMM <- emmeans(Results_growth$model,~ Species)
pairs(EMM, simple = "Species")
# P. meandrina vs. P. tuahiniensis log(size)
predict(Results_growth$model,
        newdat=expand.grid(Herbivory_trt=c("High"),
                           Enrichment_trt="Ambient",
                           Species=c("P. meandrina","P. tuahiniensis"),
                           log_start_size=c(3),
                           Block=NA),type='link')
# 2.796301 - 2.963112
# log(Pm size) - log(Pt size) = -0.167
# exp(2.796301);exp(2.963112)
(exp(-0.167)-1)*100 # the size of P. meandrina was 15% (10 – 21, 95% CI) smaller than that for P. tuahiniensis   
(exp(-0.167 + c(-2,2)*0.0323)-1)*100


EMM <- emmeans(Results_growth$model,~ Enrichment_trt)
pairs(EMM, simple = "Enrichment_trt")
# Ambient vs Enriched nutrients
predict(Results_growth$model,
        newdat=expand.grid(Herbivory_trt=c("High"),
                           Enrichment_trt=c("Ambient","Enriched"),
                           Species=c("P. meandrina"),
                           log_start_size=c(3),
                           Block=NA),type='link')
# 2.796301 - 2.874584 = -0.078283
(exp(-0.0783)-1)*100 # Sizes were 7% (1 – 13, 95% CI) lower under elevated nutrients compared to ambient conditions 
(exp(-0.0783 + c(-2,2)*0.0326)-1)*100

EMM <- emmeans(Results_growth$model,~ Herbivory_trt)
pairs(EMM, simple = "Herbivory_trt")
Results_start_size_2019$Anova_table
EMM <- emmeans(Results_start_size_2019$model,~ Species)
pairs(EMM, simple = "Species")



# 
Results_start_size_2019$Anova_table # 2019: χ^2=10.27, df = 3, p = 0.02
Results_start_size_2020$Anova_table # 2020: χ^2=9.91, df = 3, p = 0.02
EMM <- emmeans(Results_start_size_2019$model,~ Species)
pairs(EMM, simple = "Species")
(exp(-0.114)-1)*100 
(exp(-0.114 + c(-2,2)*-1.935)-1)*100


Results_final_size_Nov2021$Anova_table
EMM <- emmeans(Results_final_size_Nov2021$model,~ Species)
predict(Results_final_size_Nov2021$model,
        newdat=expand.grid(Herbivory_trt=c("High"),
                           Enrichment_trt=c("Ambient"),
                           Species=c("P. meandrina","P. tuahiniensis"),
                           Block=NA),type='link')

Results_start_size_2019$aic_table
Results_start_size_2020$aic_table
Results_final_size_Nov2021$aic_table


############################################################################




###########################################################################################
# Function to add points and predicted lines to a plot
add_data <- function(dat,predictions,Sp,Herbivory,Enrichment){
  abline(0,1,lty=2,col="grey")
  
  dat <- data.frame(dat)
  
  tmp1 <- dat %>% filter(Species == Sp[1],
                         Herbivory_trt == Herbivory,
                         Enrichment_trt == Enrichment,
                         !is.na(RGR))

  xmin1 <- min(tmp1$log_start_size,ra.rm=T)
  xmax1 <- max(tmp1$log_start_size,ra.rm=T)

  tmp2 <- dat %>% filter(Species == Sp[2],
                         Herbivory_trt == Herbivory,
                         Enrichment_trt == Enrichment,
                         !is.na(RGR))

  xmin2 <- min(tmp2$log_start_size,ra.rm=T)
  xmax2 <- max(tmp2$log_start_size,ra.rm=T)


# add points
  with(tmp1, points(start_size,
                    end_size,
                    pch=19,
                    cex=1,
                    col=adjustcolor(cols[cols$Species==Sp[1],2],alpha.f=0.2)))	

  with(tmp2, points(start_size,
                    end_size,
                    pch=19,
                    cex=1,
                    col=adjustcolor(cols[cols$Species==Sp[2],2],alpha.f=0.2)))	

# add predictions
  y1 <- predictions %>% filter(Species==Sp[1],
                               Herbivory_trt == Herbivory,
                               Enrichment_trt == Enrichment,
                               log_start_size > xmin1 & log_start_size < xmax1)

  with(y1, polygon(c(exp(log_start_size),rev(exp(log_start_size))), c(exp(lwr95),rev(exp(upr95))),
                   col=adjustcolor(cols[cols$Species==Sp[1],2],alpha.f=0.6),lwd=2,border=NA))
  with(y1, lines(exp(log_start_size), exp(fit) ,pch=19,cex=1.5,
                 col=cols[cols$Species==Sp[1],2]))	

  y2 <- predictions %>% filter(Species==Sp[2],Herbivory_trt == Herbivory,
                               Enrichment_trt == Enrichment,
                               log_start_size > xmin2 & log_start_size < xmax2)

  with(y2, polygon(c(exp(log_start_size),rev(exp(log_start_size))), c(exp(lwr95),rev(exp(upr95))),
                   col=adjustcolor(cols[cols$Species==Sp[2],2],alpha.f=0.6),lwd=2,border=NA))
  with(y2, lines(exp(log_start_size), exp(fit) ,pch=19,cex=1.5,
                 col=cols[cols$Species==Sp[2],2]))	

}
###########################################################################################


# set the colors for plotting
cols <- cbind.data.frame(Species = c("P. tuahiniensis",
                                     "P. meandrina"),
                         Color = c("#D55E00",
                                   "#0072B2"))

################## Figure 4 ###############################
quartz(width=5,height=7)
par(mfrow=c(4,2),mar=c(4,2,0,1),oma=c(0,2,4,2))

row_labs <- unique(coral_sizes_wide$Herbivory_trt)

# range(coral_sizes_wide$start_size,na.rm=T)
# range(coral_sizes_wide$end_size,na.rm=T)
xlims <- c(0,6.2)
ylims <- c(0,10.5)

sizes <- seq(0,15,1)

# a) Ambient - Very Low
plot(xlims,ylims, type="n", xaxt="n", yaxt="n",bty="l",ylab="",xlab="")
add_data(dat=coral_sizes_wide,
         predictions=Results_growth$pred,
         Sp=c("P. meandrina","P. tuahiniensis"),
         Herbivory="Very low",
         Enrichment="Ambient")
axis(side=1,sizes)
axis(side=2,sizes,las=1)
mtext(side=3,"Ambient",line=1.5)
mtext(side=3,"a) Very low",adj=0)

# b) Enriched - Very Low
plot(xlims,ylims, type="n", xaxt="n", yaxt="n",bty="l",ylab="",xlab="")
add_data(dat=coral_sizes_wide,
         predictions=Results_growth$pred,
         Sp=c("P. meandrina","P. tuahiniensis"),
         Herbivory="Very low",
         Enrichment="Enriched")
axis(side=1,sizes)
axis(side=2,sizes,las=1)
mtext(side=3,"Enriched",line=1.5)
mtext(side=3,"b) Very low",adj=0)

legend(2.8,3.9,legend=c(
  expression(paste(italic("P. tuahiniensis"))),
  expression(paste(italic("P. meandrina")))),
  lwd=2,
  pt.cex=1,
  cex=1,
  bty="n",
  xpd=T,
  col=c(cols[1,2],
        cols[2,2]))


# c) Ambient - Low
plot(xlims,ylims, type="n", xaxt="n", yaxt="n",bty="l",ylab="",xlab="")
add_data(dat=coral_sizes_wide,
         predictions=Results_growth$pred,
         Sp=c("P. meandrina","P. tuahiniensis"),
         Herbivory="Low",
         Enrichment="Ambient")
axis(side=1,sizes)
axis(side=2,sizes,las=1)
mtext(side=3,"c) Low",adj=0)

# d) Enriched - Low
plot(xlims,ylims, type="n", xaxt="n", yaxt="n",bty="l",ylab="",xlab="")
add_data(dat=coral_sizes_wide,
         predictions=Results_growth$pred,
         Sp=c("P. meandrina","P. tuahiniensis"),
         Herbivory="Low",
         Enrichment="Enriched")
axis(side=1,sizes)
axis(side=2,sizes,las=1)
mtext(side=3,"d) Low",adj=0)

# e) Ambient - Medium
plot(xlims,ylims, type="n", xaxt="n", yaxt="n",bty="l",ylab="",xlab="")
add_data(dat=coral_sizes_wide,
         predictions=Results_growth$pred,
         Sp=c("P. meandrina","P. tuahiniensis"),
         Herbivory="Medium",
         Enrichment="Ambient")
axis(side=1,sizes)
axis(side=2,sizes,las=1)
mtext(side=3,"e) Medium",adj=0)

# f) Enriched - Medium
plot(xlims,ylims, type="n", xaxt="n", yaxt="n",bty="l",ylab="",xlab="")
add_data(dat=coral_sizes_wide,
         predictions=Results_growth$pred,
         Sp=c("P. meandrina","P. tuahiniensis"),
         Herbivory="Medium",
         Enrichment="Enriched")
axis(side=1,sizes)
axis(side=2,sizes,las=1)
mtext(side=3,"f) Medium",adj=0)

# g) Ambient - High
plot(xlims,ylims, type="n", xaxt="n", yaxt="n",bty="l",ylab="",xlab="")
add_data(dat=coral_sizes_wide,
         predictions=Results_growth$pred,
         Sp=c("P. meandrina","P. tuahiniensis"),
         Herbivory="High",
         Enrichment="Ambient")
axis(side=1,sizes)
axis(side=2,sizes,las=1)
mtext(side=3,"g) High",adj=0)

# h) Enriched - High
plot(xlims,ylims, type="n", xaxt="n", yaxt="n",bty="l",ylab="",xlab="")
add_data(dat=coral_sizes_wide,
         predictions=Results_growth$pred,
         Sp=c("P. meandrina","P. tuahiniensis"),
         Herbivory="High",
         Enrichment="Enriched")
axis(side=1,sizes)
axis(side=2,sizes,las=1)
mtext(side=3,"h) High",adj=0)


mtext("Size at recruitment (longest diameter, cm)", side=1, line=-1, cex=1,outer=T)
mtext("Size after 1 year (longest diameter, cm)", side=2.5, line=0, cex=1,outer=T)
########################################################################################### 









###########################################################################################
# Function to add points and predicted lines to a plot
add_data <- function(dat,predictions,Sp,Herbivory,Enrichment){

  dat <- data.frame(dat)
  
  if(names(dat)[6]=="Size_T10"){dat$start_size <- dat$Size_T10}
  
  tmp1 <- dat %>% filter(Species == Sp[1],
                         Enrichment_trt == Enrichment,
                         !is.na(start_size))
  tmp1$x <- as.numeric(tmp1$Herbivory_trt)-0.1
  
  tmp2 <- dat %>% filter(Species == Sp[2],
                         Enrichment_trt == Enrichment,
                         !is.na(start_size))
  tmp2$x <- as.numeric(tmp2$Herbivory_trt)+0.1
  
  # add points
  set.seed(99);with(tmp1, points(jitter(x,0.2),
                    start_size,
                    pch=19,
                    cex=1,
                    col=adjustcolor(cols[cols$Species==Sp[1],2],alpha.f=0.2)))	
  
  set.seed(999);with(tmp2, points(jitter(x,0.2),
                    start_size,
                    pch=19,
                    cex=1,
                    col=adjustcolor(cols[cols$Species==Sp[2],2],alpha.f=0.2)))	
  
  # add predictions
  y1 <- predictions %>% filter(Species==Sp[1],
                               Enrichment_trt == Enrichment)
  
  with(y1, segments(as.numeric(Herbivory_trt)-0.1,
                    lwr95,
                    as.numeric(Herbivory_trt)-0.1,
                    upr95,
                    lend=2,lwd=2,
                    col=adjustcolor(cols[cols$Species==Sp[1],2],alpha.f=0.6)))
  with(y1, segments(as.numeric(Herbivory_trt)-0.1,
                    lwr85,
                    as.numeric(Herbivory_trt)-0.1,
                    upr85,
                    lend=2,lwd=4,
                    col=adjustcolor(cols[cols$Species==Sp[1],2],alpha.f=0.6)))
  with(y1, points(as.numeric(Herbivory_trt)-0.1,fit, pch=19,cex=2,
                  col=adjustcolor(cols[cols$Species==Sp[1],2],alpha.f=0.6)))	
  
  y2 <- predictions %>% filter(Species==Sp[2],
                               Enrichment_trt == Enrichment)
  with(y2, segments(as.numeric(Herbivory_trt)+0.1,
                    lwr95,
                    as.numeric(Herbivory_trt)+0.1,
                    upr95,
                    lend=2,lwd=2,
                    col=adjustcolor(cols[cols$Species==Sp[2],2],alpha.f=0.6)))
  with(y2, segments(as.numeric(Herbivory_trt)+0.1,
                    lwr85,
                    as.numeric(Herbivory_trt)+0.1,
                    upr85,
                    lend=2,lwd=4,
                    col=adjustcolor(cols[cols$Species==Sp[2],2],alpha.f=0.6)))
  with(y2, points(as.numeric(Herbivory_trt)+0.1,fit, pch=19,cex=2,
                  col=adjustcolor(cols[cols$Species==Sp[2],2],alpha.f=0.6)))	
}
########################################################################################### 



################## Figure 5 ###############################
quartz(width=5,height=6)
par(mfrow=c(3,2),mar=c(4,2,0,1),oma=c(3,2,3,2))

xlabs <- unique(coral_sizes_wide$Herbivory_trt)

# range(coral_sizes_wide$start_size,na.rm=T)
# range(coral_sizes_wide$end_size,na.rm=T)
# range(coral_sizesT10$Size_T10,na.rm=T)
xlims <- c(0.8,4.2)
ylims1 <- c(0,6.5)
ylims2 <- c(0,20)

sizes <- seq(0,15,1)
sizes2 <- seq(0,30,2)

# a) Ambient 2019
plot(xlims,ylims1, type="n", xaxt="n", yaxt="n",bty="l",ylab="",xlab="")
add_data(dat=coral_sizes_wide,
         predictions=Results_start_size_2019$pred,
         Sp=c("P. meandrina","P. tuahiniensis"),
         Enrichment="Ambient")
axis(side=1,at=1:4,labels=NA)
axis(side=2,sizes,las=1,cex=1.2)
mtext(side=2,"Size of recruits (2019)",line=2.5,adj=0.2)
mtext(side=3,"a) Ambient",adj=0)

# b) Enriched 2019
plot(xlims,ylims1, type="n", xaxt="n", yaxt="n",bty="l",ylab="",xlab="")
add_data(dat=coral_sizes_wide,
         predictions=Results_start_size_2019$pred,
         Sp=c("P. meandrina","P. tuahiniensis"),
         Enrichment="Enriched")
axis(side=1,at=1:4,labels=NA)
axis(side=2,sizes,las=1,cex=1.2)
mtext(side=3,"b) Enriched",adj=0)

# c) Ambient 2020
plot(xlims,ylims1, type="n", xaxt="n", yaxt="n",bty="l",ylab="",xlab="")
add_data(dat=coral_sizes_wide,
         predictions=Results_start_size_2020$pred,
         Sp=c("P. meandrina","P. tuahiniensis"),
         Enrichment="Ambient")
axis(side=1,at=1:4,labels=NA)
axis(side=2,sizes,las=1,cex=1.2)
mtext(side=2,"Size of recruits (2020)",line=2.5,adj=0)
mtext(side=3,"c) Ambient",adj=0)

# d) Enriched 2020
plot(xlims,ylims1, type="n", xaxt="n", yaxt="n",bty="l",ylab="",xlab="")
add_data(dat=coral_sizes_wide,
         predictions=Results_start_size_2020$pred,
         Sp=c("P. meandrina","P. tuahiniensis"),
         Enrichment="Enriched")
axis(side=1,at=1:4,labels=NA)
axis(side=2,sizes,las=1,cex=1.2)
mtext(side=3,"d) Enriched",adj=0)

# e) Ambient Nov 2021
plot(xlims,ylims2, type="n", xaxt="n", yaxt="n",bty="l",ylab="",xlab="")
add_data(dat=coral_sizesT10,
         predictions=Results_final_size_Nov2021$pred,
         Sp=c("P. meandrina","P. tuahiniensis"),
         Enrichment="Ambient")
axis(side=1,at=1:4,labels=NA)
axis(side=2,sizes2,las=1,cex=1.2)
mtext(side=2,"Size of all colonies (2021)",line=2.5,adj=0.4)
mtext(side=3,"e) Ambient",adj=0)
text(1:4,rep(-2.5,4),xlabs,xpd=T,srt=45,adj=1,cex=1.2)

# f) Enriched Nov 2021
plot(xlims,ylims2, type="n", xaxt="n", yaxt="n",bty="l",ylab="",xlab="")
add_data(dat=coral_sizesT10,
         predictions=Results_final_size_Nov2021$pred,
         Sp=c("P. meandrina","P. tuahiniensis"),
         Enrichment="Enriched")
axis(side=1,at=1:4,labels=NA)
axis(side=2,sizes2,las=1,cex=1.2)
mtext(side=3,"f) Enriched",adj=0)
text(1:4,rep(-2.5,4),xlabs,xpd=T,srt=45,adj=1,cex=1.2)

legend(2,22,legend=c(
  expression(paste(italic("P. tuahiniensis"))),
  expression(paste(italic("P. meandrina")))),
  pch=19,
  lwd=2,
  pt.cex=1.2,
  cex=1,
  bty="n",
  xpd=T,
  col=c(adjustcolor(cols[1,2],alpha.f = 0.6),
        adjustcolor(cols[2,2],alpha.f = 0.6)))

mtext("Consumer pressure treatment", side=1, line=1, cex=1.2, outer=T)
########################################################################################### 




