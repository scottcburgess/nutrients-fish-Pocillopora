# Code to analyse results from the 'RECHARGE' plots
# Code written by Scott Burgess July 2022

# R version 4.2.3 (2023-03-15) -- "Shortstop Beagle"
library('dplyr') # v1.1.2
library('glmmTMB') # v1.1.7
library('AICcmodavg')
library('car')
library('emmeans')
library('DHARMa') # 0.4.6


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



# Function to compare models using AICc 
AnalyzeAICc <- function(dat,yname,family){
  dat$yvar <- dat[,names(dat)==yname,]
  
  m1 <- glmmTMB(yvar ~ Herbivory_trt * Enrichment_trt + (1|Block), family=family, data=dat)
  m2 <- glmmTMB(yvar ~ Herbivory_trt + Enrichment_trt + (1|Block), family=family, data=dat)

  aic_table <- aictab(list(m1=m1,m2=m2))
  Anova_table <- Anova(get(aic_table$Modnames[1]),type="III")
  
  pred <- expand.grid(Herbivory_trt=unique(dat$Herbivory_trt),
                      Enrichment_trt=unique(dat$Enrichment_trt),
                      Block=NA)
  p <- predict(get(aic_table$Modnames[1]),newdat=pred,se.fit=T)
  pred$fit <- (p$fit)
  pred$lwr <- (p$fit - 2*p$se.fit)
  pred$upr <- (p$fit + 2*p$se.fit)
  
  DHARMaOutput <- simulateResiduals(fittedModel = get(aic_table$Modnames[1]), plot = F)
  
  list(pred=pred,
       model=get(aic_table$Modnames[1]),
       DHARMaOutput=DHARMaOutput,
       Anova_table=Anova_table,
       aic_table=aic_table)
}


# Annual survival of Pocillopora in the plots from 2018 to 2021
Times <- cbind.data.frame(time.point=c(c("T0","T3","T6","T9")),
                          Year = c("2018","2019","2020","2021"))
# T0 = July/Aug 2018
# T3 = July/Aug 2019
# T6 = July/Aug 2020
# T9 = July/Aug 2021
# T10 = November 2021

Herb.labels <- c("Very low","Low","Medium","High")

# Import data
poc_recruit_mort <- read.csv('poc_recruit_mort.csv')

# Add a blocking factor
poc_recruit_mort$Block <- substr(poc_recruit_mort$Block_plot,1,1)

# Rename the herbivory treatments
poc_recruit_mort$Herbivory_trt <- ifelse(poc_recruit_mort$Herbivory_trt=="1X1",Herb.labels[1],
                            ifelse(poc_recruit_mort$Herbivory_trt=="2X2",Herb.labels[2],
                                   ifelse(poc_recruit_mort$Herbivory_trt=="3X3",Herb.labels[3],Herb.labels[4])))
poc_recruit_mort$Herbivory_trt <- factor(poc_recruit_mort$Herbivory_trt, 
                            levels=Herb.labels)

# Check number of species recruited in each year
Rec_year <- poc_recruit_mort %>% filter(Recruitment=="1")
with(Rec_year,table(Timepoint))


# Recruitment by year (this includes corals that died before being sampled at T10)
tmp <- poc_recruit_mort %>% subset(Recruitment=="1")
# Calculate the number of recruits in each Plot.ID in each year.
Recruit_by_year <- aggregate(Recruitment ~ Timepoint + Block_plot + Herbivory_trt + Enrichment_trt, 
                             FUN=sum,
                             data=tmp,
                             drop=F)
# View(Recruit_by_year %>% filter(Block_plot=="A3"))

# Drop the false unwanted treatment combinations created by previous operation
Recruit_by_year <- Recruit_by_year %>% filter((Block_plot=="A3" & Enrichment_trt=="Enriched") |
                                    (Block_plot=="A4" & Enrichment_trt=="Ambient") |
                                    (Block_plot=="B1" & Enrichment_trt=="Enriched") |
                                    (Block_plot=="B2" & Enrichment_trt=="Ambient") |
                                    (Block_plot=="C2" & Enrichment_trt=="Ambient") |
                                    (Block_plot=="C4" & Enrichment_trt=="Enriched") |
                                    (Block_plot=="D1" & Enrichment_trt=="Ambient") |
                                    (Block_plot=="D3" & Enrichment_trt=="Enriched"))
Recruit_by_year$Recruitment[is.na(Recruit_by_year$Recruitment)] <- 0
# View(Recruit_by_year %>% filter(Block_plot=="A3"))

# Add a blocking factor
Recruit_by_year$Block <- substr(Recruit_by_year$Block_plot,1,1)

# tmp <- Recruit_by_year[Recruit_by_year$Timepoint=="T0",]
# aggregate(Recruitment~Herbivory_trt+Enrichment_trt,data=tmp,FUN=sum) # Note there were no recruits in 1x1 Enriched

Recruit_by_year_T0_results <- AnalyzeAICc(dat=Recruit_by_year[Recruit_by_year$Timepoint=="T0",], yname="Recruitment",family="poisson")
Recruit_by_year_T3_results <- AnalyzeAICc(dat=Recruit_by_year[Recruit_by_year$Timepoint=="T3",], yname="Recruitment",family="poisson")
Recruit_by_year_T6_results <- AnalyzeAICc(dat=Recruit_by_year[Recruit_by_year$Timepoint=="T6",], yname="Recruitment",family="poisson")
Recruit_by_year_T9_results <- AnalyzeAICc(dat=Recruit_by_year[Recruit_by_year$Timepoint=="T9",], yname="Recruitment",family="poisson")
Recruit_by_year_results <- list(Recruit_by_year_T0_results=Recruit_by_year_T0_results,
                                Recruit_by_year_T3_results=Recruit_by_year_T3_results,
                                Recruit_by_year_T6_results=Recruit_by_year_T6_results,
                                Recruit_by_year_T9_results=Recruit_by_year_T9_results)

samples <- poc_recruit_mort %>% filter(Rec_timepoint=="T3", Recruitment=="1") %>% select(Un_ID)
tmp <- poc_recruit_mort %>% filter(Un_ID %in% samples$Un_ID, Timepoint=="T6")
Survival_by_year_T3_results <- AnalyzeAICc(dat=tmp, yname="Alive", family="binomial")

samples <- poc_recruit_mort %>% filter(Rec_timepoint=="T6", Recruitment=="1") %>% select(Un_ID)
tmp <- poc_recruit_mort %>% filter(Un_ID %in% samples$Un_ID, Timepoint=="T9")
Survival_by_year_T6_results <- AnalyzeAICc(dat=tmp, yname="Alive", family="binomial")

round(range(plogis(Survival_by_year_T3_results$pred$fit))*100,0) # Results
round(range(plogis(Survival_by_year_T6_results$pred$fit))*100,0) # Results

Survival_by_year_T3_results$Anova_table
Survival_by_year_T6_results$Anova_table

Survival_by_year_results <- list(Survival_by_year_T3_results=Survival_by_year_T3_results,
                                 Survival_by_year_T6_results=Survival_by_year_T6_results)


offsets <- 0.1
quartz(width=4,height=3)
par(mfcol=c(1,2),mar=c(5,3,2,1),oma=c(0,1,0,0))

for(i in 1:2){
# Survival
  tmp.preds <- Survival_by_year_results[i][[1]]$pred
  plot(c(0.7,4.3),c(0,1),type="n",ylab="",xlab="",xaxt="n",yaxt="n",bty="n")
  axis(side=1,at=1:4,labels=NA)
  text(1:4,rep(-0.15,4),Herb.labels,xpd=T,srt=45,adj=1)
  axis(side=2,at=seq(0,1,0.2),las=1)
  
  if(i==1){mtext(side=3,"a) 2019 to 2020",line=0)
    mtext(side=2,"Probability of survival",line=3)}
      
  if(i==2){mtext(side=3,"b) 2020 to 2021",line=0)}
  
  # Add model predictions
  with(tmp.preds[tmp.preds$Enrichment_trt=="Ambient",], 
       points(as.numeric(factor(Herbivory_trt))-offsets,
              plogis(fit),
              pch=19,
              col=adjustcolor(alpha.f = 0.6,"#98c1d9")))
  with(tmp.preds[tmp.preds$Enrichment_trt=="Ambient",], 
       segments(1:4-offsets,plogis(lwr),1:4-offsets,plogis(upr),
                col=adjustcolor(alpha.f = 0.6,"#98c1d9")))
  
  with(tmp.preds[tmp.preds$Enrichment_trt=="Enriched",], 
       points(as.numeric(factor(Herbivory_trt))+offsets,
              plogis(fit),
              pch=19,
              col=adjustcolor(alpha.f = 0.6,"#3d5a80")))
  with(tmp.preds[tmp.preds$Enrichment_trt=="Enriched",], 
       segments(1:4+offsets,plogis(lwr),1:4+offsets,plogis(upr),
                col=adjustcolor(alpha.f = 0.6,"#3d5a80")))
  
  }
legend(1,0.4,c("Ambient","Enriched"),
       col=c(adjustcolor(alpha.f = 0.6,"#98c1d9"),
             adjustcolor(alpha.f = 0.6,"#3d5a80")),
       pch=19,
       cex=0.9,
       bty="n")
mtext(side=1,"Fish consumer pressure",line=-1,outer=T)


