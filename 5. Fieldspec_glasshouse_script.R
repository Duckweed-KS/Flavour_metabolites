#flavour hyperspec script, stick together
setwd("C:\\Users\\kelli\\OneDrive - The University of Nottingham\\Nottingham\\Year 2\\UK Lemna - Bradford - Hull\\Flavour analysis\\HS")

#one <- read.csv("Batch1_labelled.csv")
#two <- read.csv("Batch2_labelled.csv")

#with species
#one <- read.csv("Batch1_remWR+sp.csv")
#two <- read.csv("Batch2_remWR+sp.csv")

#with species corrected
one <- read.csv("Batch1_remWR+sp_new.csv")
two <- read.csv("Batch2_remWR+sp_new.csv")


hs <- rbind(one, two)

tail(hs)
head(hs)
names(hs)

#check classes of all cols
sapply(hs, class)
hs$Rep <- as.numeric(hs$Rep)

#factors?
hs$Accession_rep <- as.factor(hs$Accession_rep)
hs$Accession <- as.factor(hs$Accession)
hs$Plant_rep <- as.factor(hs$Plant_rep)
hs$Batch <- as.factor(hs$Batch)
hs$Species <- as.factor(hs$Species)

library(dplyr)
library(stringr)

str(hs)
class(hs)

names(hs) <- str_replace_all(names(hs), 'X', 'R')
length(hs)
names(hs)

#remove wr and blank
hs2 <- hs[!(row.names(hs) %in% c("1","98", "168")),]

#cut out reps that are below a certain value e.g bkg or water
#hs remove row 8, 20, 29, 30, 31, 44, 90, 117, 153, 162,
#233, 234, 292, 293, 294,
hs2 <- hs[!(row.names(hs) %in% c("8","20", "29", "30", "31",
                                 "44", "90", "117", "153", "162",
                                 "233", "234", "292", "293", "294")),]

#subselecting
hs2 %>% select(Accession, Plant_rep) #just displays them
hs_ord <- hs2 %>% arrange(Accession, Plant_rep)
unique(hs2$Accession, hs2$Plant_rep) 
length(unique(hs2$Accession,order = ascending))
hs2 %>% filter(Accession != "WR")

hs %>% select(Species, Accession, Plant_rep)
unique(hs$Species) #5 levels

hs <- hs2

#accessions
#models from data for each individual
hsmodsp <- hs %>% group_by(Accession) %>% mutate(PRI = (R531-R570)/(R531+R570))
hs %>% group_by(Accession) %>% mutate(PRI = (R531-R570)/(R531+R570))
hsmod <- hsmodsp[, which(names(hsmodsp) %in% c("Accession", "Accession_rep", "Batch", "PRI"))] #prints 2 output columns

hsmodsp1 <- hs %>% group_by(Accession) %>% mutate(GM = (R750/R550)-1)
hs %>% group_by(Accession) %>% mutate(GM = (R750/R550)-1)
hsmod1 <- hsmodsp1[, which(names(hsmodsp1) %in% c("Accession", "Accession_rep", "Batch", "GM"))]

hsmodsp2 <- hs %>% group_by(Accession) %>% mutate(GI = (R554/R677))
hs %>% group_by(Accession) %>% mutate(GI = (R554/R677))
hsmod2 <- hsmodsp2[, which(names(hsmodsp2) %in% c("Accession", "Accession_rep", "Batch", "GI"))]

hsmodsp3 <- hs %>% group_by(Accession) %>% mutate(NDVI = (R800-R680)/(R800+R680))
hs %>% group_by(Accession) %>% mutate(NDVI = (R800-R680)/(R800+R680))
hsmod3 <- hsmodsp3[, which(names(hsmodsp3) %in% c("Accession", "Accession_rep", "Batch", "NDVI"))]

#make model four
hsmodsp4 <- hs %>% group_by(Accession) %>% mutate(WI = (R900/R970))
hs %>% group_by(Accession) %>% mutate(WI = (R900/R970))
hsmod4 <- hsmodsp4[, which(names(hsmodsp4) %in% c("Accession", "Accession_rep", "Batch", "WI"))]

#mod 5

hsmodsp5 <- hs %>% group_by(Accession) %>% mutate(OVI = (R760/R730))
hs %>% group_by(Accession) %>% mutate(OVI = (R760/R730))
hsmod5 <- hsmodsp5[, which(names(hsmodsp5) %in% c("Accession", "Accession_rep", "Batch", "OVI"))]

#mod 6

hsmodsp6 <- hs %>% group_by(Accession) %>% mutate(SRPI = (R430/R680))
hs %>% group_by(Accession) %>% mutate(SRPI = (R430/R680))
hsmod6 <- hsmodsp6[, which(names(hsmodsp6) %in% c("Accession", "Accession_rep", "Batch", "SRPI"))]

#mod 7

hsmodsp7 <- hs %>% group_by(Accession) %>% mutate(NDWI1 = (R970-R900)/(R970+R900))
hs %>% group_by(Accession) %>% mutate(NDWI1 = (R970-R900)/(R970+R900))
hsmod7 <- hsmodsp7[, which(names(hsmodsp7) %in% c("Accession", "Accession_rep", "Batch", "NDWI1"))]

#mod 8

hsmodsp8 <- hs %>% group_by(Accession) %>% mutate(NDWI = (R860-R1240)/(R860+R1240))
hs %>% group_by(Accession) %>% mutate(NDWI = (R860-R1240)/(R860+R1240))
hsmod8 <- hsmodsp8[, which(names(hsmodsp8) %in% c("Accession", "Accession_rep", "Batch", "Accession_rep", "NDWI"))]

#ADD mod 9
hsmodsp9 <- hs %>% group_by(Accession) %>% mutate(ARI = (1/R550)-(1/R700))
hs %>% group_by(Accession) %>% mutate(ARI = (1/R550)-(1/R700))
hsmod9 <- hsmodsp9[, which(names(hsmodsp9) %in% c("Accession", "Accession_rep", "Batch", "Accession_rep", "ARI"))]

#parameters per tray
models1 <- cbind(hsmod, hsmod1, hsmod2, hsmod3, hsmod4,
                hsmod5, hsmod6, hsmod7, hsmod8, hsmod9)

#species
#models from data for each individual
hsmodsp <- hs %>% group_by(Species) %>% mutate(PRI = (R531-R570)/(R531+R570))
hs %>% group_by(Species) %>% mutate(PRI = (R531-R570)/(R531+R570))
hsmod <- hsmodsp[, which(names(hsmodsp) %in% c("Species", "Accession_rep", "Batch", "PRI"))] #prints 2 output columns

hsmodsp1 <- hs %>% group_by(Species) %>% mutate(GM = (R750/R550)-1)
hs %>% group_by(Species) %>% mutate(GM = (R750/R550)-1)
hsmod1 <- hsmodsp1[, which(names(hsmodsp1) %in% c("Species", "Accession_rep", "Batch", "GM"))]

hsmodsp2 <- hs %>% group_by(Species) %>% mutate(GI = (R554/R677))
hs %>% group_by(Species) %>% mutate(GI = (R554/R677))
hsmod2 <- hsmodsp2[, which(names(hsmodsp2) %in% c("Species", "Accession_rep", "Batch", "GI"))]

hsmodsp3 <- hs %>% group_by(Species) %>% mutate(NDVI = (R800-R680)/(R800+R680))
hs %>% group_by(Species) %>% mutate(NDVI = (R800-R680)/(R800+R680))
hsmod3 <- hsmodsp3[, which(names(hsmodsp3) %in% c("Species", "Accession_rep", "Batch", "NDVI"))]

#make model four
hsmodsp4 <- hs %>% group_by(Species) %>% mutate(WI = (R900/R970))
hs %>% group_by(Species) %>% mutate(WI = (R900/R970))
hsmod4 <- hsmodsp4[, which(names(hsmodsp4) %in% c("Species", "Accession_rep", "Batch", "WI"))]

#mod 5

hsmodsp5 <- hs %>% group_by(Species) %>% mutate(OVI = (R760/R730))
hs %>% group_by(Species) %>% mutate(OVI = (R760/R730))
hsmod5 <- hsmodsp5[, which(names(hsmodsp5) %in% c("Species", "Accession_rep", "Batch", "OVI"))]

#mod 6

hsmodsp6 <- hs %>% group_by(Species) %>% mutate(SRPI = (R430/R680))
hs %>% group_by(Species) %>% mutate(SRPI = (R430/R680))
hsmod6 <- hsmodsp6[, which(names(hsmodsp6) %in% c("Species", "Accession_rep", "Batch", "SRPI"))]

#mod 7

hsmodsp7 <- hs %>% group_by(Species) %>% mutate(NDWI1 = (R970-R900)/(R970+R900))
hs %>% group_by(Species) %>% mutate(NDWI1 = (R970-R900)/(R970+R900))
hsmod7 <- hsmodsp7[, which(names(hsmodsp7) %in% c("Species", "Accession_rep", "Batch", "NDWI1"))]

#mod 8

hsmodsp8 <- hs %>% group_by(Species) %>% mutate(NDWI = (R860-R1240)/(R860+R1240))
hs %>% group_by(Species) %>% mutate(NDWI = (R860-R1240)/(R860+R1240))
hsmod8 <- hsmodsp8[, which(names(hsmodsp8) %in% c("Species", "Accession_rep", "Batch", "Accession_rep", "NDWI"))]

#NEW mod 9
hsmodsp9 <- hs %>% group_by(Species) %>% mutate(ARI = (1/R550)-(1/R700))
hs %>% group_by(Species) %>% mutate(ARI = (1/R550)-(1/R700))
hsmod9 <- hsmodsp9[, which(names(hsmodsp9) %in% c("Species", "Accession_rep", "Batch", "Accession_rep", "ARI"))]


#accession reps
#models from data for each individual
hsmodsp <- hs %>% group_by(Accession) %>% mutate(PRI = (R531-R570)/(R531+R570))
hs %>% group_by(Accession) %>% mutate(PRI = (R531-R570)/(R531+R570))
hsmod <- hsmodsp[, which(names(hsmodsp) %in% c("Species", "Accession_rep", "Batch", "PRI"))] #prints 2 output columns

hsmodsp1 <- hs %>% group_by(Accession_rep) %>% mutate(GM = (R750/R550)-1)
hs %>% group_by(Accession_rep) %>% mutate(GM = (R750/R550)-1)
hsmod1 <- hsmodsp1[, which(names(hsmodsp1) %in% c("Species", "Accession_rep", "Batch", "GM"))]

hsmodsp2 <- hs %>% group_by(Accession_rep) %>% mutate(GI = (R554/R677))
hs %>% group_by(Accession_rep) %>% mutate(GI = (R554/R677))
hsmod2 <- hsmodsp2[, which(names(hsmodsp2) %in% c("Species", "Accession_rep", "Batch", "GI"))]

hsmodsp3 <- hs %>% group_by(Accession_rep) %>% mutate(NDVI = (R800-R680)/(R800+R680))
hs %>% group_by(Accession_rep) %>% mutate(NDVI = (R800-R680)/(R800+R680))
hsmod3 <- hsmodsp3[, which(names(hsmodsp3) %in% c("Species", "Accession_rep", "Batch", "NDVI"))]

#make model four
hsmodsp4 <- hs %>% group_by(Accession_rep) %>% mutate(WI = (R900/R970))
hs %>% group_by(Accession_rep) %>% mutate(WI = (R900/R970))
hsmod4 <- hsmodsp4[, which(names(hsmodsp4) %in% c("Species", "Accession_rep", "Batch", "WI"))]

#mod 5

hsmodsp5 <- hs %>% group_by(Accession_rep) %>% mutate(OVI = (R760/R730))
hs %>% group_by(Accession_rep) %>% mutate(OVI = (R760/R730))
hsmod5 <- hsmodsp5[, which(names(hsmodsp5) %in% c("Species", "Accession_rep", "Batch", "OVI"))]

#mod 6

hsmodsp6 <- hs %>% group_by(Accession_rep) %>% mutate(SRPI = (R430/R680))
hs %>% group_by(Accession_rep) %>% mutate(SRPI = (R430/R680))
hsmod6 <- hsmodsp6[, which(names(hsmodsp6) %in% c("Species", "Accession_rep", "Batch", "SRPI"))]

#mod 7

hsmodsp7 <- hs %>% group_by(Accession_rep) %>% mutate(NDWI1 = (R970-R900)/(R970+R900))
hs %>% group_by(Accession_rep) %>% mutate(NDWI1 = (R970-R900)/(R970+R900))
hsmod7 <- hsmodsp7[, which(names(hsmodsp7) %in% c("Species", "Accession_rep", "Batch", "NDWI1"))]

#mod 8

hsmodsp8 <- hs %>% group_by(Accession_rep) %>% mutate(NDWI = (R860-R1240)/(R860+R1240))
hs %>% group_by(Accession_rep) %>% mutate(NDWI = (R860-R1240)/(R860+R1240))
hsmod8 <- hsmodsp8[, which(names(hsmodsp8) %in% c("Species", "Accession_rep", "Batch", "Accession_rep", "NDWI"))]

#NEW mod 9
hsmodsp9 <- hs %>% group_by(Accession_rep) %>% mutate(ARI = (1/R550)-(1/R700))
hs %>% group_by(Accession_rep) %>% mutate(ARI = (1/R550)-(1/R700))
hsmod9 <- hsmodsp9[, which(names(hsmodsp9) %in% c("Species", "Accession_rep", "Batch", "Accession_rep", "ARI"))]

#notes
PRI <- (R531-R570)/(R531+R570)
GM <- (R750/R550)-1 
GI <- (R554/R677)
WI <- (R900/R970)
OVI2 <- (R760/R730)
SPRI <- (R430/R680)
NDWI1 <- (R970-R900)/(R970+R900)
NDWI <- (R860-R1240)/(R860+R1240)
ARI <- (1/R550)-(1/R700)

models <- cbind(hsmod, hsmod1, hsmod2, hsmod3, hsmod4,
                hsmod5, hsmod6, hsmod7, hsmod8, hsmod9)

#boxplots results of paramteres
boxplot(GM ~ Accession_rep, models)
boxplot(GI ~ Accession_rep, models)
boxplot(NDVI ~ Accession_rep, models)
boxplot(WI ~ Accession_rep, models)
boxplot(PRI ~ Accession_rep, models)
boxplot(OVI ~ Accession_rep, models)
boxplot(SRPI ~ Accession_rep, models)
boxplot(NDWI1 ~ Accession_rep, models)
boxplot(NDWI ~ Accession_rep, models)


#still incluing WR and blank, need to remove
#not working
#models$Accession <- models$Accession[2:26,]

write.csv(models, "Flavour_glasshouse_params.csv")
write.csv(models, "Flavour_glasshouse_params_red+ARI.csv")
#reduced form as removed blanks, water, min cols
write.csv(models, "Flavour_glasshouse_params_sp.csv")
write.csv(models1, "Flavour_glasshouse_params_accrep+ARI.csv")
#when saving wr and blanks not included, read back in
#and re-run graphs now correct

models <- read.csv("Flavour_glasshouse_params.csv")
models <- read.csv("Flavour_glasshouse_params_red.csv")
models <- read.csv("Flavour_glasshouse_params_red+ARI.csv")

models <- read.csv("Flavour_glasshouse_params_sp.csv")
models$Batch <- as.factor(models$Batch)

#for species
##POST HOC tukey
#diff by species and access rep
model <- aov(OVI~Species+Accession_rep+Batch, data=models)
summary(model) # summary of anova
TukeyHSD(model, conf.level=.95)

#diff by all batch *
model1 <- aov(GI~Species+Accession_rep+Batch, data=models)
summary(model1) # summary of anova
TukeyHSD(model1, conf.level=.95)

#diff by all
model2 <- aov(GM~Species+Accession_rep+Batch, data=models)
summary(model2) # summary of anova
TukeyHSD(model2, conf.level=.95)

#diff by all
model3 <- aov(NDVI~Species+Accession_rep+Batch, data=models)
summary(model3) # summary of anova
TukeyHSD(model3, conf.level=.95)

#diff by Species and accession rep
model4 <- aov(PRI~Species+Accession_rep+Batch, data=models)
summary(model4) # summary of anova
TukeyHSD(model4, conf.level=.95)

#diff by all batch *
model5 <- aov(NDWI~Species+Accession_rep+Batch, data=models)
summary(model5) # summary of anova
TukeyHSD(model5, conf.level=.95)

#no diff
model6 <- aov(NDWI1~Species+Accession_rep+Batch, data=models)
summary(model6) # summary of anova
TukeyHSD(model6, conf.level=.95)

#by species and access rep but not batch
model7 <- aov(SRPI~Species+Accession_rep+Batch, data=models)
summary(model7) # summary of anova
TukeyHSD(model7, conf.level=.95)

#no diff
model8 <- aov(WI~Species+Accession_rep+Batch, data=models)
summary(model8) # summary of anova
TukeyHSD(model8, conf.level=.95)

#no diff
model9 <- aov(ARI~Species+Accession_rep+Batch, data=models)
summary(model9) # summary of anova
TukeyHSD(model9, conf.level=.95)

param <- models1

#summarise data into single values for accessions
library(dplyr)
Summary <- param %>% group_by(Accession) %>% summarise(NDVI_mean = mean(NDVI), NDVI_stdev = sd(NDVI), NDVI_n= n(), NDVI_maximum = max(NDVI))
Summary

Summary1 <- param %>% group_by(Accession) %>% summarise(GI_mean = mean(GI), GI_stdev = sd(GI), GI_n= n(), GI_maximum = max(GI))
Summary1

Summary2 <- param %>% group_by(Accession) %>% summarise(GM_mean = mean(GM), GM_stdev = sd(GM), GM_n= n(), GM_maximum = max(GM))
Summary2

Summary3 <- param %>% group_by(Accession) %>% summarise(PRI_mean = mean(PRI), PRI_stdev = sd(PRI), PRI_n= n(), PRI_maximum = max(PRI))
Summary3

Summary4 <- param %>% group_by(Accession) %>% summarise(WI_mean = mean(WI), WI_stdev = sd(WI), WI_n= n(), WI_maximum = max(WI))
Summary4

Summary5 <- param %>% group_by(Accession) %>% summarise(OVI_mean = mean(OVI), OVI_stdev = sd(OVI), OVI_n= n(), OVI_maximum = max(OVI))
Summary5

Summary6 <- param %>% group_by(Accession) %>% summarise(SRPI_mean = mean(SRPI), SRPI_stdev = sd(SRPI), SRPI_n= n(), SRPI_maximum = max(SRPI))
Summary6

Summary7 <- param %>% group_by(Accession) %>% summarise(NDWI1_mean = mean(NDWI1), NDWI1_stdev = sd(NDWI1), NDWI1_n= n(), NDWI1_maximum = max(NDWI1))
Summary7

Summary8 <- param %>% group_by(Accession) %>% summarise(NDWI_mean = mean(NDWI), NDWI_stdev = sd(NDWI), NDWI_n= n(), NDWI_maximum = max(NDWI))
Summary8

Summary9 <- param %>% group_by(Accession...2) %>% summarise(ARI_mean = mean(ARI), ARI_stdev = sd(ARI), ARI_n= n(), ARI_maximum = max(ARI))
Summary9

Summ <- cbind(Summary, Summary1, Summary2, Summary3, Summary4, 
              Summary5, Summary6, Summary7, Summary8)

#summaries per accessiom from params from each accession rep
#boxplots results of paramteres
boxplot(GM ~ Accession, models1)
boxplot(GI ~ Accession, models1)
boxplot(NDVI ~ Accession, models1)
boxplot(WI ~ Accession, models1)
boxplot(PRI ~ Accession, models1)
boxplot(OVI ~ Accession, models1)
boxplot(SRPI ~ Accession, models1)
boxplot(NDWI1 ~ Accession, models1)
boxplot(NDWI ~ Accession, models1)

#summary by species
#summarise data into single values for accessions
param <- models

Summary <- param %>% group_by(Species) %>% summarise(NDVI_mean = mean(NDVI), NDVI_stdev = sd(NDVI), NDVI_n= n(), NDVI_maximum = max(NDVI))
Summary

Summary1 <- param %>% group_by(Species) %>% summarise(GI_mean = mean(GI), GI_stdev = sd(GI), GI_n= n(), GI_maximum = max(GI))
Summary1

Summary2 <- param %>% group_by(Species) %>% summarise(GM_mean = mean(GM), GM_stdev = sd(GM), GM_n= n(), GM_maximum = max(GM))
Summary2

Summary3 <- param %>% group_by(Species) %>% summarise(PRI_mean = mean(PRI), PRI_stdev = sd(PRI), PRI_n= n(), PRI_maximum = max(PRI))
Summary3

Summary4 <- param %>% group_by(Species) %>% summarise(WI_mean = mean(WI), WI_stdev = sd(WI), WI_n= n(), WI_maximum = max(WI))
Summary4

Summary5 <- param %>% group_by(Species) %>% summarise(OVI_mean = mean(OVI), OVI_stdev = sd(OVI), OVI_n= n(), OVI_maximum = max(OVI))
Summary5

Summary6 <- param %>% group_by(Species) %>% summarise(SRPI_mean = mean(SRPI), SRPI_stdev = sd(SRPI), SRPI_n= n(), SRPI_maximum = max(SRPI))
Summary6

Summary7 <- param %>% group_by(Species) %>% summarise(NDWI1_mean = mean(NDWI1), NDWI1_stdev = sd(NDWI1), NDWI1_n= n(), NDWI1_maximum = max(NDWI1))
Summary7

Summary8 <- param %>% group_by(Species) %>% summarise(NDWI_mean = mean(NDWI), NDWI_stdev = sd(NDWI), NDWI_n= n(), NDWI_maximum = max(NDWI))
Summary8

Summsp <- cbind(Summary, Summary1, Summary2, Summary3, Summary4, 
              Summary5, Summary6, Summary7, Summary8)

#boxplots results of paramteres
boxplot(GM ~ Species, models)
boxplot(GI ~ Species, models)
boxplot(NDVI ~ Species, models)
boxplot(WI ~ Species, models)
boxplot(PRI ~ Species, models)
boxplot(OVI ~ Species, models)
boxplot(SRPI ~ Species, models)
boxplot(NDWI1 ~ Species, models)
boxplot(NDWI ~ Species, models)

#accession reps
param <- models

#summarise data into single values for accessions
Summary <- param %>% group_by(Accession_rep) %>% summarise(NDVI_mean = mean(NDVI), NDVI_stdev = sd(NDVI), NDVI_n= n(), NDVI_maximum = max(NDVI))
Summary

Summary1 <- param %>% group_by(Accession_rep) %>% summarise(GI_mean = mean(GI), GI_stdev = sd(GI), GI_n= n(), GI_maximum = max(GI))
Summary1

Summary2 <- param %>% group_by(Accession_rep) %>% summarise(GM_mean = mean(GM), GM_stdev = sd(GM), GM_n= n(), GM_maximum = max(GM))
Summary2

Summary3 <- param %>% group_by(Accession_rep) %>% summarise(PRI_mean = mean(PRI), PRI_stdev = sd(PRI), PRI_n= n(), PRI_maximum = max(PRI))
Summary3

Summary4 <- param %>% group_by(Accession_rep) %>% summarise(WI_mean = mean(WI), WI_stdev = sd(WI), WI_n= n(), WI_maximum = max(WI))
Summary4

Summary5 <- param %>% group_by(Accession_rep) %>% summarise(OVI_mean = mean(OVI), OVI_stdev = sd(OVI), OVI_n= n(), OVI_maximum = max(OVI))
Summary5

Summary6 <- param %>% group_by(Accession_rep) %>% summarise(SRPI_mean = mean(SRPI), SRPI_stdev = sd(SRPI), SRPI_n= n(), SRPI_maximum = max(SRPI))
Summary6

Summary7 <- param %>% group_by(Accession_rep) %>% summarise(NDWI1_mean = mean(NDWI1), NDWI1_stdev = sd(NDWI1), NDWI1_n= n(), NDWI1_maximum = max(NDWI1))
Summary7

Summary8 <- param %>% group_by(Accession_rep) %>% summarise(NDWI_mean = mean(NDWI), NDWI_stdev = sd(NDWI), NDWI_n= n(), NDWI_maximum = max(NDWI))
Summary8

Summacrep <- cbind(Summary, Summary1, Summary2, Summary3, Summary4, 
              Summary5, Summary6, Summary7, Summary8)


#boxplots results of paramteres accessions
barplot(GM_mean ~ Accession, Summ)
barplot(GI_mean ~ Accession, Summ)
barplot(NDVI_mean ~ Accession, Summ)
barplot(WI_mean ~ Accession, Summ)
barplot(PRI_mean ~ Accession, Summ)
barplot(OVI_mean ~ Accession, Summ)
barplot(SRPI_mean ~ Accession, Summ)
barplot(NDWI1_mean ~ Accession, Summ)
barplot(NDWI_mean ~ Accession, Summ)

#species
#boxplots results of paramteres accessions
barplot(GM_mean ~ Species, Summsp)
barplot(GI_mean ~ Species, Summsp)
barplot(NDVI_mean ~ Species, Summsp)
barplot(WI_mean ~ Species, Summsp)
barplot(PRI_mean ~ Species, Summsp)
barplot(OVI_mean ~ Species, Summsp)
barplot(SRPI_mean ~ Species, Summsp)
barplot(NDWI1_mean ~ Species, Summsp)
barplot(NDWI_mean ~ Species, Summsp)

#old form #write.csv(Summ, "HS_param_summ.csv")
#save reduced form
write.csv(Summ, "HS_param_summ.csv")
write.csv(Summsp, "HS_param_summ_sp.csv")
#write.csv(Summacrep, "HS_param_summ_accrep.csv")

#stats tests

##POST HOC tukey
#diff by accession and accession rep
model <- aov(OVI~Accession+Accession_rep+Batch, data=models1)
summary(model) # summary of anova
TukeyHSD(model, conf.level=.95)

#diff by accessions and accession reps batch *
model1 <- aov(GI~Accession+Accession_rep+Batch, data=models1)
summary(model1) # summary of anova
TukeyHSD(model1, conf.level=.95)

#diff by all *** green model differed by batches
model2 <- aov(GM~Accession+Accession_rep+Batch, data=models1)
summary(model2) # summary of anova
TukeyHSD(model2, conf.level=.95)

#diff by all nedvi differed by batches
model3 <- aov(NDVI~Accession+Accession_rep+Batch, data=models1)
summary(model3) # summary of anova
TukeyHSD(model3, conf.level=.95)

#diff by accession and accession rep
model4 <- aov(PRI~Accession+Accession_rep+Batch, data=models1)
summary(model4) # summary of anova
TukeyHSD(model4, conf.level=.95)

#diff by all batch *
model5 <- aov(NDWI~Accession+Accession_rep+Batch, data=models1)
summary(model5) # summary of anova
TukeyHSD(model5, conf.level=.95)

#ns
model6 <- aov(NDWI1~Accession+Accession_rep+Batch, data=models1)
summary(model6) # summary of anova
TukeyHSD(model6, conf.level=.95)

#by accession and access rep but not batch
model7 <- aov(SRPI~Accession+Accession_rep+Batch, data=models1)
summary(model7) # summary of anova
TukeyHSD(model7, conf.level=.95)

#ns
model8 <- aov(WI~Accession+Accession_rep+Batch, data=models1)
summary(model8) # summary of anova
TukeyHSD(model8, conf.level=.95)

#water params dont vary, greenness params vary between batches
#ovi, pri, spri differ between reps and accessions
#gm and ndvi and gi less so between all 3

#species
##POST HOC tukey
#NS BY SPECIES
model <- aov(OVI~Species+Accession_rep+Batch, data=models)
summary(model) # summary of anova
TukeyHSD(model, conf.level=.95)

#*** species
model1 <- aov(GI~Species+Accession_rep+Batch, data=models)
summary(model1) # summary of anova
TukeyHSD(model1, conf.level=.95)

#species **
model2 <- aov(GM~Species+Accession_rep+Batch, data=models)
summary(model2) # summary of anova
TukeyHSD(model2, conf.level=.95)

#sig species
model3 <- aov(NDVI~Species+Accession_rep+Batch, data=models)
summary(model3) # summary of anova
TukeyHSD(model3, conf.level=.95)

#sig species
model4 <- aov(PRI~Species+Accession_rep+Batch, data=models)
summary(model4) # summary of anova
TukeyHSD(model4, conf.level=.95)

#species ***
model5 <- aov(NDWI~Species+Accession_rep+Batch, data=models)
summary(model5) # summary of anova
TukeyHSD(model5, conf.level=.95)

#ns species
model6 <- aov(NDWI1~Species+Accession_rep+Batch, data=models)
summary(model6) # summary of anova
TukeyHSD(model6, conf.level=.95)

#species ***
model7 <- aov(SRPI~Species+Accession_rep+Batch, data=models)
summary(model7) # summary of anova
TukeyHSD(model7, conf.level=.95)

#ns
model8 <- aov(WI~Species+Accession_rep+Batch, data=models)
summary(model8) # summary of anova
TukeyHSD(model8, conf.level=.95)
