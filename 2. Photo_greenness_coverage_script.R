#photo analysis params, stick together
setwd("C:\\Users\\kelli\\OneDrive - The University of Nottingham\\Nottingham\\Year 2\\UK Lemna - Bradford - Hull\\Flavour analysis\\HS")

#hs <- read.csv("Coverage greenness from photos.csv")
#hs <- read.csv("Coverage greenness from photos + sp.csv") wrong sp labels
hs <- read.csv("Coverage greenness from photos + spchange.csv")
hs <- read.csv("Coverage greenness from photos + spchange + Ljp.csv")

tail(hs)
head(hs)
names(hs)

#STILL TO NORMALISE AREA TO 3 STARTER COLS
# not done as area as a percentage and green area done by pixel avg

#check classes of all cols
sapply(hs, class)
hs$Rep <- as.numeric(hs$Rep)

#factors?
hs$Accession_rep <- as.factor(hs$Accession_rep)
hs$Accession <- as.factor(hs$Accession)
hs$Month <- as.factor(hs$Month)
hs$Batch <- as.factor(hs$Batch)
hs$Pic_rep <- as.factor(hs$Pic_rep)

hs$Species <- as.factor(hs$Species)

class(hs$Accession)
class(hs$Accession_rep)
class(hs$Month)
class(hs$Batch)
class(hs$Pic_rep)
class(hs$Species)

levels(hs$Species) #=4

library(dplyr)
library(stringr)

str(hs)
class(hs)

#subselecting
hs %>% select(Accession, Rep) #just displays them
hs_ord <- hs %>% arrange(Accession, Rep)
unique(hs$Accession, hs$Rep) 
length(unique(hs$Accession,order = ascending)) #25

#boxplots results of paramteres
boxplot(Green_area ~ Accession, hs) #ks12 highest, ks03, 66a low
boxplot(X.Area ~ Accession, hs) #percentage areas
#ks 12 high, ks03, ks13, ks66a low, 
#09, 21, nuff1 low avg but variable

#boxplots results of species
boxplot(Green_area ~ Species, hs) #green area equal between sp
boxplot(X.Area ~ Species, hs) #percentage areas
#L minor lower avg but variable

#variable relationships
#xy plot higher coverage = higher green area
plot(Green_area ~ X.Area, hs)
corr <- cor.test(hs$Green_area, hs$X.Area,
                 method = "pearson"
)
corr$estimate #0.79 strong
corr$p.value #sig

#high coverage = high green area 0.80 v sig

#stats tests

##POST HOC tukey
#diff by all
model <- aov(Green_area~Accession+Rep+Batch, data=hs)
summary(model) # summary of anova
TukeyHSD(model, conf.level=.95)
#ks03 ks66a bad, ks12 good, ks17 good

#diff by all
model1 <- aov(X.Area~Accession+Rep+Batch, data=hs)
summary(model1) # summary of anova
TukeyHSD(model1, conf.level=.95)
#ks03 bad ks66a bad, ks12 good, ks17 good

#diff by all
model2 <- aov(Green_area~Accession+Rep+Month, data=hs)
summary(model2) # summary of anova
TukeyHSD(model2, conf.level=.95)

#diff by all
modelsp1 <- aov(X.Area~Species+Rep+Batch, data=hs)
summary(modelsp1) # summary of anova
TukeyHSD(modelsp1, conf.level=.95)
#ks03 bad ks66a bad, ks12 good, ks17 good

#diff by all
modelsp2 <- aov(Green_area~Species+Rep+Batch, data=hs)
summary(modelsp2) # summary of anova
TukeyHSD(modelsp2, conf.level=.95)

#diff by all
model3 <- aov(X.Area~Accession+Rep+Month, data=hs)
summary(model3) # summary of anova
TukeyHSD(model3, conf.level=.95)

#consider reps seperately as sig?
model1 <- aov(X.Area~Accession_rep, data=hs)
summary(model1) # summary of anova
TukeyHSD(model1, conf.level=.95)
#ks03 bad ks66a bad, ks12 good, ks17 good

param <- hs

library(dplyr)

#summarise data into single values for accession_reps
Summarysp <- param %>% group_by(Species) %>% summarise(Green_area_mean = mean(Green_area), Green_area_stdev = sd(Green_area), Green_area_maximum = max(Green_area))
Summarysp

Summary1sp <- param %>% group_by(Species) %>% summarise(X.Area_mean = mean(X.Area), X.Area_stdev = sd(X.Area), X.Area_n= n(), X.Area_maximum = max(X.Area))
Summary1sp

Summ_sp <- cbind(Summarysp, Summary1sp)

#summarise data into single values for accession_reps
Summary <- param %>% group_by(Accession_rep) %>% summarise(Green_area_mean = mean(Green_area), Green_area_stdev = sd(Green_area), Green_area_n= n(), Green_area_maximum = max(Green_area))
Summary

Summary1 <- param %>% group_by(Accession_rep) %>% summarise(X.Area_mean = mean(X.Area), X.Area_stdev = sd(X.Area), X.Area_n= n(), X.Area_maximum = max(X.Area))
Summary1

Summ_rep <- cbind(Summary, Summary1)

Summ_rep <- Summ_rep[,-6]
Summ_rep %>% arrange(desc(X.Area_mean)) %>% select(Accession_rep, X.Area_mean) %>% top_n(50)

#too much to label plots
barplot(Green_area_mean ~ Accession_rep, Summ_rep)
barplot(X.Area_mean ~ Accession_rep, Summ_rep)

#how do 1 2 3 4 correspond to glasshouse positions, shelves

#see differences between each grouping variable
library(plyr)
#PERCENT COVERAGE
mu <- ddply(hs, "Rep", summarise, grp.mean=mean(X.Area))
head(mu) #1 and 2 higher 49and 53 vs 16 and 15
mu <- ddply(hs, "Accession", summarise, grp.mean=mean(X.Area))
head(mu) #ks12 62 highest, 17 53, 25 51 top 3. bottom = 1.45 66a, ks03 1.9
mu <- ddply(hs, "Accession_rep", summarise, grp.mean=mean(X.Area))
head(mu) #1 2 often better, ks12 did better in all inc 3 and 4
mu <- ddply(hs, "Species", summarise, grp.mean=mean(X.Area))
head(mu) #s poly highest, l minor lowest

library(ggplot2)

#boxplot accession
ETR_Genotype<- ggplot(hs, aes(x=reorder(Accession,X.Area, FUN = median), y=X.Area)) +
  geom_boxplot(aes(fill = Species))+
  scale_color_viridis_d()+
  #facet_wrap(. ~ Bench, nrow = 2, scales = "free_y") +
  #ylim(0,0.750)+
  theme(legend.position="right")+
  theme_classic()+
  theme(panel.spacing=unit(0,"pt"),
        panel.border=element_rect(colour="black", fill=NA),
        panel.background = element_rect(size = 1, linetype = "solid"))+
  theme(axis.text.x = element_text(color="#000000",size=14, angle=45, hjust = 1),
        axis.text.y = element_text(color="#000000",size=14, angle=0),
        title = element_text(size = 14))+
  theme(legend.text=element_text(size=12))+  
  labs(title = expression("(a) Coverage"),
       subtitle = expression(italic("By Genotype")),
       x=expression("Genotype"),
       y=expression("Coverage"))
ETR_Genotype

ETR_Genotype1<- ggplot(hs, aes(x=reorder(Species,X.Area, FUN = median), y=X.Area)) +
  geom_boxplot(aes(fill = Species))+
  scale_color_viridis_d()+
  #facet_wrap(. ~ Bench, nrow = 2, scales = "free_y") +
  #ylim(0,0.750)+
  theme(legend.position="right")+
  theme_classic()+
  theme(panel.spacing=unit(0,"pt"),
        panel.border=element_rect(colour="black", fill=NA),
        panel.background = element_rect(size = 1, linetype = "solid"))+
  theme(axis.text.x = element_text(color="#000000",size=14, angle=45, hjust = 1),
        axis.text.y = element_text(color="#000000",size=14, angle=0),
        title = element_text(size = 14))+
  theme(legend.text=element_text(size=12))+  
  labs(title = expression("(b) Coverage"),
       subtitle = expression(italic("By Genotype")),
       x=expression("Genotype"),
       y=expression("Coverage"))
ETR_Genotype1

#boxplot
ETR_Genotype2<- ggplot(hs, aes(x=reorder(Accession,Green_area, FUN = median), y=Green_area)) +
  geom_boxplot(aes(fill = Species))+
  scale_color_viridis_d()+
  #facet_wrap(. ~ Bench, nrow = 2, scales = "free_y") +
  #ylim(0,0.750)+
  theme(legend.position="right")+
  theme_classic()+
  theme(panel.spacing=unit(0,"pt"),
        panel.border=element_rect(colour="black", fill=NA),
        panel.background = element_rect(size = 1, linetype = "solid"))+
  theme(axis.text.x = element_text(color="#000000",size=14, angle=45, hjust = 1),
        axis.text.y = element_text(color="#000000",size=14, angle=0),
        title = element_text(size = 14))+
  theme(legend.text=element_text(size=12))+  
  labs(title = expression("(b) Greenness"),
       subtitle = expression(italic("By Genotype")),
       x=expression("Genotype"),
       y=expression("Greenness"))
ETR_Genotype2

#stick together
#looks bad, bad method or very variable
library(gridExtra)
grid.arrange(ETR_Genotype, ETR_Genotype2, ncol = 2)

param <- hs

#summarise data into single values for accessions
#not working
Summary <- param %>% group_by(Accession) %>% summarise(Green_area_mean = mean(Green_area), Green_area_stdev = sd(Green_area), Green_area_n= n(), Green_area_maximum = max(Green_area))
Summary

Summary1 <- param %>% group_by(Accession) %>% summarise(X.Area_mean = mean(X.Area), X.Area_stdev = sd(X.Area), X.Area_n= n(), X.Area_maximum = max(X.Area))
Summary1

#summarise data into 4 values for each accession rep
Summary <- hs %>% group_by(Accession_rep) %>% summarise(Green_area_mean = mean(Green_area), Green_area_stdev = sd(Green_area), Green_area_n= n(), Green_area_maximum = max(Green_area))
Summary

Summary1 <- hs %>% group_by(Accession_rep) %>% summarise(X.Area_mean = mean(X.Area), X.Area_stdev = sd(X.Area), X.Area_n= n(), X.Area_maximum = max(X.Area))
Summary1

Summ <- cbind(Summary, Summary1)

#boxplots results of paramteres
#narrowed down to 4 results instead of multiples within reps
barplot(Green_area_mean ~ Accession_rep, Summ)
barplot(X.Area_mean ~ Accession_rep, Summ)

write.csv(Summ, "Green_area_%cov_summary_Ljp.csv")

#order accessions by green area and area mean
hs <- read.csv("Green_area_%cov_summary_addcols.csv")
library(ggplot2)
#change order of species
hs$Species <- factor(hs$Species, levels = c("L. minor", "L. japonica", "L. turionifera", "L. minuta", "S. polyrhiza"))

#redo boxplots
ETR_Genotypeslim<- ggplot(hs, aes(x=reorder(Accession,X.Area_mean, FUN = median), y=X.Area_mean)) +
  geom_boxplot(aes(fill = Species))+
  scale_fill_manual(values = c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"))+
  #scale_color_viridis_d()+
  #facet_wrap(. ~ Bench, nrow = 2, scales = "free_y") +
  #ylim(0,0.750)+
  theme(legend.position="right")+
  theme_classic()+
  theme(panel.spacing=unit(0,"pt"),
        panel.border=element_rect(colour="black", fill=NA),
        panel.background = element_rect(size = 1, linetype = "solid"))+
  theme(axis.text.x = element_text(color="#000000",size=12, angle=90, hjust = 1),
        axis.text.y = element_text(color="#000000",size=14, angle=0),
        title = element_text(size = 14))+
  theme(legend.text=element_text(size=12))+  
  labs(title = expression("(a) Coverage"),
       #subtitle = expression(italic("By Genotype")),
       x=expression("Genotype"),
       y=expression("Coverage"))
ETR_Genotypeslim

ETR_Genotype2slim<- ggplot(hs, aes(x=reorder(Accession,Green_area_mean, FUN = median), y=Green_area_mean)) +
  geom_boxplot(aes(fill = Species))+
  scale_fill_manual(values = c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"))+
  #scale_color_viridis_d()+
  #facet_wrap(. ~ Bench, nrow = 2, scales = "free_y") +
  #ylim(0,0.750)+
  theme(legend.position="right")+
  theme_classic()+
  theme(panel.spacing=unit(0,"pt"),
        panel.border=element_rect(colour="black", fill=NA),
        panel.background = element_rect(size = 1, linetype = "solid"))+
  theme(axis.text.x = element_text(color="#000000",size=12, angle=90, hjust = 1),
        axis.text.y = element_text(color="#000000",size=14, angle=0),
        title = element_text(size = 14))+
  theme(legend.text=element_text(size=12))+  
  labs(title = expression("(b) Greenness"),
       #subtitle = expression(italic("By Genotype")),
       x=expression("Genotype"),
       y=expression("Greenness"))
ETR_Genotype2slim

library(gridExtra)
pdf('Coverage_greenness_boxplotsperaccession+spcoloring_Ljp.pdf', width=22, height=13)

#best way to save
tiff('Coverage_greenness_boxplotsperaccession+spcoloring_Ljp.tiff', units="in", width=22, height=13, res=300, compression = 'lzw')
grid.arrange(ETR_Genotypeslim, ETR_Genotype2slim, ncol=2, nrow=2)
dev.off()
#dev.new()

#summarise data into inidivual value per accssion?
#NOT WORKING?
as.factor(hs$Accession) #25 levels
#still calling as accession_reps
Summary <- hs %>% group_by(Accession) %>% summarise(Green_area_meana = mean(Green_area_mean), Green_area_stdeva = sd(Green_area_mean), Green_area_maximuma = max(Green_area_mean))
Summary

Summary1 <- hs %>% group_by(Accession) %>% summarise(X.Area_mean = mean(X.Area_mean), X.Area_stdev = sd(X.Area_mean), X.Area_maximum = max(X.Area_mean))
Summary1

Summ <- cbind(Summary, Summary1)

Summr <- hs[,-8]
#not working?
Summr %>% arrange(desc(Green_area_mean)) %>% select(Accession, Green_area_mean) %>% top_n(25)
#1       KS12        160.8170
#2      LY01A        146.1442
#3      KS06B        145.4677
#23      KS13        117.9136
#24      KS03        101.7279
#25     KS66A         95.4045
Summr %>% arrange(desc(X.Area_mean)) %>% select(Accession, X.Area_mean) %>% top_n(25)
#1       KS12    62.03130
#2       KS17    53.93815
#3       KS25    51.65062
#23      KS13    19.02372
#24      KS03     1.96400
#25     KS66A     1.45830

#not normalised area to original size of 3fr col area