#flavour light data from glasshouse
setwd("C:\\Users\\kelli\\OneDrive - The University of Nottingham\\Nottingham\\Year 2\\UK Lemna - Bradford - Hull\\Flavour analysis\\HS")

li <- read.csv("Light_spectro.csv") #wrong species
li <- read.csv("Light_spectro_spchange.csv")

library(dplyr)
#summarise light by accessions
Summary <- li %>% group_by(Accession) %>% summarise(Dec_Light_mean = mean(Dec_Light), Dec_Light_stdev = sd(Dec_Light), Dec_Light_n= n())
Summary 
Summary1 <- li %>% group_by(Accession) %>% summarise(Dec_PPFD_mean = mean(Dec_PPFD), Dec_PPFD_stdev = sd(Dec_PPFD))
Summary1
Summary2 <- li %>% group_by(Accession) %>% summarise(Dec_PFD_mean = mean(Dec_PFD), Dec_PFD_stdev = sd(Dec_PFD))
Summary2 
Summary3 <- li %>% group_by(Accession) %>% summarise(Dec_B_mean = mean(Dec_B), Dec_B_stdev = sd(Dec_B))
Summary3
Summary4 <- li %>% group_by(Accession) %>% summarise(Dec_G_mean = mean(Dec_G), Dec_G_stdev = sd(Dec_G))
Summary4
Summary5 <- li %>% group_by(Accession) %>% summarise(Dec_R_mean = mean(Dec_R), Dec_R_stdev = sd(Dec_R))
Summary5
Summary6 <- li %>% group_by(Accession) %>% summarise(Dec_UV_mean = mean(Dec_UV), Dec_UV_stdev = sd(Dec_UV))
Summary6
Summary7 <- li %>% group_by(Accession) %>% summarise(Dec_FR_mean = mean(Dec_FR), Dec_FR_stdev = sd(Dec_FR))
Summary7
Summary8 <- li %>% group_by(Accession) %>% summarise(Dec_LambdaP_mean = mean(Dec_LambdaP), Dec_LambdaP_stdev = sd(Dec_LambdaP))
Summary8
Summary9 <- li %>% group_by(Accession) %>% summarise(Feb1_Light_mean = mean(Feb1_Light, na.rm = TRUE), Feb1_Light_stdev = sd(Feb1_Light, na.rm = TRUE))
Summary9
Summary10 <- li %>% group_by(Accession) %>% summarise(Feb2_Light_mean = mean(Feb2_Light, na.rm = TRUE), Feb2_Light_stdev = sd(Feb2_Light, na.rm = TRUE))
Summary10
Summary11 <- li %>% group_by(Accession) %>% summarise(Feb2_Digi.fc._mean = mean(Feb2_Digi.fc.), Feb2_Digi_stdev = sd(Feb2_Digi.fc.))
Summary11
Summary12 <- li %>% group_by(Accession) %>% summarise(Feb2_PPFD_mean = mean(Feb2_PPFD), Feb2_PPFD_stdev = sd(Feb2_PPFD))
Summary12
Summary13 <- li %>% group_by(Accession) %>% summarise(Feb2_PFD_mean = mean(Feb2_PFD), Feb2_PFD_stdev = sd(Feb2_PFD))
Summary13 
Summary14 <- li %>% group_by(Accession) %>% summarise(Feb2_B_mean = mean(Feb2_B), Feb2_B_stdev = sd(Feb2_B))
Summary14
Summary15 <- li %>% group_by(Accession) %>% summarise(Feb2_G_mean = mean(Feb2_G), Feb2_G_stdev = sd(Feb2_G))
Summary15
Summary16 <- li %>% group_by(Accession) %>% summarise(Feb2_R_mean = mean(Feb2_R), Feb2_R_stdev = sd(Feb2_R))
Summary16
Summary17 <- li %>% group_by(Accession) %>% summarise(Feb2_UV_mean = mean(Feb2_UV), Feb2_UV_stdev = sd(Feb2_UV))
Summary17
Summary18 <- li %>% group_by(Accession) %>% summarise(Feb2_FR_mean = mean(Feb2_FR), Feb2_FR_stdev = sd(Feb2_FR))
Summary18
Summary19 <- li %>% group_by(Accession) %>% summarise(Feb2_LambdaP_mean = mean(Feb2_LambdaP), Feb2_LambdaP_stdev = sd(Feb2_LambdaP))
Summary19

Summ <- cbind(Summary, Summary1, Summary2, Summary3, Summary4, 
              Summary5, Summary6, Summary7, Summary8, Summary9,
              Summary10, Summary11, Summary12, Summary13, Summary14, 
              Summary15, Summary16, Summary17, Summary18, Summary19)

write.csv(Summ,"Light_spectro_summ.csv")

#now grouped by accession so 4 reps combined, just compare diffs between accessions?

names(Summ)

names(li)
#see differences between dec and feb light levels
#xy plot higher coverage = higher green area
plot(Dec_PPFD ~ Feb2_PPFD, col="black", pch=19, cex=1.5, data=li)
#text(li$Dec_PPFD, li$Feb2_PPFD, labels=li$Accession)
text(Dec_PPFD ~ Feb2_PPFD, labels=Accession_rep,data=li, cex=0.9, font=2)

plot(Dec_PFD ~ Feb2_PFD, li)
plot(Dec_B ~ Feb2_B, li)
plot(Dec_G ~ Feb2_G, li)
plot(Dec_R ~ Feb2_R, li)
plot(Dec_UV ~ Feb2_UV, li)
plot(Dec_FR ~ Feb2_FR, li)
plot(Dec_LambdaP ~ Feb2_LambdaP, li) #doesnt work well
plot(Feb2_Digi ~ Feb2_Light, li)
plot(Feb2_Digi ~ Feb1_Light, li)
plot(Feb2_Digi ~ Dec_Light, li)

corr <- cor.test(li$Dec_PPFD, li$Feb2_PPFD,
                 method = "pearson"
)
corr$estimate #weak
corr$p.value #ns

#are all components stronger in Feb than Dec?
corr <- cor.test(li$Dec_PFD, li$Feb2_PFD,
                 method = "pearson"
)
corr$estimate #weak
corr$p.value #ns

corr <- cor.test(li$Dec_R, li$Feb2_R,
                 method = "pearson"
)
corr$estimate #weak
corr$p.value #ns

corr <- cor.test(li$Dec_B, li$Feb2_B,
                 method = "pearson"
)
corr$estimate #weak
corr$p.value #ns

corr <- cor.test(li$Dec_G, li$Feb2_G,
                 method = "pearson"
)
corr$estimate #weak
corr$p.value #ns

corr <- cor.test(li$Dec_UV, li$Feb2_UV,
                 method = "pearson"
)
corr$estimate #weak
corr$p.value #ns

corr <- cor.test(li$Dec_FR, li$Feb2_FR,
                 method = "pearson"
)
corr$estimate #weak
corr$p.value #ns

corr <- cor.test(li$Dec_LambdaP, li$Feb2_LambdaP,
                 method = "pearson"
)
corr$estimate #weak neg correlation
corr$p.value #sig

#treat shelves differently
model3 <- aov(Dec_Light~Shelf, data=li)
summary(model3) # summary of anova
TukeyHSD(model3, conf.level=.95)
model3 <- aov(Feb2_Light~Shelf, data=li)
summary(model3) # summary of anova
TukeyHSD(model3, conf.level=.95)

#batch not different in dec feb1, but is in feb2
model3 <- aov(Dec_Light~Batch, data=li)
summary(model3) # summary of anova
TukeyHSD(model3, conf.level=.95)
model3 <- aov(Feb1_Light~Batch, data=li)
summary(model3) # summary of anova
TukeyHSD(model3, conf.level=.95)
model3 <- aov(Feb2_Light~Batch, data=li)
summary(model3) # summary of anova
TukeyHSD(model3, conf.level=.95)

#not sig diff between accessions
model3 <- aov(Dec_Light~Accession, data=li)
summary(model3) # summary of anova
TukeyHSD(model3, conf.level=.95)
model3 <- aov(Dec_PPFD~Accession, data=li)
summary(model3) # summary of anova
TukeyHSD(model3, conf.level=.95)
model3 <- aov(Dec_PFD~Accession, data=li)
summary(model3) # summary of anova
TukeyHSD(model3, conf.level=.95)
model3 <- aov(Dec_B~Accession, data=li)
summary(model3) # summary of anova
TukeyHSD(model3, conf.level=.95)
model3 <- aov(Dec_G~Accession, data=li)
summary(model3) # summary of anova
TukeyHSD(model3, conf.level=.95)
model3 <- aov(Dec_R~Accession, data=li)
summary(model3) # summary of anova
TukeyHSD(model3, conf.level=.95)
model3 <- aov(Dec_UV~Accession, data=li)
summary(model3) # summary of anova
TukeyHSD(model3, conf.level=.95)
model3 <- aov(Dec_FR~Accession, data=li)
summary(model3) # summary of anova
TukeyHSD(model3, conf.level=.95)
model3 <- aov(Dec_LambdaP~Accession, data=li)
summary(model3) # summary of anova
TukeyHSD(model3, conf.level=.95)
model3 <- aov(Feb2_Light~Accession, data=li)
summary(model3) # summary of anova
TukeyHSD(model3, conf.level=.95)
model3 <- aov(Feb2_PPFD~Accession, data=li)
summary(model3) # summary of anova
TukeyHSD(model3, conf.level=.95)
model3 <- aov(Feb2_PFD~Accession, data=li)
summary(model3) # summary of anova
TukeyHSD(model3, conf.level=.95)
model3 <- aov(Feb2_B~Accession, data=li)
summary(model3) # summary of anova
TukeyHSD(model3, conf.level=.95)
model3 <- aov(Feb2_G~Accession, data=li)
summary(model3) # summary of anova
TukeyHSD(model3, conf.level=.95)
model3 <- aov(Feb2_R~Accession, data=li)
summary(model3) # summary of anova
TukeyHSD(model3, conf.level=.95)
model3 <- aov(Feb2_UV~Accession, data=li)
summary(model3) # summary of anova
TukeyHSD(model3, conf.level=.95)
model3 <- aov(Feb2_FR~Accession, data=li)
summary(model3) # summary of anova
TukeyHSD(model3, conf.level=.95)
model3 <- aov(Feb2_LambdaP~Accession, data=li)
summary(model3) # summary of anova
TukeyHSD(model3, conf.level=.95)

#not sig diff between accessions dec
kruskal.test(Dec_Light ~ Accession, data = li)
kruskal.test(Dec_PPFD ~ Accession, data = li)
kruskal.test(Dec_PFD ~ Accession, data = li)
kruskal.test(Dec_B ~ Accession, data = li)
kruskal.test(Dec_G ~ Accession, data = li)
kruskal.test(Dec_R ~ Accession, data = li)
kruskal.test(Dec_UV ~ Accession, data = li)
kruskal.test(Dec_FR ~ Accession, data = li)
kruskal.test(Dec_LambdaP ~ Accession, data = li)

#not sig diff between accessions dec
kruskal.test(Feb2_Light ~ Accession, data = li)
kruskal.test(Feb2_PPFD ~ Accession, data = li)
kruskal.test(Feb2_PFD ~ Accession, data = li)
kruskal.test(Feb2_B ~ Accession, data = li)
kruskal.test(Feb2_G ~ Accession, data = li)
kruskal.test(Feb2_R ~ Accession, data = li)
kruskal.test(Feb2_UV ~ Accession, data = li)
kruskal.test(Feb2_FR ~ Accession, data = li)
kruskal.test(Feb2_LambdaP ~ Accession, data = li)

#by Shelf
kruskal.test(Dec_Light ~ Shelf, data = li)
kruskal.test(Dec_PPFD ~ Shelf, data = li)
kruskal.test(Dec_PFD ~ Shelf, data = li)
kruskal.test(Dec_B ~ Shelf, data = li)
kruskal.test(Dec_G ~ Shelf, data = li)
kruskal.test(Dec_R ~ Shelf, data = li)
kruskal.test(Dec_UV ~ Shelf, data = li)
kruskal.test(Dec_FR ~ Shelf, data = li)
kruskal.test(Dec_LambdaP ~ Shelf, data = li)

kruskal.test(Feb2_Light ~ Shelf, data = li)
kruskal.test(Feb2_PPFD ~ Shelf, data = li)
kruskal.test(Feb2_PFD ~ Shelf, data = li)
kruskal.test(Feb2_B ~ Shelf, data = li)
kruskal.test(Feb2_G ~ Shelf, data = li)
kruskal.test(Feb2_R ~ Shelf, data = li)
kruskal.test(Feb2_UV ~ Shelf, data = li)
kruskal.test(Feb2_FR ~ Shelf, data = li)
kruskal.test(Feb2_LambdaP ~ Shelf, data = li)

# all light vars vary between 4 shelves at both time points
#how, which shelves and which variables?
#load library
library(FSA)

#perform Dunn's Test with Bonferroni correction for p-values
dunnTest(Dec_Light ~ Shelf,
         data=li,
         method="bonferroni")

boxplot(Dec_Light ~ Shelf, data=li) #2 and 3 higher
boxplot(Dec_PPFD ~ Shelf, data=li) #1 and 2 higher
boxplot(Dec_PFD ~ Shelf, data=li) #1 and 2 higher
boxplot(Dec_FR ~ Shelf, data=li) # 1 and 2 higher
boxplot(Dec_G ~ Shelf, data=li) #1 and 2 higher
boxplot(Dec_B ~ Shelf, data=li) #2 higher
boxplot(Dec_G ~ Shelf, data=li) #1 and 2 higher
boxplot(Dec_UV ~ Shelf, data=li) #r2 higher
boxplot(Dec_LambdaP ~ Shelf, data=li) #1 and 2 higher

boxplot(Feb2_Light ~ Shelf, data=li) #l1 lower
boxplot(Feb2_PPFD ~ Shelf, data=li) #l1 lower
boxplot(Feb2_PFD ~ Shelf, data=li) #l1 lower
boxplot(Feb2_FR ~ Shelf, data=li) # l1 l2 l3 lower
boxplot(Feb2_G ~ Shelf, data=li) #l1 lower
boxplot(Feb2_B ~ Shelf, data=li) #l1 l2 l3 lower
boxplot(Feb2_G ~ Shelf, data=li) #l1 lower
boxplot(Feb2_UV ~ Shelf, data=li) #l1 l2 l3 lower
boxplot(Feb2_LambdaP ~ Shelf, data=li) #l4 higher

#average table for shelves for each parameter
library(plyr)
#PERCENT COVERAGE
mu <- ddply(li, "Shelf", summarise, grp.mean=mean(Feb2_FR), grp.sd=sd(Feb2_FR))
mu1 <- ddply(li, "Shelf", summarise, grp.mean=mean(Feb2_UV), grp.sd=sd(Feb2_UV))
mu2 <- ddply(li, "Shelf", summarise, grp.mean=mean(Feb2_Light), grp.sd=sd(Feb2_Light))
mu3 <- ddply(li, "Shelf", summarise, grp.mean=mean(Feb2_PPFD), grp.sd=sd(Feb2_PPFD))
mu4 <- ddply(li, "Shelf", summarise, grp.mean=mean(Feb2_R), grp.sd=sd(Feb2_R))
mu5 <- ddply(li, "Shelf", summarise, grp.mean=mean(Feb2_B), grp.sd=sd(Feb2_B))
mu6 <- ddply(li, "Shelf", summarise, grp.mean=mean(Feb2_G), grp.sd=sd(Feb2_G))
mu7 <- ddply(li, "Shelf", summarise, grp.mean=mean(Feb2_LambdaP), grp.sd=sd(Feb2_LambdaP))
mu8 <- ddply(li, "Shelf", summarise, grp.mean=mean(Feb2_PFD), grp.sd=sd(Feb2_PFD))

Feb <- cbind(mu, mu1, mu2, mu3, mu4, mu5, mu6, mu7, mu8)
colnames(Feb) <- c("Shelf", "FR", "FR", "Shelf", "UV", "UV",
                   "Shelf", "LI", "LI", "Shelf", "PPFD", "PFFD",
                   "Shelf", "R", "R", "Shelf", "B", "B",
                   "Shelf", "G", "G", "Shelf", "LA", "LA",
                   "Shelf", "PFD", "PFD")

lu <- ddply(li, "Shelf", summarise, grp.mean=mean(Dec_FR), grp.sd=sd(Dec_FR))
lu1 <- ddply(li, "Shelf", summarise, grp.mean=mean(Dec_UV), grp.sd=sd(Dec_UV))
lu2 <- ddply(li, "Shelf", summarise, grp.mean=mean(Dec_Light), grp.sd=sd(Dec_Light))
lu3 <- ddply(li, "Shelf", summarise, grp.mean=mean(Dec_PPFD), grp.sd=sd(Dec_PPFD))
lu4 <- ddply(li, "Shelf", summarise, grp.mean=mean(Dec_R), grp.sd=sd(Dec_R))
lu5 <- ddply(li, "Shelf", summarise, grp.mean=mean(Dec_B), grp.sd=sd(Dec_B))
lu6 <- ddply(li, "Shelf", summarise, grp.mean=mean(Dec_G), grp.sd=sd(Dec_G))
lu7 <- ddply(li, "Shelf", summarise, grp.mean=mean(Dec_LambdaP), grp.sd=sd(Dec_LambdaP))
lu8 <- ddply(li, "Shelf", summarise, grp.mean=mean(Dec_PFD), grp.sd=sd(Dec_PFD))

Dec <- cbind(lu, lu1, lu2, lu3, lu4, lu5, lu6, lu7, lu8)
colnames(Dec) <- c("Shelf", "FR", "FR", "Shelf", "UV", "UV",
                    "Shelf", "LI", "LI", "Shelf", "PPFD", "PFFD",
                    "Shelf", "R", "R", "Shelf", "B", "B",
                    "Shelf", "G", "G", "Shelf", "LA", "LA",
                    "Shelf", "PFD", "PFD")

#ignoring reps and looking at species instead
boxplot(Dec_Light ~ Accession, data=li) #2 and 3 higher
boxplot(Dec_PPFD ~ Accession, data=li) #1 and 2 higher
boxplot(Dec_PFD ~ Accession, data=li) #1 and 2 higher
boxplot(Dec_FR ~ Accession, data=li) # 1 and 2 higher
boxplot(Dec_G ~ Accession, data=li) #1 and 2 higher
boxplot(Dec_B ~ Accession, data=li) #2 higher
boxplot(Dec_G ~ Accession, data=li) #1 and 2 higher
boxplot(Dec_UV ~ Accession, data=li) #r2 higher
boxplot(Dec_LambdaP ~ Accession, data=li) #1 and 2 higher

boxplot(Feb2_Light ~ Accession, data=li) #l1 lower
boxplot(Feb2_PPFD ~ Accession, data=li) #l1 lower
boxplot(Feb2_PFD ~ Accession, data=li) #l1 lower
boxplot(Feb2_FR ~ Accession, data=li) # l1 l2 l3 lower
boxplot(Feb2_G ~ Accession, data=li) #l1 lower
boxplot(Feb2_B ~ Accession, data=li) #l1 l2 l3 lower
boxplot(Feb2_G ~ Accession, data=li) #l1 lower
boxplot(Feb2_UV ~ Accession, data=li) #l1 l2 l3 lower
boxplot(Feb2_LambdaP ~ Accession, data=li) #l4 higher

library(plyr)
#by accession
mu <- ddply(li, "Accession", summarise, grp.mean=mean(Feb2_FR), grp.sd=sd(Feb2_FR))
mu1 <- ddply(li, "Accession", summarise, grp.mean=mean(Feb2_UV), grp.sd=sd(Feb2_UV))
mu2 <- ddply(li, "Accession", summarise, grp.mean=mean(Feb2_Light, na.rm=TRUE), grp.sd=sd(Feb2_Light, na.rm=TRUE))
mu3 <- ddply(li, "Accession", summarise, grp.mean=mean(Feb2_PPFD), grp.sd=sd(Feb2_PPFD))
mu4 <- ddply(li, "Accession", summarise, grp.mean=mean(Feb2_R), grp.sd=sd(Feb2_R))
mu5 <- ddply(li, "Accession", summarise, grp.mean=mean(Feb2_B), grp.sd=sd(Feb2_B))
mu6 <- ddply(li, "Accession", summarise, grp.mean=mean(Feb2_G), grp.sd=sd(Feb2_G))
mu7 <- ddply(li, "Accession", summarise, grp.mean=mean(Feb2_LambdaP), grp.sd=sd(Feb2_LambdaP))
mu8 <- ddply(li, "Accession", summarise, grp.mean=mean(Feb2_PFD), grp.sd=sd(Feb2_PFD))
mu9 <- ddply(li, "Accession", summarise, grp.mean=mean(Feb1_Light, na.rm=TRUE), grp.sd=sd(Feb1_Light, na.rm=TRUE))
mu10 <- ddply(li, "Accession", summarise, grp.mean=mean(Feb2_Digi.fc., na.rm=TRUE), grp.sd=sd(Feb2_Digi.fc., na.rm=TRUE))
Feb <- cbind(mu, mu1, mu2, mu3, mu4, mu5, mu6, mu7, mu8, mu9, mu10)
colnames(Feb) <- c("Accession", "FR", "FR", "Accession", "UV", "UV",
                   "Accession", "LI", "LI", "Accession", "PPFD", "PFFD",
                   "Accession", "R", "R", "Accession", "B", "B",
                   "Accession", "G", "G", "Accession", "LA", "LA",
                   "Accession", "PFD", "PFD", "Accession", "LI1", "LI1",
                   "Accession", "Digi", "Digi")

lu <- ddply(li, "Accession", summarise, grp.mean=mean(Dec_FR), grp.sd=sd(Dec_FR))
lu1 <- ddply(li, "Accession", summarise, grp.mean=mean(Dec_UV), grp.sd=sd(Dec_UV))
lu2 <- ddply(li, "Accession", summarise, grp.mean=mean(Dec_Light), grp.sd=sd(Dec_Light))
lu3 <- ddply(li, "Accession", summarise, grp.mean=mean(Dec_PPFD), grp.sd=sd(Dec_PPFD))
lu4 <- ddply(li, "Accession", summarise, grp.mean=mean(Dec_R), grp.sd=sd(Dec_R))
lu5 <- ddply(li, "Accession", summarise, grp.mean=mean(Dec_B), grp.sd=sd(Dec_B))
lu6 <- ddply(li, "Accession", summarise, grp.mean=mean(Dec_G), grp.sd=sd(Dec_G))
lu7 <- ddply(li, "Accession", summarise, grp.mean=mean(Dec_LambdaP), grp.sd=sd(Dec_LambdaP))
lu8 <- ddply(li, "Accession", summarise, grp.mean=mean(Dec_PFD), grp.sd=sd(Dec_PFD))

Dec <- cbind(lu, lu1, lu2, lu3, lu4, lu5, lu6, lu7, lu8)
colnames(Dec) <- c("Accession", "FR", "FR", "Accession", "UV", "UV",
                   "Accession", "LI", "LI", "Accession", "PPFD", "PFFD",
                   "Accession", "R", "R", "Accession", "B", "B",
                   "Accession", "G", "G", "Accession", "LA", "LA",
                   "Accession", "PFD", "PFD")

write.csv(Feb, "Light_Feb_AccSumm.csv")
write.csv(Dec, "Light_Dec_AccSumm.csv")
