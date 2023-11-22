#script to explore duckweed metabolite data

setwd("C:\\Users\\kelli\\OneDrive - The University of Nottingham\\Nottingham\\Year 2\\Metabolomics")

#gr <- read.csv("Metabolite_all.csv")
#to run pca and heat map with no l. gibba
#gr <- read.csv("Metabolite_all_spnew.csv")
gr <- read.csv("Metabolite_all_spnew+Ljp.csv")

tail(gr)
head(gr)
names(gr)

sapply(gr, class) #for factors
#gr1 <- lapply(gr,as.numeric)
#sapply(gr1, class) #for numric cols
#gr1$Accession_rep <- gr$Accession_rep
#gr1$Accession <- gr$Accession
#gr1$Rep<- gr$Rep
#gr1$Sample.no <- gr$Sample.no
#gr1$Species <- gr$Species

library(dplyr)
library(stringr)
library(ggplot2)

#gr1 <- gr

gr %>% select(Accession) #just displays them #96
unique(gr$Accession) #24
length(unique(gr$Accession,order = ascending)) #24
length(unique(gr$Species,order = ascending)) #5

as.factor(gr$Accession_rep)
as.factor(gr$Accession)
as.factor(gr$Rep)
as.factor(gr$Species)

#messes up order in  col so wrong labels
#levels(gr$Species) <- c("L. minor", "L. minuta", "L. gibba",
                        #"L. turionifera", "S. polyrhiza")  

names(gr)

#summarise sds
sum <- gr %>% summarise_all(mean)

#visual of data for each column for accession
#ks03 and ks66a lower in alanine, arginine, asparagine,
#aspartic acid, isoleucine, leucine, phenylalanine,
#proline, serine, threonine, glycine, methionine, tyrosine

par(mfrow = c(1,1)) 
par(mfrow = c(3,5)) 
boxplot(L.Alanine~Accession,data=gr)
boxplot(L.Arginine~Accession,data=gr)
boxplot(L.Asparagine~Accession,data=gr)
boxplot(L.Aspartic.Acid~Accession,data=gr)
boxplot(L.Glutamic.acid~Accession,data=gr)
boxplot(L.Glutamine~Accession,data=gr)
boxplot(L.Isoleucine_1~Accession,data=gr)
boxplot(L.Leucine_2~Accession,data=gr)
boxplot(L.Phenylalanine~Accession,data=gr)
boxplot(L.Proline~Accession,data=gr)
boxplot(L.Serine~Accession,data=gr)
boxplot(L.Threonine~Accession,data=gr)
boxplot(L.Tryptophane~Accession,data=gr)
boxplot(L.Valine~Accession,data=gr)
boxplot(Glycine~Accession,data=gr)
boxplot(L.Histidine~Accession,data=gr)
boxplot(L.Methionine~Accession,data=gr)
boxplot(L.Tyrosine~Accession,data=gr)
boxplot(Tyramine~Accession,data=gr)
boxplot(Tryptamine~Accession,data=gr)
boxplot(Shikimic.acid~Accession,data=gr)

#interesting ks03 ks66a ks12 diff
par(mfrow = c(3,3)) 
boxplot(L.Alanine~Accession,data=gr)
boxplot(L.Isoleucine_1~Accession,data=gr)
boxplot(L.Leucine_2~Accession,data=gr)
boxplot(L.Phenylalanine~Accession,data=gr)
boxplot(L.Proline~Accession,data=gr)
boxplot(L.Serine~Accession,data=gr)
boxplot(L.Threonine~Accession,data=gr)
boxplot(L.Methionine~Accession,data=gr)
boxplot(L.Tyrosine~Accession,data=gr)

#make by ggplot
one <- ggplot(gr, aes(x=reorder(Accession,L.Alanine), y=L.Alanine)) + 
  labs(y= "Alanine") +
  labs(x="") +
  geom_boxplot() +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8))
two <- ggplot(gr, aes(x=reorder(Accession,L.Isoleucine_1), y=L.Isoleucine_1)) + 
  labs(y= "Isoleucine") +
  labs(x="") +
  geom_boxplot() +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8))
three <- ggplot(gr, aes(x=reorder(Accession,L.Leucine_2), y=L.Leucine_2)) + 
  labs(y= "Leucine") +
  labs(x="") +
  geom_boxplot() +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8))
four <- ggplot(gr, aes(x=reorder(Accession,L.Phenylalanine), y=L.Phenylalanine)) + 
  labs(y= "Phenylalanine") +
  labs(x="") +
  geom_boxplot() +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8))
five <- ggplot(gr, aes(x=reorder(Accession,L.Proline), y=L.Proline)) + 
  labs(y= "Proline") +
  labs(x="") +
  geom_boxplot() +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8))
six <- ggplot(gr, aes(x=reorder(Accession,L.Serine), y=L.Serine)) + 
  labs(y= "Serine") +
  labs(x="") +
  geom_boxplot() +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8))
sev <- ggplot(gr, aes(x=reorder(Accession,L.Threonine), y=L.Threonine)) + 
  labs(y= "Threonine") +
  labs(x="") +
  geom_boxplot() +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8))
eight <- ggplot(gr, aes(x=reorder(Accession,L.Methionine), y=L.Methionine)) + 
  labs(y= "Methionine", x="Accession") +
  geom_boxplot() +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8))
nine <- ggplot(gr, aes(x=reorder(Accession,L.Tyrosine), y=L.Tyrosine)) + 
  labs(y= "Tyrosine") +
  labs(x="") +
  geom_boxplot() +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

library(gridExtra)
grid.arrange(one,two,three,four,five,six,sev,eight,nine, ncol=3)
#dev.off()

boxplot(Chlorogenic.acid~Accession,data=gr)
boxplot(putative.Chlorogenic.acid.Isomer~Accession,data=gr)
boxplot(Cyanidine.3.Glc~Accession,data=gr)
boxplot(Luteoline.8.C.Glc~Accession,data=gr)
boxplot(Apigenine.8.C.Glc~Accession,data=gr)
boxplot(Luteoline.7.O.Glc~Accession,data=gr)
boxplot(Apigenine.7.O.Glc~Accession,data=gr)

boxplot(Apigenine~Accession,data=gr)
boxplot(Luteoline~Accession,data=gr)
boxplot(Glucose~Accession,data=gr)
boxplot(Fructose~Accession,data=gr)
boxplot(Sucrose~Accession,data=gr)
boxplot(Starch..Glu.~Accession,data=gr)
boxplot(Starch.mg...g~Accession,data=gr)
boxplot(Cyanidine.3.Glc.1~Accession,data=gr)
boxplot(Cyanidin.Mal.Glc~Accession,data=gr)
boxplot(Chlorogenic.acid.1~Accession,data=gr)
boxplot(Luteolin.7.O.Glc~Accession,data=gr)
boxplot(Luteolin.8.C.Glc~Accession,data=gr)
boxplot(Apigenin.7.O.Glc~Accession,data=gr)
boxplot(Apigenin.8.C.Glc~Accession,data=gr)
boxplot(Luteolin~Accession,data=gr)
boxplot(Apigenin~Accession,data=gr)
boxplot(putative.Chlorogenic.acid.Isomere~Accession,data=gr)

boxplot(Glucose~Accession_rep,data=gr)
boxplot(Fructose~Accession_rep,data=gr)
boxplot(Sucrose~Accession_rep,data=gr)
boxplot(Starch..Glu.~Accession_rep,data=gr)
boxplot(Starch.mg...g~Accession_rep,data=gr)
aov <- aov(Glucose~Accession,data=gr)
summary(aov)
tuk_out <- TukeyHSD(aov, "Accession", conf.level=.95)
str(tuk_out)
stripchart(gr, Glucose ~Accession, pch="|", ylim=c(.5,2.5))

#visual of data for each column by species
par(mfrow = c(1,1)) 
boxplot(L.Alanine~Species,data=gr)
boxplot(L.Arginine~Species,data=gr)
boxplot(L.Asparagine~Species,data=gr)
boxplot(L.Aspartic.Acid~Species,data=gr)
boxplot(L.Glutamic.acid~Species,data=gr)
boxplot(L.Glutamine~Species,data=gr)
boxplot(L.Isoleucine_1~Species,data=gr)
boxplot(L.Leucine_2~Species,data=gr)
boxplot(L.Phenylalanine~Species,data=gr)
boxplot(L.Proline~Species,data=gr)
boxplot(L.Serine~Species,data=gr)
boxplot(L.Threonine~Species,data=gr)
boxplot(L.Tryptophane~Species,data=gr)
boxplot(L.Valine~Species,data=gr)
boxplot(Chlorogenic.acid~Species,data=gr)
boxplot(putative.Chlorogenic.acid.Isomer~Species,data=gr)
boxplot(Cyanidine.3.Glc~Species,data=gr)
boxplot(Luteoline.8.C.Glc~Species,data=gr)
boxplot(Apigenine.8.C.Glc~Species,data=gr)
boxplot(Luteoline.7.O.Glc~Species,data=gr)
boxplot(Apigenine.7.O.Glc~Species,data=gr)
boxplot(Glycine~Species,data=gr)
boxplot(L.Histidine~Species,data=gr)
boxplot(L.Methionine~Species,data=gr)
boxplot(L.Tyrosine~Species,data=gr)
boxplot(Tyramine~Species,data=gr)
boxplot(Tryptamine~Species,data=gr)
boxplot(Shikimic.acid~Species,data=gr)

boxplot(Apigenine~Species,data=gr)
boxplot(Luteoline~Species,data=gr)
boxplot(Glucose~Species,data=gr)
boxplot(Fructose~Species,data=gr)
boxplot(Sucrose~Species,data=gr)
boxplot(Starch..Glu.~Species,data=gr)
boxplot(Starch.mg...g~Species,data=gr)
boxplot(Cyanidine.3.Glc.1~Species,data=gr)
boxplot(Cyanidin.Mal.Glc~Species,data=gr)
boxplot(Chlorogenic.acid.1~Species,data=gr)
boxplot(Luteolin.7.O.Glc~Species,data=gr)
boxplot(Luteolin.8.C.Glc~Species,data=gr)
boxplot(Apigenin.7.O.Glc~Species,data=gr)
boxplot(Apigenin.8.C.Glc~Species,data=gr)
boxplot(Luteolin~Species,data=gr)
boxplot(Apigenin~Species,data=gr)
boxplot(putative.Chlorogenic.acid.Isomere~Species,data=gr)

#label outliers sugar
library(ggplot2)
ggplot(gr, aes(x = factor(Species), y = Glucose, fill = factor(Species))) + 
  geom_boxplot() +
  stat_summary(
    aes(label = round(stat(y), 1)),
    geom = "text", 
    fun.y = function(y) { o <- boxplot.stats(y)$out; if(length(o) == 0) NA else o },
    hjust = -1
  )

names(gr)
#sig effect of accession and species on metabolite?
aov <- aov(L.Alanine~Species*Accession,data=gr)
summary(aov) #n.s, sig by accession
aov <- aov(L.Arginine~Accession,data=gr)
summary(aov) #n.s, sig by accession
aov <- aov(L.Asparagine~Accession,data=gr)
summary(aov)
aov <- aov(Aspartic.Acid~Accession,data=gr)
summary(aov)
aov <- aov(Glutamic.Acid~Accession,data=gr)
summary(aov)
aov <- aov(L.Alanine~Accession,data=gr)
summary(aov)
aov <- aov(L.Alanine~Accession,data=gr)
summary(aov)
aov <- aov(L.Alanine~Accession,data=gr)
summary(aov)
tuk_out <- TukeyHSD(aov, "Accession", conf.level=.95)
str(tuk_out)
tuk_out

aov <- aov(L.Arginine~Species+Accession,data=gr)
summary(aov) #species ***, acc ns
aov <- aov(L.Asparagine~Species+Accession,data=gr)
summary(aov) #species ns, acc ***
aov <- aov(L.Aspartic.Acid~Species+Accession,data=gr)
summary(aov) #sp ** acc *
aov <- aov(L.Glutamic.acid~Species+Accession,data=gr)
summary(aov) #sp ** acc *
aov <- aov(L.Glutamine~Species+Accession,data=gr)
summary(aov) #n.s
aov <- aov(L.Isoleucine_1~Species+Accession,data=gr)
summary(aov) #n.s sp, acc *
aov <- aov(L.Leucine_2~Species+Accession,data=gr)
summary(aov) #ns sp, * accession
aov <- aov(L.Phenylalanine~Species+Accession,data=gr)
summary(aov) #ns sp, * accession
aov <- aov(L.Proline~Species+Accession,data=gr)
summary(aov) #ns sp, * accession
aov <- aov(L.Serine~Species+Accession,data=gr)
summary(aov) #ns sp, ** access
aov <- aov(L.Threonine~Species+Accession,data=gr)
summary(aov) #ns sp, * accession
aov <- aov(L.Tryptophane~Species+Accession,data=gr)
summary(aov) #*** sp, ns acc
aov <- aov(L.Valine~Species+Accession,data=gr)
summary(aov) #sp *, ns acc
aov <- aov(Chlorogenic.acid~Species+Accession,data=gr)
summary(aov) #sp ***, * acc
aov <- aov(putative.Chlorogenic.acid.Isomer~Species+Accession,data=gr)
summary(aov) #sp ns, ns acc
aov <- aov(Cyanidine.3.Glc~Species+Accession,data=gr)
summary(aov) #ns sp, ns accession
aov <- aov(Luteoline.8.C.Glc~Species+Accession,data=gr)
summary(aov) #*** sp, ns accession
aov <- aov(Apigenine.8.C.Glc~Species+Accession,data=gr)
summary(aov) #*** species, ns access
aov <- aov(Luteoline.7.O.Glc~Species+Accession,data=gr)
summary(aov) #*** #sp, access ***
aov <- aov(Apigenine.7.O.Glc~Species+Accession,data=gr)
summary(aov) #*** #species, access***
aov <- aov(Glycine~Species+Accession,data=gr)
summary(aov) #ns species, * access
aov <- aov(L.Histidine~Species+Accession,data=gr)
summary(aov) #**species, ** access
aov <- aov(L.Methionine~Species+Accession,data=gr)
summary(aov) #***species, * access
aov <- aov(L.Tyrosine~Species+Accession,data=gr)
summary(aov) #** species, ** acc
aov <- aov(Tyramine~Species+Accession,data=gr)
summary(aov) #*** species, ns acc
aov <- aov(Shikimic.acid~Species+Accession,data=gr)
summary(aov) #*** species, ns acc
aov <- aov(Apigenine~Species+Accession,data=gr)
summary(aov) #*** species, ns acc
aov <- aov(Luteoline~Species+Accession,data=gr)
summary(aov) #*** species, ns acc
aov <- aov(Glucose~Species+Accession,data=gr)
summary(aov) #ns
aov <- aov(Fructose~Species+Accession,data=gr)
summary(aov) #ns
aov <- aov(Sucrose~Species+Accession,data=gr)
summary(aov) #ns
aov <- aov(Starch..Glu.~Species+Accession,data=gr)
summary(aov) #ns
aov <- aov(Starch.mg...g~Species+Accession,data=gr)
summary(aov) #ns
aov <- aov(Cyanidine.3.Glc.1~Species+Accession,data=gr)
summary(aov) #*** species, *** acc
aov <- aov(Cyanidin.Mal.Glc~Species+Accession,data=gr)
summary(aov) #*** species, *** acc
aov <- aov(Chlorogenic.acid.1~Species+Accession,data=gr)
summary(aov) #ns
aov <- aov(Luteolin.7.O.Glc~Species+Accession,data=gr)
summary(aov) #*** species, *** acc
aov <- aov(Luteolin.8.C.Glc~Species+Accession,data=gr)
summary(aov) #*** species, * acc
aov <- aov(Luteolin.8.C.Glc~Species+Accession,data=gr)
summary(aov) #*** species, * acc
aov <- aov(Apigenin.7.O.Glc~Species+Accession,data=gr)
summary(aov) #*** species, *** acc
aov <- aov(Apigenin.8.C.Glc~Species+Accession,data=gr)
summary(aov) #*** species, *** acc
aov <- aov(Luteolin~Species+Accession,data=gr)
summary(aov) #*** species, ns acc
aov <- aov(Apigenin~Species+Accession,data=gr)
summary(aov) #*** species, ns acc

aov <- aov(Glucose~Accession,data=gr)
summary(aov) #ns
aov <- aov(Fructose~Accession,data=gr)
summary(aov) #ns
aov <- aov(Sucrose~Accession,data=gr)
summary(aov) #* species
aov <- aov(Starch..Glu.~Accession,data=gr)
summary(aov) #ns
aov <- aov(Starch.mg...g~Accession,data=gr)
summary(aov) #ns

#tryptamine 0 not plotted
#hplc chlorogenic isomer 0 not plotted

#remove tryptamine as all 0s 
gr <- gr[,-33]
#avg values and SD per accession
#avgs
sum <- gr %>% group_by(Accession) %>% select(-Accession_rep, -Rep, -Sample.no, -Sample.Position, -Species) %>% summarise_all((mean))
sum
#sds
sum1 <- gr %>% group_by(Accession) %>% select(-Accession_rep, -Rep, -Sample.no, -Sample.Position, -Species) %>% summarise_all(sd)
sum1

#not working, try diff way
library(plyr)
sum_mu <- ddply(gr, "Accession", summarise, grp.mean=mean(L.Alanine), grp.sd=sd(L.Alanine))

#with na's in
sum_sp<- ddply(gr, "Species", summarise, grp.mean=mean(Sucrose, na.rm=TRUE), grp.sd=sd(Sucrose, na.rm=TRUE))
sum_sp1<- ddply(gr, "Species", summarise, grp.mean=mean(Chlorogenic.acid, na.rm=TRUE), grp.sd=sd(Chlorogenic.acid, na.rm=TRUE))
sum_sp2<- ddply(gr, "Species", summarise, grp.mean=mean(putative.Chlorogenic.acid.Isomer, na.rm=TRUE), grp.sd=sd(putative.Chlorogenic.acid.Isomer, na.rm=TRUE))
sum_sp3<- ddply(gr, "Species", summarise, grp.mean=mean(Cyanidine.3.Glc, na.rm=TRUE), grp.sd=sd(Cyanidine.3.Glc, na.rm=TRUE))

#test, see if avgs match old method - does
sum_sp4<- ddply(gr, "Species", summarise, grp.mean=mean(L.Alanine, na.rm=TRUE), grp.sd=sd(L.Alanine, na.rm=TRUE))


acc <- cbind(sum, sum1)

#avgs sds per species
#not working, dply problem?
sum2 <- gr %>% group_by(Species) %>% select(-Accession_rep, -Rep, -Sample.no, -Sample.Position, -Accession) %>% summarise_all(mean)
sum2
#sds
sum3 <- gr %>% group_by(Species) %>% select(-Accession_rep, -Rep, -Sample.no, -Sample.Position, -Accession) %>% summarise_all(sd)
sum3

sp <- cbind(sum2, sum3)

sum2 <- read.csv("Met_summ_sp_newsp.csv")

#write.csv(acc, "Met_summ_acc.csv")
#write.csv(sp, "Met_summ_sp.csv")

#new sp groups
write.csv(acc, "Met_summ_acc_newsp.csv")
write.csv(sum2, "Met_summ_sp_newsp_Ljp.csv")

#remove na's in summaries as did not detect compound
library(tidyr)
gr <- gr %>% mutate_all(funs(replace_na(.,0)))
sum <- sum %>% mutate_all(funs(replace_na(.,0)))
sum2 <- sum2 %>% mutate_all(funs(replace_na(.,0)))

#heat map attempt currently with all raw data for accession reps
library(RColorBrewer)
#scale numeric data
#to 20 in a.as, 27:32
#21:27, 34:36 is sec mets
#36 to 38, 40 sugars
#hplc sec mets 41:49
aas1 <- gr[7:20]
aas2 <- gr[28:32]
sm1 <- gr[21:27]
sm2 <- gr[34:35]
su <- gr[36:38]
su1 <- gr[40]
sm_hplc <- gr[42:49]
access <- gr[1]

aas <- cbind(access, aas1, aas2)
sug <- cbind(access, su, su1)

boxplot(Sucrose~Accession,data=sug)
boxplot(Fructose~Accession,data=sug)
boxplot(Glucose~Accession,data=sug)
boxplot(Starch.mg...g~Accession,data=sug)

aov <- aov(Glucose~Accession,data=sug)
summary(aov) #ns
aov <- aov(Fructose~Accession,data=sug)
summary(aov) #ns
aov <- aov(Sucrose~Accession,data=sug)
summary(aov) #ns
aov <- aov(Starch..Glu.~Accession,data=sug)
summary(aov) #ns
aov <- aov(Starch.mg...g~Accession,data=sug)
summary(aov) #ns

boxplot(Sucrose~Species,data=gr)
boxplot(Fructose~Species,data=gr)
boxplot(Glucose~Species,data=gr)
boxplot(Starch.mg...g~Species,data=gr)

#summarise aas per accession
aa_sum <- aas %>% group_by(Accession) %>% summarise_all(mean)
aa_sum
#sds
#aa_sum1 <- aas %>% group_by(Species) summarise_all(sd)
#aa_sum1

names(aas)

puttog <- cbind(aas1, aas2, su, sm1, sm2, sm_hplc)

ions_all <- aas[2:19]
ions_all <- aa_sum[2:19]
ions_all <- sug[2:5]
ions_all <- sum_sug[2:5]

#remove tryptamine as all 0s 
#gr <- gr[,-33]
ions_all <- gr[7:49]
ions_all <- ions_all[-33] #remove double starch
ions_sc <- scale(ions_all,center=T,scale=T)
ions_sc <- as.matrix(ions_sc)

class(ions_sc)
as.numeric(ions_sc)

#assign row names from df
rownames(ions_sc) <- gr$Accession
#or
names(aa_sum)
names(aa_sum) <- c("Accession", "Alanine", "Arginine", "Asparagine",
                  "Aspartic acid", "Glutamic acid", "Glutamine",
                  "Isoleucine", "Leucine", "Phenylalanine", "Proline",
                  "Serine", "Threonine", "Tryptophane", "Valine", "Glycine",
                  "Histidine", "Methionine", "Tyrosine", "Thyramine")
#remove L. and _1 and 2_ in names
rownames(ions_sc) <- aa_sum$Accession
rownames(ions_sc) <- sug$Accession
rownames(ions_sc) <- sum_sug$Accession

heatmap(ions_sc)
heatmap(ions_sc, Colv = NA, Rowv = NA, scale="column")
#change color scheme yellow blue
colfunc <- colorRampPalette(c("blue", "yellow"))
heatmap(ions_sc,col=colfunc(11),scale="row")

#use avgs per accessions to make heatmap
#24 levels
#subset so numeric
sum <- sum[,-18] #remove cyanidine as all nas
sume <- sum[,2:43]
ions_sc <- scale(sume,center=T,scale=T)
ions_sc <- as.matrix(ions_sc)
#assign row names from df
rownames(ions_sc) <- sum$Accession
colfunc <- colorRampPalette(c("blue", "yellow"))
heatmap(ions_sc,col=colfunc(11),scale="column")
heatmap(ions_sc,col=colfunc(11),scale="row")
#no obvious groupings of accession

#BEST HEATPLOT BY SPECIES
#group by species
#sum2 <- sum2[,-18] #remove cyanidine as all nas
library(tidyr)
#try read in new sum2 adjusted with ones that had nas, sumsp1, sumsp2 etc
sum2 <- read.csv("Met_summ_sp_newsp_Ljp_renames.csv")
#sum2 <- sum2 %>% mutate_all(funs(replace_na(.,0)))
sumx <- sum2[,2:44]
ions_sc <- scale(sumx,center=T,scale=T)
ions_sc <- as.matrix(ions_sc)
#assign row names from df
rownames(ions_sc) <- sum2$Species
ions_sc <- ions_sc[,-27]

#png cutting half data off
par(mar=c(7,4,4,2)+0.1)
#png(filename='Heatmap_avgs_species.png', width=800, height=500)
png(filename='Heatmap_avgs_species_newsp.png', width=800, height=500)
colfunc <- colorRampPalette(c("blue", "yellow"))
#heatmap(ions_sc,col=colfunc(15),scale="row",cexCol=0.9,margins=c(12,8))
heatmap(ions_sc,col=colfunc(15),scale="column",cexCol=0.9,margins=c(12,8))
dev.off()

#tiff full dataset shown
par(mar=c(7,4,4,2)+0.1)
#tiff('Heatmap_species_flavours_1.tiff', units="in", width=14, height=12, res=300, compression = 'lzw')
tiff('Heatmap_avgs_species_newsp_rname+Ljp.tiff', units="in", width=14, height=12, res=300, compression = 'lzw')
par(mfrow = c(3, 3),  mar=c(5,4.5,4,2))
colfunc <- colorRampPalette(c("blue", "yellow"))
heatmap(ions_sc,col=colfunc(15),scale="column", cexCol=1.2,margins=c(12,8))
#heatmap(ions_sc,col=colfunc(15),scale="column",cexCol=0.9,margins=c(12,8))
dev.off()

#install.packages("ComplexHeatmap") not for this version r
#library(ComplexHeatmap)

#try to split
aas1 <- gr[7:20]
aas2 <- gr[28:33]
sm1 <- gr[21:27]
sm2 <- gr[34:36]
su <- gr[37:40]
sm_hplc <- gr[42:49]

puttog <- cbind(aas1, aas2, sug, sm1, sm2, sm_hplc, gr$Species)

#reorder species
sum2 <- read.csv("Met_summ_sp_newsp_Ljp_renames.csv")
#try to split
aas1 <- sum2[2:15] #correct
aas2 <- sum2[23:28] #correct
sm1 <- sum2[16:22] #correct
sm2 <- sum2[29:31] # correct
su <- sum2[32:35] #correct
sm_hplc <- sum2[36:43] #correct

puttog <- cbind(aas1, aas2, su, sm1, sm2, sm_hplc, sum2$Species)


#avgs sds per species
colnames(puttog)[43] ="Species"
puttog$Species <- as.factor(puttog$Species)
#sum2 <- puttog %>% group_by(Species) %>% summarise_all(mean)
#sum2
colnames(puttog)[26] <- "Chlorogenic.acid.Isomer"


#sumx <- sum2[-1]
sumx <- puttog[-43]

ions_sc <- scale(sumx,center=T,scale=T)
ions_sc <- as.matrix(ions_sc)
#assign row names from df
rownames(ions_sc) <- puttog$Species


#reordered data
par(mar=c(7,4,4,2)+0.1)
#tiff('Heatmap_species_flavours_1.tiff', units="in", width=14, height=12, res=300, compression = 'lzw')
tiff('Heatmap_avgs_species_newsp_rname_rord.tiff', units="in", width=14, height=12, res=300, compression = 'lzw')
par(mfrow = c(3, 3),  mar=c(5,4.5,4,2))
colfunc <- colorRampPalette(c("blue", "yellow"))
heatmap(ions_sc,col=colfunc(15),scale="column", Colv = NA, cexCol=1.2,margins=c(12,8))
#heatmap(ions_sc,col=colfunc(15),scale="column",cexCol=0.9,margins=c(12,8))
dev.off()

#improve heatmap with grouping gaps + side scale
#install.packages("pheatmap")
library(pheatmap)
display.brewer.all(colorblindFriendly = TRUE)

cols <- 

pheatmap(ions_sc,
         color = hcl.colors(75, "BluYl"),
         scale='column',
         gaps_row=c(1,3),
         gaps_col=c(3,6,9),
         #annotation_row = rownames(ions_sc),
         annotation_names_row=F,
         #cluster_rows = TRUE,
         #cluster_cols = TRUE,
         show_colnames = TRUE,
         show_rownames = TRUE)

#still basic try this:
#not available for R version
#install.packages("ComplexHeatmap")
#library(ComplexHeatmap)

#pca
ions_all <- hl_llmeans[,27:47]

#using raw data
ions_all <- hl_llmeans[,3:45]
#define factors
ions_access <- hl_llmeans[,2]
ions_species <- hl_llmeans[,5]
str(ions_all)

library(FactoMineR)
ions_acc <- ions_all[ ,c(1:43)] # selecting columns from csv
ions_acc.pca <- PCA(ions_acc, quali.sup=42) 
print(ions_acc.pca)
head(ions_acc.pca)
print(summary(ions_acc.pca)) #s

#plot cos 2 as bar graph #high = good representation on pc
library(FactoMineR)
library(factoextra)

#basic plot
biplot(ions_acc.pca) #not working

fviz_pca_biplot(ions_acc.pca, repel=TRUE, pointsize=6, pointshape=21, col.var="black", arrowsize=0.6, labelsize=5, col.ind=ions_species, palette=c("green2", "gold", "blue", "purple"), addEllipses=TRUE, ellipse.type="confidence")

fviz_cos2(ions_acc.pca, choice = "var", axes = 1:2) #K and S biggest
#top phenylalanine, sec mets

fviz_pca_var(ions_acc.pca, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
)

ions_acc.pca<-prcomp(ions_acc[ ,1:43],center=T,scale=T)
str(ions_acc.pca)
mypc <- ions_acc.pca$x #define new variable x = pcs just need to plot these on xy graph

#plot accession, cant tell which which accession n =25 but all
#reps so not colored same for each rep
dev.off()
plot(mypc[,1], mypc[,2], col = hl_llmeans$Species,
     las=1, xlab="PC1 (32%)", ylab="PC2 (19%)",
     pch=16, cex=1.5, xlim=c(-10,15), ylim=c(-10,15)) #pcs as columns, produce xy plot, las is rotation of axis numbers, pch plot shape, ylim expand out so legend room
abline(v=0, lty=2, col="lightgrey") #draw line, lty is segmented
abline(h=0, lty=2, col="lightgrey") #0 lines dashed
#legend("right", pch=16, col=my_pal, cex=1, c("KS02", "KS03", "KS04", "KS06A", "KS06B",
#                                             "KS09", "KS12", "KS13", "KS14", "KS15",
#                                             "KS16", "KS17","KS18", "KS20", "KS21", "KS22", "KS25",
#                                             "KS28", "KS29", "KS66A", "KS77A", "KS78A",
#                                             "LY01A", "LY01B", "Nuff1")) #customise legend seperately cex=txtsize cols 1-4 stnd, concat order of places as want to display
#too busy #can show outliers
#text(x=mypc[,1],y=mypc[,2], labels =ions_access, pos=2) #txt and pos to define labels on plot, 2 = top of plot, 1 bottom

#need more than 12 colors as repeating
my_pal <- scico::scico(length(unique(ions_species)), palette = "batlow")
my_pal <- scico::scico(length(unique(ions_access)), palette = "lisbon")

#species all raw data
plot(mypc[,1], mypc[,2], col = hl_llmeans$Species,
     las=1, xlab="PC1 (32%)", ylab="PC2 (19%)",
     pch=16, cex=1.5, xlim=c(-10,10), ylim=c(-10,10)) #pcs as columns, produce xy plot, las is rotation of axis numbers, pch plot shape, ylim expand out so legend room
abline(v=0, lty=2, col="lightgrey") #draw line, lty is segmented
abline(h=0, lty=2, col="lightgrey") #0 lines dashed
legend("right", pch=16, col=unique(ions_species), cex=1, c("L. minor", "L. minuta", "S. polyrhiza", "L. turionifera")) #customise legend seperately cex=txtsize cols 1-4 stnd, concat order of places as want to display
#to show if coloring/grouping properly
text(x=mypc[,1],y=mypc[,2], labels =ions_access, pos=2) #txt and pos to define labels on plot, 2 = top of plot, 1 bottom
#dev.new()

#run pca with summary data so only 25 points
ions_acc <- sum[ ,c(2:43)] # selecting columns from csv
ions_acc.pca <- PCA(ions_acc, quali.sup=42) #cant exceed max col number
print(ions_acc.pca)
head(ions_acc.pca)
print(summary(ions_acc.pca))

#new sp
sum_species <- c("L. minor", "L. minor", "L. minor", "L. minuta",
                 "L. minor", "S. polyrhiza", "L. minor", "L. minor", "L. minor",
                 "L. turionifera", "L. minor","L. minor", "L. minuta", "L. minor", "L. turionifera", "L. minuta",
                 "L. minor", "L. minor", "L. minor", "S. polyrhiza", "S. polyrhiza",
                 "L. minor", "L. minuta", "L. minor") 

as.factor(sum_species)
#too busy #can show outliers

biplot(ions_acc.pca) #not working

fviz_pca_biplot(ions_acc.pca, repel=TRUE, pointsize=6, pointshape=21, col.var="black", arrowsize=0.6, labelsize=5, col.ind=sum_species, palette=c("green2", "gold", "blue", "purple"), addEllipses=TRUE, ellipse.type="confidence")
#not working as accessions

fviz_cos2(ions_acc.pca, choice = "var", axes = 1:2) #K and S biggest
#sec mets, proline, isoleuc

fviz_pca_var(ions_acc.pca, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
)

ions_acc.pca<-prcomp(ions_acc[ ,1:42],center=T,scale=T) #list can include mix mat and df
str(ions_acc.pca)
mypc <- ions_acc.pca$x #define new variable x = pcs just need to plot these on xy graph

library(vegan) #needed for ellipses
#plot species, can see some species seperation
plot(mypc[,1], mypc[,2], col = factor(hl_llmeans$Species),
     las=1, xlab="PC1 (26%)", ylab="PC2 (15%)",
     pch=16, cex=1.5, xlim=c(-15,15), ylim=c(-15,15)) #pcs as columns, produce xy plot, las is rotation of axis numbers, pch plot shape, ylim expand out so legend room
abline(v=0, lty=2, col="lightgrey") #draw line, lty is segmented
abline(h=0, lty=2, col="lightgrey") #0 lines dashed
legend("right", pch=16, col=unique(factor(hl_llmeans$Specises)), cex=1, c("L. minor", "L. minuta", "S. polyrhiza", "L. turionifera")) #customise legend seperately cex=txtsize cols 1-4 stnd, concat order of places as want to display
#too busy #can show outliers
text(x=mypc[,1],y=mypc[,2], labels =hl_llmeans$Accession, pos=2) #txt and pos to define labels on plot, 2 = top of plot, 1 bottom
#l minor and l minuta apart, l gibba most varied
#ks03 quite near ks12? ly01a near potential l gibbas
ordiellipse(ions_acc.pca,hl_llmeans$Species,conf=0.95)
#dev.new()

#attempt ggplot pca
#best one for raw data colored by sp
fviz_pca_ind(ions_acc.pca, geom.ind = "point",
             label = "all",
             pointshape = 21,
             pointsize = 3, 
             fill=hl_llmeans$Species,
             addEllipses=TRUE, legend.title = "Species", title="")+
  scale_color_manual(values=c("red", "blue", "darkgreen", "purple"))+
  theme_classic()
ggsave("PCA_Metabolites_spnew_col.tiff", units="in", width=5, height=4, dpi=300, compression = 'lzw')


#with ggplot
#install.packages("ggforce") #used to plot all elipses <3 points
library(ggforce)
ions_acc.pca
PC1<-ions_acc.pca$x[,1]
PC2<-ions_acc.pca$x[,2]
ggplot(hl_llmeans, 
       aes(x = PC1, 
           y = PC2, 
           color = ions_species,
           repel = TRUE),
       invisible='quali') +
  geom_point(size=5) +
  # Don't use default Ellipses!!!!
  # addEllipses = TRUE,
  # ADD ggforce's ellipses
  ggforce::geom_mark_ellipse(aes(fill = ions_species,
                                 color = ions_species)) +
  theme(legend.position = 'right') +
  theme_classic()
coord_equal()

#just for l minor l minuta elipses
pca <- ggplot(sum, aes(x = PC1, y = PC2, color = sum_species)) +
  geom_point(size=5) +
  ggforce::geom_mark_ellipse(aes(fill = sum_species,
                                 color = sum_species)) +
  theme_classic() +
  geom_hline(yintercept=0, color="gray", linetype="dashed") + 
  geom_vline(xintercept=0, color="gray", linetype="dashed") +
  scale_x_continuous(limits = c(-15, 15), name="PC1 (31%)") +
  scale_y_continuous(limits = c(-15, 15), name="PC2 (22%)") +
  theme(legend.title=element_blank()) +
  theme(legend.text = element_text(face = "italic"))
ggsave("pca_spnew.tiff", dpi = 300, width = 17, height = 15 , units = "cm") 

shapiro.test(gr$L.Alanine)
#normality tests shapiro wilk not normal

#summary data
hl_llmeans <- read.csv("Met_summ_acc_newsp.csv")
hl_llmeans <- read.csv("Met_summ_sp_newsp.csv")
hl_llmeans <- read.csv("Met_summ_sp_newsp_Ljp_renames.csv")

#not summary data
hl_llmeans <- read.csv("Metabolite_all_spnew+Ljp.csv")
#raw
hl_llmeans <- gr[c(1,4,7:49)]

shapiro.test(hl_llmeans$L.Alanine)
#normal
shapiro.test(hl_llmeans$L.Arginine)
#sig so not normal
shapiro.test(hl_llmeans$L.Asparagine)
#normal
shapiro.test(hl_llmeans$L.Aspartic.Acid)
#sig so not normal
shapiro.test(hl_llmeans$L.Glutamic.acid)
#sig so not normal
shapiro.test(hl_llmeans$L.Glutamine)
#normal
shapiro.test(hl_llmeans$L.Isoleucine_1)
#normal
shapiro.test(hl_llmeans$L.Leucine_2)
#normal
shapiro.test(hl_llmeans$L.Phenylalanine)
#normal
shapiro.test(hl_llmeans$L.Proline)
#normal
shapiro.test(hl_llmeans$L.Serine)
#normal
shapiro.test(hl_llmeans$L.Threonine)
#normal
shapiro.test(hl_llmeans$L.Tryptophane)
#normal
shapiro.test(hl_llmeans$L.Valine)
#not normal
shapiro.test(hl_llmeans$Chlorogenic.acid)
#not normal
shapiro.test(hl_llmeans$putative.Chlorogenic.acid.Isomer)
#doesnt work
shapiro.test(hl_llmeans$Cyanidine.3.Glc)
#doesnt work
shapiro.test(hl_llmeans$Luteoline.8.C.Glc)
#not normal
shapiro.test(hl_llmeans$Apigenine.8.C.Glc)
#not normal
shapiro.test(hl_llmeans$Luteoline.7.O.Glc)
#not normal
shapiro.test(hl_llmeans$Apigenine.7.O.Glc)
#not normal
shapiro.test(hl_llmeans$Glycine)
#normal
shapiro.test(hl_llmeans$L.Histidine)
shapiro.test(hl_llmeans$L.Methionine)
shapiro.test(hl_llmeans$L.Tyrosine)
shapiro.test(hl_llmeans$Tyramine)
shapiro.test(hl_llmeans$Shikimic.acid)
#all normal

shapiro.test(hl_llmeans$Apigenine)
shapiro.test(hl_llmeans$Luteoline)
shapiro.test(hl_llmeans$Glucose)
shapiro.test(hl_llmeans$Fructose)
shapiro.test(hl_llmeans$Sucrose)
#all not normal

shapiro.test(hl_llmeans$Starch..Glu.)
shapiro.test(hl_llmeans$Starch.mg...g)
#normal

shapiro.test(hl_llmeans$Cyanidine.3.Glc.1)
shapiro.test(hl_llmeans$Cyanidin.Mal.Glc)
shapiro.test(hl_llmeans$Chlorogenic.acid.1)
shapiro.test(hl_llmeans$Luteolin.7.O.Glc)
#all not normal 
shapiro.test(hl_llmeans$Luteolin.8.C.Glc)
shapiro.test(hl_llmeans$Luteolin.8.C.Glc)
#normal

shapiro.test(hl_llmeans$Apigenin.7.O.Glc)
shapiro.test(hl_llmeans$Apigenin.8.C.Glc)
shapiro.test(hl_llmeans$Luteolin)
shapiro.test(hl_llmeans$Apigenin)
#not normal

#for accession dataset
#count normal 20
#count not normal 22
#generally a.as more normal, sec mets not normal
#either split analysis or use non parametric test

#all not normal when use raw dataset

#kruskal wallis accessions
kruskal.test(L.Alanine ~ Accession, data = gr)
#sig
#pairwise wilcox test a.as
pairwise.wilcox.test(gr$L.Alanine, gr$Accession,
                     p.adjust.method = "BH")
pairwise.wilcox.test(gr$L.Arginine, gr$Accession,
                     p.adjust.method = "BH")
pairwise.wilcox.test(gr$L.Asparagine, gr$Accession,
                     p.adjust.method = "BH")
pairwise.wilcox.test(gr$L.Aspartic.Acid, gr$Accession,
                     p.adjust.method = "BH")
pairwise.wilcox.test(gr$L.Glutamine, gr$Accession,
                     p.adjust.method = "BH")
pairwise.wilcox.test(gr$L.Isoleucine_1, gr$Accession,
                     p.adjust.method = "BH")
pairwise.wilcox.test(gr$L.Leucine_2, gr$Accession,
                     p.adjust.method = "BH")
pairwise.wilcox.test(gr$Phenylalanine, gr$Accession,
                     p.adjust.method = "BH")
pairwise.wilcox.test(gr$Proline, gr$Accession,
                     p.adjust.method = "BH")
pairwise.wilcox.test(gr$L.Serine, gr$Accession,
                     p.adjust.method = "BH")
pairwise.wilcox.test(gr$L.Threonine, gr$Accession,
                     p.adjust.method = "BH")
pairwise.wilcox.test(gr$L.Tryptophane, gr$Accession,
                     p.adjust.method = "BH")
pairwise.wilcox.test(gr$L.Valine, gr$Accession,
                     p.adjust.method = "BH")

#kruskal wallis raw data
kruskal.test(L.Alanine ~ Species, data = hl_llmeans)
#n.s
kruskal.test(L.Arginine~Species, data=hl_llmeans)
 #0.01
kruskal.test(L.Asparagine~Species, data=hl_llmeans)
 #n.s
kruskal.test(L.Aspartic.Acid~Species, data=hl_llmeans)
 #0.05
kruskal.test(L.Glutamic.acid~Species, data=hl_llmeans)
 #sp 0.005 
kruskal.test(L.Glutamine~Species, data=hl_llmeans)
 #n.s
kruskal.test(L.Isoleucine_1~Species, data=hl_llmeans)
 #n.s
kruskal.test(L.Leucine_2~Species, data=hl_llmeans)
 #n.s
kruskal.test(L.Phenylalanine~Species, data=hl_llmeans)
 #n.s
kruskal.test(L.Proline~Species, data=hl_llmeans)
 #n.s
kruskal.test(L.Serine~Species, data=hl_llmeans)
 #n.s
kruskal.test(L.Threonine~Species, data=hl_llmeans)
 #n.s
kruskal.test(L.Tryptophane~Species, data=hl_llmeans)
 #sp 0.0001 
kruskal.test(L.Valine~Species, data=hl_llmeans)
 #n.s
kruskal.test(Chlorogenic.acid~Species, data=hl_llmeans)
 #<0.0000001
kruskal.test(putative.Chlorogenic.acid.Isomer~Species, data=hl_llmeans)
 #<0.0000001
kruskal.test(Cyanidine.3.Glc~Species, data=hl_llmeans)
 #*<0.0000001
kruskal.test(Luteoline.8.C.Glc~Species, data=hl_llmeans)
 #<0.0000001
kruskal.test(Apigenine.8.C.Glc~Species, data=hl_llmeans)
 #<0.0000001
kruskal.test(Luteoline.7.O.Glc~Species, data=hl_llmeans)
 #<0.0000001
kruskal.test(Apigenine.7.O.Glc~Species, data=hl_llmeans)
 #<0.0000001
kruskal.test(Glycine~Species, data=hl_llmeans)
 #n.s
kruskal.test(L.Histidine~Species, data=hl_llmeans)
 #0.01 sp
kruskal.test(L.Methionine~Species, data=hl_llmeans)
 #0.01 sp
kruskal.test(L.Tyrosine~Species, data=hl_llmeans)
 #0.003 sp
kruskal.test(Tyramine~Species, data=hl_llmeans)
 #0.000005 
kruskal.test(Shikimic.acid~Species, data=hl_llmeans)
 #<0.0000001
kruskal.test(Apigenine~Species, data=hl_llmeans)
 #<0.0000001
kruskal.test(Luteoline~Species, data=hl_llmeans)
 #*** species, ns acc
kruskal.test(Glucose~Species, data=hl_llmeans)
 #ns
kruskal.test(Fructose~Species, data=hl_llmeans)
 #ns
kruskal.test(Sucrose~Species, data=hl_llmeans)
 #ns
kruskal.test(Starch..Glu.~Species, data=hl_llmeans)
 #ns
kruskal.test(Starch.mg...g~Species, data=hl_llmeans)
 #ns
kruskal.test(Cyanidine.3.Glc.1~Species, data=hl_llmeans)
 #<0.0000001 sp
kruskal.test(Cyanidin.Mal.Glc~Species, data=hl_llmeans)
 #<0.0000001 sp
kruskal.test(Chlorogenic.acid.1~Species, data=hl_llmeans)
 #sp 0.03
kruskal.test(Luteolin.7.O.Glc~Species, data=hl_llmeans)
 #sp 0.000005
kruskal.test(Luteolin.8.C.Glc~Species, data=hl_llmeans)
 #sp 0.0002
kruskal.test(Luteolin.8.C.Glc~Species, data=hl_llmeans)
 #sp 0.0002
kruskal.test(Apigenin.7.O.Glc~Species, data=hl_llmeans)
 #<0.0000001
kruskal.test(Apigenin.8.C.Glc~Species, data=hl_llmeans)
 #<0.0000001
kruskal.test(Luteolin~Species, data=hl_llmeans)
 #<0.0000001
kruskal.test(Apigenin~Species, data=hl_llmeans)
 #<0.0000001

#a.as barely different, sugars not different, sec mets diff
#sig a.as: = 8
#L.Arginine, Glutamic.acid, L.Tryptophane,
#L.Histidine, L.Methionine, L.Tyrosine, Tyramine, Shikimic.acid

#sig a.as raw data with Ljp grouping species = 8 (or 9 if inc shik)
#tyramine most different, tryptophan, glutamic acid, tyrosine, methionine, histidine, arginine, asp acid
#different: include aspartic acid, shikimic class as sec met or precursor?

#essential: histidine, methionine, tryptophan
#non essential: tyramine (derived from tyrosine), glutamic acid, tryrosine, aspartic acid, arginine

#pairwise wilcox test a.as
pairwise.wilcox.test(gr$L.Arginine, gr$Species,
                     p.adjust.method = "BH")
#l turion diff l.minuta, almost sig diff from all others
pairwise.wilcox.test(gr$L.Glutamic.acid, gr$Species,
                     p.adjust.method = "BH")
#l.turion diff to l.minu, s. poly diff from l. minu
pairwise.wilcox.test(gr$L.Tryptophane, gr$Species,
                     p.adjust.method = "BH")
#l mino diff to ljp, l minu diff to lmino and ljp, s. poly diff to l.mino and ljp,
pairwise.wilcox.test(gr$L.Histidine, gr$Species,
                     p.adjust.method = "BH")
#L minu diff L. minor, S poly diff L. mino
pairwise.wilcox.test(gr$L.Methionine, gr$Species,
                     p.adjust.method = "BH")
#ns
pairwise.wilcox.test(gr$L.Tyrosine, gr$Species,
                     p.adjust.method = "BH")
#s poly diff all
pairwise.wilcox.test(gr$Tyramine, gr$Species,
                     p.adjust.method = "BH")
#l minu diff L. mino Ljp, L turion diff L. mino Ljp, S poly diff L. minu and L. turion
pairwise.wilcox.test(gr$L.Aspartic.Acid, gr$Species,
                     p.adjust.method = "BH")
#ns
pairwise.wilcox.test(gr$Shikimic.acid, gr$Species,
                     p.adjust.method = "BH")
#l minu diff l mino ljp, l turion diff l minu, s poly diff all

#aspartic acid and methionine not diff by species

#change order of species groups
gr$Species <- factor(gr$Species, levels = c("L. minor", "L. japonica", "L. turionifera", "L. minuta","S. polyrhiza"))
gr$Species

#7 a.as sig to plot as boxplots
#tiff('AAs_species_boxplots_2_KWtest_spnew_newcol_Ljp.tiff', units="in", width=14, height=12, res=300, compression = 'lzw')
pdf('AAs_species_boxplots_2_KWtest_spnew_newcol_Ljp.pdf', width=14, height=12)
par(mfrow = c(3,3),  mar=c(5,4.5,4,2))
boxplot(L.Tryptophane~Species,data=gr, col = c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        main= ("Essential amino acids"), ylab = expression(paste("Tryptophan  (", mu, "mol g/DW)")), xlab =  "", las=2,
        ylim=c(0,12),
        cex.lab=1.5, cex.axis=1.5)#***
text(1, 11.5, "a", cex=1.5)
text(2, 11.5, "b", cex=1.5)
text(3, 11.5, "abc", cex=1.5)
text(4, 11.5, "c", cex=1.5)
text(5, 11.5, "c", cex=1.5)
#l mino diff to ljp, l minu diff to lmino and ljp, s. poly diff to l.mino and ljp,
boxplot(L.Histidine~Species,data=gr, col = c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        ylab = expression(paste("Histidine  (", mu, "mol g/DW)")), xlab =  "", las=2,
        ylim=c(0,5),
        cex.lab=1.5, cex.axis=1.5)#***
text(1, 4.8, "a", cex=1.5)
text(2, 4.8, "ab", cex=1.5)
text(3, 4.8, "ab", cex=1.5)
text(4, 4.8, "b", cex=1.5)
text(5, 4.8, "b", cex=1.5)
#boxplot()
#boxplot()
boxplot(L.Arginine~Species,data=gr, col = c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        ylab = expression(paste("Arginine  (", mu, "mol g/DW)")), xlab =  "", las=2,
        ylim=c(0,150),
        cex.lab=1.5, cex.axis=1.5)#***
text(1, 148, "ab", cex=1.5)
text(2, 148, "ab", cex=1.5)
text(3, 148, "a", cex=1.5)
text(4, 148, "b", cex=1.5)
text(5, 148, "ab", cex=1.5)
boxplot(L.Glutamic.acid~Species,data=gr, col = c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        ylab = expression(paste("Glutamic acid  (", mu, "mol g/DW)")), xlab =  "", las=2,
        ylim=c(0,60),
        cex.lab=1.5, cex.axis=1.5)#***
text(1, 58, "ab", cex=1.5)
text(2, 58, "ab", cex=1.5)
text(3, 58, "a", cex=1.5)
text(4, 58, "b", cex=1.5)
text(5, 58, "a", cex=1.5)
boxplot(L.Tyrosine~Species,data=gr, col = c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        ylab = expression(paste("Tyrosine  (", mu, "mol g/DW)")), xlab =  "", xaxt="n",  las=2,
        ylim=c(0,12),
        cex.lab=1.5, cex.axis=1.5)#***
text(1, 11.8, "a", cex=1.5)
text(2, 11.8, "a", cex=1.5)
text(3, 11.8, "a", cex=1.5)
text(4, 11.8, "a", cex=1.5)
text(5, 11.8, "b", cex=1.5)
boxplot(Tyramine~Species,data=gr, col = c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        ylab = expression(paste("Tyramine  (", mu, "mol g/DW)")), xlab =  "", xaxt="n",  las=2,
        ylim=c(0,0.5),
        cex.lab=1.5, cex.axis=1.5)#***
text(1, 0.48, "a", cex=1.5)
text(2, 0.48, "a", cex=1.5)
text(3, 0.48, "b", cex=1.5)
text(4, 0.48, "b", cex=1.5)
text(5, 0.48, "a", cex=1.5)
boxplot(Shikimic.acid~Species,data=gr, col = c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        ylab = expression(paste("Shikimic acid  (", mu, "AU g/DW)")), xlab =  "", xaxt="n", las=2,
        ylim=c(0,3),
        cex.lab=1.5, cex.axis=1.5)#***
text(1, 2.9, "b", cex=1.5)
text(2, 2.9, "b", cex=1.5)
text(3, 2.9, "b", cex=1.5)
text(4, 2.9, "c", cex=1.5)
text(5, 2.9, "a", cex=1.5)
#l minu diff l mino ljp, l turion diff l minu, s poly diff all
dev.off()

#Sec metabolites
#lcms data
pairwise.wilcox.test(gr$Chlorogenic.acid, gr$Species,
                     p.adjust.method = "BH")
#s poly diff all
pairwise.wilcox.test(gr$putative.Chlorogenic.acid.Isomer, gr$Species,
                     p.adjust.method = "BH")
#s poly diff ljp
pairwise.wilcox.test(gr$Cyanidine.3.Glc, gr$Species,
                     p.adjust.method = "BH")
#s poly diff all
pairwise.wilcox.test(gr$Luteoline.8.C.Glc, gr$Species,
                     p.adjust.method = "BH")
#s poly diff l minor ljp l minu, l. turion diff l.mino and ljp, lminu diff ltu, l mino diff ljp
pairwise.wilcox.test(gr$Apigenine.8.C.Glc, gr$Species,
                     p.adjust.method = "BH")
#s poly diff all, l. turion diff l.mino ljp, minu diff ljp ltur
pairwise.wilcox.test(gr$Luteoline.7.O.Glc, gr$Species,
                     p.adjust.method = "BH")
#s poly diff all, l. minuta diff all, l mino l jp diff
pairwise.wilcox.test(gr$Apigenine.7.O.Glc, gr$Species,
                     p.adjust.method = "BH")
#s poly diff all
pairwise.wilcox.test(gr$Apigenine, gr$Species,
                     p.adjust.method = "BH")
#s poly diff all
pairwise.wilcox.test(gr$Luteoline, gr$Species,
                     p.adjust.method = "BH")
#s poly diff all, ltu diff to lmo ljp, l minu diff to l mino ljp

#Chlorogenic.acid, putative.Chlorogenic.acid.Isomer, 
#Cyanidine.3.Glc, Apigenine.8.C.Glc, 
#Luteoline.7.O.Glc, Apigenine.7.O.Glc, Luteoline.8.C.Glc

#9 a.as sig to plot as boxplots
#tiff('SecMetLCMS_species_boxplots_2_KWtest_spnew_newcol+Ljp.tiff', units="in", width=14, height=12, res=300, compression = 'lzw')
pdf('SecMetLCMS_species_boxplots_2_KWtest_spnew_newcol+Ljp.pdf', width=14, height=12)
par(mfrow = c(3,3),  mar=c(5,4.5,4,2))
boxplot(Chlorogenic.acid~Species,data=gr, col = c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        main="Secondary metabolites (LCMS)", cex.main=2, ylab = "Chlorogenic.acid (AU /g)", xlab =  "", las=2,
        ylim=c(0,1),
        cex.lab=1.5, cex.axis=1.5)#***
text(1, 0.95, "b", cex=1.5)
text(2, 0.95, "b", cex=1.5)
text(3, 0.95, "b", cex=1.5)
text(4, 0.95, "b", cex=1.5)
text(5, 0.95, "a", cex=1.5)
#s poly diff all
boxplot(putative.Chlorogenic.acid.Isomer~Species,data=gr, col = c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        ylab = "Chlorogenic acid isomer (AU /g)", xlab =  "", las=2,
        ylim=c(0,1.5),
        cex.lab=1.5, cex.axis=1.5)#***
text(1, 1.5, "ab", cex=1.5)
text(2, 1.5, "b", cex=1.5)
text(3, 1.5, "ab", cex=1.5)
text(4, 1.5, "ab", cex=1.5)
text(5, 1.5, "a", cex=1.5)
#s poly diff ljp
boxplot(Cyanidine.3.Glc~Species,data=gr, col = c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        ylab = "Cyanidin-3-Glucoside (AU/ g)", xlab =  "", las=2,
        ylim=c(0,11),
        cex.lab=1.5, cex.axis=1.5)#***
text(1, 10.8, "b", cex=1.5)
text(2, 10.8, "b", cex=1.5)
text(3, 10.8, "b", cex=1.5)
text(4, 10.8, "b", cex=1.5)
text(5, 10.8, "a", cex=1.5)
boxplot(Luteoline~Species,data=gr, col = c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        ylab = "Luteolin (AU/ g)", xlab =  "", las=2,
        ylim=c(0,13),
        cex.lab=1.5, cex.axis=1.5)#***
text(1, 12.8, "c", cex=1.5)
text(2, 12.8, "c", cex=1.5)
text(3, 12.8, "b", cex=1.5)
text(4, 12.8, "b", cex=1.5)
text(5, 12.8, "a", cex=1.5)
#s poly diff all, ltu diff to lmo ljp, l minu diff to l mino ljp
boxplot(Apigenine~Species,data=gr, col = c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        ylab = "Apigenin (AU/ g)", xlab =  "", las=2,
        ylim=c(0,2),
        cex.lab=1.5, cex.axis=1.5)#***
text(1, 1.9, "b", cex=1.5)
text(2, 1.9, "b", cex=1.5)
text(3, 1.9, "b", cex=1.5)
text(4, 1.9, "b", cex=1.5)
text(5, 1.9, "a", cex=1.5)
#s poly diff all
boxplot(Luteoline.7.O.Glc~Species,data=gr, col = c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        ylab = "Luteolin-7-O-Glucoside (AU/ g)", xlab =  "", las=2,
        ylim=c(0,40),
        cex.lab=1.5, cex.axis=1.5)#***
text(1, 39.6, "e", cex=1.5)
text(2, 39.6, "d", cex=1.5)
text(3, 39.6, "cde", cex=1.5)
text(4, 39.6, "b", cex=1.5)
text(5, 39.6, "a", cex=1.5)
#s poly diff all, l. minuta diff all, l mino l jp diff
boxplot(Apigenine.7.O.Glc~Species,data=gr, col = c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        ylab = "Apigenin-7-O-Glucoside (AU/ g)", xlab =  "", xaxt = "n", las=2,
        ylim=c(0,10),
        cex.lab=1.5, cex.axis=1.5)#***
text(1, 9.8, "b", cex=1.5)
text(2, 9.8, "b", cex=1.5)
text(3, 9.8, "b", cex=1.5)
text(4, 9.8, "b", cex=1.5)
text(5, 9.8, "a", cex=1.5)
#s poly diff all
boxplot(Luteoline.8.C.Glc~Species,data=gr, col = c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        ylab = "Luteolin-8-C-Glucoside (AU/ g)", xlab =  "", xaxt = "n", las=2,
        ylim=c(0,42),
        cex.lab=1.5, cex.axis=1.5)#***
text(1, 39.6, "d", cex=1.5)
text(2, 39.6, "c", cex=1.5)
text(3, 39.6, "ab", cex=1.5)
text(4, 39.6, "cd", cex=1.5)
text(5, 39.6, "a", cex=1.5)
#s poly diff l minor ljp l minu, l. turion diff l.mino and ljp, lminu diff ltu, l mino diff ljp
boxplot(Apigenine.8.C.Glc~Species,data=gr, col = c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        ylab = "Apigenin-8-C-Glucoside (AU/ g)", xlab =  "", xaxt = "n", las=2,
        ylim=c(0,15),
        cex.lab=1.5, cex.axis=1.5)#***
text(1, 14.8, "cd", cex=1.5)
text(2, 14.8, "cd", cex=1.5)
text(3, 14.8, "a", cex=1.5)
text(4, 14.8, "ad", cex=1.5)
text(5, 14.8, "b", cex=1.5)
#s poly diff all, l. turion diff l.mino ljp, minu diff ljp ltur
dev.off()

#Sec mets
#HPLC data
#Sec metabolites
#lcms data
pairwise.wilcox.test(gr$Cyanidine.3.Glc.1, gr$Species,
                     p.adjust.method = "BH")
#s poly diff all, l turion diff l. jp
pairwise.wilcox.test(gr$Cyanidin.Mal.Glc, gr$Species,
                     p.adjust.method = "BH")
#s poly diff all, l turion diff to l jp
pairwise.wilcox.test(gr$Chlorogenic.acid.1, gr$Species,
                     p.adjust.method = "BH")
#l turion nearly diff from l minor l jap
pairwise.wilcox.test(gr$Luteolin.7.O.Glc, gr$Species,
                     p.adjust.method = "BH")
#s poly diff all, l. turion diff l. minu
pairwise.wilcox.test(gr$Luteolin.8.C.Glc, gr$Species,
                     p.adjust.method = "BH")
#l minu diff l. minor l. jap and l. turion, s poly
pairwise.wilcox.test(gr$Apigenin.7.O.Glc, gr$Species,
                     p.adjust.method = "BH")
#l turion diff l jp l mino, l minu diff all, s poly diff all
pairwise.wilcox.test(gr$Apigenine.8.C.Glc, gr$Species,
                     p.adjust.method = "BH")
#s poly diff all, l. turion diff to ljp lmino, l minu diff to l jp l turi, l mino diff l jp
pairwise.wilcox.test(gr$Luteolin, gr$Species,
                     p.adjust.method = "BH")
#l turi diff lmino ljp, l minu diff to lmino ljp, a poly diff to all but ltu
pairwise.wilcox.test(gr$Apigenin, gr$Species,
                     p.adjust.method = "BH")
#s poly diff all

#Cyanidine.3.Glc.1,Cyanidin.Mal.Glc,Chlorogenic.acid.1,
#Luteolin.7.O.Glc,Luteolin.8.C.Glc,Apigenin.7.O.Glc,
#Apigenin.8.C.Glc,Luteolin,Apigenin

#9 sec mets by HPLC
#tiff('SecMetHPLC_species_boxplots_2_KWtest_spnew_newcol+Ljp.tiff', units="in", width=14, height=12, res=300, compression = 'lzw')
pdf('SecMetHPLC_species_boxplots_2_KWtest_spnew_newcol+Ljp.pdf', width=14, height=12)
par(mfrow = c(3,3),  mar=c(5,4.5,4,2))
boxplot(Cyanidine.3.Glc.1~Species,data=gr, col = c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        main="Secondary metabolites (HPLC)", cex.main=2, ylab = expression(paste("Cyanidin-3-Glucoside  (", mu, "mol g/DW)")), xlab =  "", las=2,
        ylim=c(0,5),
        cex.lab=1.5, cex.axis=1.5)#***
text(1, 4.8, "bc", cex=1.5)
text(2, 4.8, "c", cex=1.5)
text(3, 4.8, "b", cex=1.5)
text(4, 4.8, "bc", cex=1.5)
text(5, 4.8, "a", cex=1.5)
#s poly diff all, l turion diff l. jp
boxplot(Cyanidin.Mal.Glc~Species,data=gr, col = c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        ylab = expression(paste("Cyanidin-malonyl glucoside  (", mu, "mol g/DW)")), xlab =  "", las=2,
        ylim=c(0,10),
        cex.lab=1.5, cex.axis=1.5)#***
text(1, 9.8, "bc", cex=1.5)
text(2, 9.8, "c", cex=1.5)
text(3, 9.8, "b", cex=1.5)
text(4, 9.8, "bc", cex=1.5)
text(5, 9.8, "a", cex=1.5)
boxplot(Luteolin~Species,data=gr, col = c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        ylab = expression(paste("Luteolin  (", mu, "mol g/DW)")), xlab =  "", las=2,
        ylim=c(0,5),
        cex.lab=1.5, cex.axis=1.5)#***
text(1, 4.8, "c", cex=1.5)
text(2, 4.8, "c", cex=1.5)
text(3, 4.8, "ab", cex=1.5)
text(4, 4.8, "ab", cex=1.5)
text(5, 4.8, "ab", cex=1.5)
#l turi diff lmino ljp, l minu diff to lmino ljp, a poly diff to all but ltu
boxplot(Apigenin~Species,data=gr, col = c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        ylab = expression(paste("Apigenin  (", mu, "mol g/DW)")), xlab =  "", las=2,
        ylim=c(0,2),
        cex.lab=1.5, cex.axis=1.5)#***
text(1, 1.8, "b", cex=1.5)
text(2, 1.8, "b", cex=1.5)
text(3, 1.8, "b", cex=1.5)
text(4, 1.8, "b", cex=1.5)
text(5, 1.8, "a", cex=1.5)
boxplot(Luteolin.7.O.Glc~Species,data=gr, col = c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        ylab = expression(paste("Luteolin-7-O-Glucoside  (", mu, "mol g/DW)")), xlab =  "", las=2,
        ylim=c(0,25),
        cex.lab=1.5, cex.axis=1.5)#***
text(1, 24.8, "bc", cex=1.5)
text(2, 24.8, "bc", cex=1.5)
text(3, 24.8, "bc", cex=1.5)
text(4, 24.8, "b", cex=1.5)
text(5, 24.8, "a", cex=1.5)
#s poly diff all, l. turion diff l. minu
boxplot(Apigenin.7.O.Glc~Species,data=gr, col = c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        ylab = expression(paste("Apigenin-7-O-Glucoside  (", mu, "mol g/DW)")), xlab =  "", las=2,
        ylim=c(0,5),
        cex.lab=1.5, cex.axis=1.5)#***
text(1, 4.6, "d", cex=1.5)
text(2, 4.6, "d", cex=1.5)
text(3, 4.6, "cd", cex=1.5)
text(4, 4.6, "b", cex=1.5)
text(5, 4.6, "a", cex=1.5)
#l turion diff l jp l mino, l minu diff all, s poly diff all
pairwise.wilcox.test(gr$Luteolin.8.C.Glc, gr$Species,
                     p.adjust.method = "BH")
boxplot(Luteolin.8.C.Glc~Species,data=gr, col = c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        ylab = expression(paste("Luteolin-8-C-Glucoside  (", mu, "mol g/DW)")), xlab =  "", xaxt = "n", las=2,
        ylim=c(0,18),
        cex.lab=1.5, cex.axis=1.5)#***
text(1, 17.8, "a", cex=1.5)
text(2, 17.8, "a", cex=1.5)
text(3, 17.8, "a", cex=1.5)
text(4, 17.8, "b", cex=1.5)
text(5, 17.8, "a", cex=1.5)
#l minu diff l. minor l. jap and l. turion, s poly
boxplot(Apigenin.8.C.Glc~Species,data=gr, col = c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        ylab = expression(paste("Apigenin-8-C-Glucoside  (", mu, "mol g/DW)")), xlab =  "", xaxt = "n", las=2,
        ylim=c(0,20),
        cex.lab=1.5, cex.axis=1.5)#***
text(1, 19.6, "d", cex=1.5)
text(2, 19.6, "c", cex=1.5)
text(3, 19.6, "ab", cex=1.5)
text(4, 19.6, "b", cex=1.5)
text(5, 19.6, "bd", cex=1.5)
#s poly diff all, l. turion diff to ljp lmino, l minu diff to l jp l turi, l mino diff l jp
boxplot(Chlorogenic.acid.1~Species,data=gr, col = c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        ylab = expression(paste("Chlorogenic acid  (", mu, "mol g/DW)")), xlab =  "", xaxt = "n", las=2,
        ylim=c(0,1),
        cex.lab=1.5, cex.axis=1.5)#***
text(1, 0.95, "a", cex=1.5)
text(2, 0.95, "a", cex=1.5)
text(3, 0.95, "a", cex=1.5)
text(4, 0.95, "a", cex=1.5)
text(5, 0.95, "a", cex=1.5)
dev.off()

#cant use hl_llmeans for sugar as missing data for sucrose
library(plyr)
sug_mu <- ddply(gr, "Accession", summarise, grp.mean=mean(Glucose))
sug_mu1 <- ddply(gr, "Accession", summarise, grp.mean=mean(Fructose))
sug_mu2 <- ddply(gr, "Accession", summarise, grp.mean=mean(Sucrose, na.rm=TRUE))
sug_mu3 <- ddply(gr, "Accession", summarise, grp.mean=mean(Starch.mg...g))
sum_sug <- cbind(sug_mu,sug_mu1,sug_mu2,sug_mu3)
sum_sug <- sum_sug[c(1,2,4,6,8)]
names(sum_sug) <- c("Accession", "Glucose", "Fructose", "Sucrose", "Starch")
names(sum_sug)

#sugars
pairwise.wilcox.test(gr$Glucose, gr$Species,
                     p.adjust.method = "BH")
pairwise.wilcox.test(gr$Fructose, gr$Species,
                     p.adjust.method = "BH")
pairwise.wilcox.test(gr$Sucrose, gr$Species,
                     p.adjust.method = "BH")
pairwise.wilcox.test(gr$Starch.mg...g, gr$Species,
                     p.adjust.method = "BH")

pairwise.wilcox.test(gr$Glucose, gr$Accession,
                     p.adjust.method = "BH")
pairwise.wilcox.test(gr$Fructose, gr$Accession,
                     p.adjust.method = "BH")
pairwise.wilcox.test(gr$Sucrose, gr$Accession,
                     p.adjust.method = "BH")
pairwise.wilcox.test(gr$Starch.mg...g, gr$Accession,
                     p.adjust.method = "BH")

#do a summary table of sugar content per species, with sd
#can conclude if high variability or not, accession level variation
#sugar
#cant use hl_llmeans for sugar as missing data for sucrose

library(dplyr)
#bulk dirty way to do it
att <- gr %>% group_by(Species) %>% summarise_each(funs(mean, sd), na.rm=T)
sug <- att[c(1,37:41,86:90)]
#sucrose not summarising
write.csv(sug, "Sugars_species_summ.csv")
#edited to add sucrose data in

#old?
library(plyr)
sug_mu <- ddply(gr, "Species", summarise, grp.mean=mean(Glucose), sd   = sd(Glucose, na.rm=TRUE),)
sug_mu1 <- ddply(gr, "Accession", summarise, grp.mean=mean(Fructose))
sug_mu2 <- ddply(gr, "Accession", summarise, grp.mean=mean(Sucrose, na.rm=TRUE))
sug_mu3 <- ddply(gr, "Accession", summarise, grp.mean=mean(Starch.mg...g))
sum_sug <- cbind(sug_mu,sug_mu1,sug_mu2,sug_mu3)

#try with just summs
#doesnt work without replicates
aov <- aov(Glucose~Accession,data=gr)
aov <- aov(Fructose~Accession,data=gr)
aov <- aov(Sucrose~Accession,data=gr)
aov <- aov(Starch.mg...g~Accession,data=gr)
summary(aov)
tuk_out <- TukeyHSD(aov, "Accession", conf.level=.95)
str(tuk_out)
pairwise.wilcox.test(gr$Glucose, gr$Accession,
                     p.adjust.method = "BH")
