#script to explore duckweed flavour data
#new data with all herbs included
#to make pca

setwd("C:\\Users\\kelli\\OneDrive - The University of Nottingham\\Nottingham\\Year 2\\UK Lemna - Bradford - Hull\\Flavour analysis\\Kellie - flavour\\Duckweed flavour analysis")

gr <- read.csv("Flavour_NEW2_mixedRT.csv")
#new with L. japonica species, subset total sp number + 1
gr <- read.csv("Flavour_NEW2_mixedRT+Ljp.csv")

#use these to filter by species type?
unique(ions_species)
#remove coriander and basil as nothing like duckweed
gr <- gr[-(109:116),]
#remove dandelion
gr <- gr[-(109:112),]
#remove spinach
gr <- gr[-(105:108),]
#rerun above gr into ions_all
ions_all <- gr

#pca with raw
ions_all <- gr[6:97]

#reduce species by number required
#ions_all <- grppb2[6:97]

#reduce species by number required with reduced compounds
#ions_all <- grppb2[6:65]

#define factors when using all species
ions_access <- gr[,2]
ions_species <- gr[,5]

#or
#define factors
#ions_access <- grppb2[,2]
#ions_species <- grppb2[,5]
str(ions_all)

print(ions_species) #prints all species every observation
levels(ions_species)

library(FactoMineR)
ions_acc <- ions_all[ ,c(1:92)] # selecting columns from csv
ions_acc.pca <- PCA(ions_acc, quali.sup=91) 
#cant exceed max col number

#for reduced no of compounds for 4 dws no basil comps
#or for raw? no changes
#ions_acc <- ions_all[ ,c(1:60)] # selecting columns from csv
#ions_acc.pca <- PCA(ions_acc, quali.sup=59) 
print(ions_acc.pca)
head(ions_acc.pca)
print(summary(ions_acc.pca)) #shows pc contirbutions of eigens

#FOR VARIABLE COS 2 GRAPH
#plot cos 2 as bar graph #high = good representation on pc
library(FactoMineR)
library(factoextra)

#basic plot
biplot(ions_acc.pca) #not working

#RUN THESE AFTER FILTERING BY SPECIES AND RE-RUN PCA
#8 species # now 9
fviz_pca_biplot(ions_acc.pca, repel=TRUE, pointsize=6, pointshape=21, col.var="black", arrowsize=0.6, labelsize=5, col.ind=ions_species, palette=c("green2", "gold", "blue", "purple", "red", "green", "pink", "orange", "grey"), addEllipses=TRUE, ellipse.type="confidence")
#now 7 species
fviz_pca_biplot(ions_acc.pca, repel=TRUE, pointsize=6, pointshape=21, col.var="black", arrowsize=0.6, labelsize=5, col.ind=ions_species, palette=c("green2", "gold", "blue", "purple", "red", "orange", "pink"), addEllipses=TRUE, ellipse.type="confidence")
#6 species
fviz_pca_biplot(ions_acc.pca, repel=TRUE, pointsize=6, pointshape=21, col.var="black", arrowsize=0.6, labelsize=5, col.ind=ions_species, palette=c("green2", "gold", "blue", "purple", "red", "orange"), addEllipses=TRUE, ellipse.type="confidence")
#5 species
fviz_pca_biplot(ions_acc.pca, repel=TRUE, pointsize=6, pointshape=21, col.var="black", arrowsize=0.6, labelsize=5, col.ind=ions_species, palette=c("green2", "blue", "purple", "red", "orange"), addEllipses=TRUE, ellipse.type="confidence")
#just dws
fviz_pca_biplot(ions_acc.pca, repel=TRUE, pointsize=6, pointshape=21, col.var="black", arrowsize=0.6, labelsize=5, col.ind=ions_species, palette=c("green2", "blue", "purple", "red"), addEllipses=TRUE, ellipse.type="confidence")
#23 and 40 outliers

fviz_cos2(ions_acc.pca, choice = "var", axes = 1:2) #K and S biggest
#top 5 penten-3-ol, penten-1-ol z, beta-cyclitrol,trans-beta-ionome,cyclohexanol 2,6 dimethyl

fviz_pca_var(ions_acc.pca, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
)

#do this first to get pc columns, then visualise variable
ions_acc.pca<-prcomp(ions_acc[ ,1:91],center=T,scale=T) #list can include mix mat and df
#or for reduced compounds no basil
#ions_acc.pca<-prcomp(ions_acc[ ,1:60],center=T,scale=T)
str(ions_acc.pca)
mypc <- ions_acc.pca$x #define new variable x = pcs just need to plot these on xy graph

#plot accession, cant tell which which accession n =25 but all
#reps so not colored same for each rep
plot(mypc[,1], mypc[,2], col = my_pal,
     las=1, xlab="PC1 (26%)", ylab="PC2 (15%)",
     pch=16, cex=1.5, xlim=c(-10,15), ylim=c(-10,15)) #pcs as columns, produce xy plot, las is rotation of axis numbers, pch plot shape, ylim expand out so legend room
abline(v=0, lty=2, col="lightgrey") #draw line, lty is segmented
abline(h=0, lty=2, col="lightgrey") #0 lines dashed
legend("right", pch=16, col=my_pal, cex=1, c("KS02", "KS03", "KS04", "KS06A", "KS06B",
                                             "KS09", "KS12", "KS13", "KS14", "KS15",
                                             "KS16", "KS17","KS18", "KS20", "KS21", "KS22", "KS25",
                                             "KS28", "KS29", "KS66A", "KS77A", "KS78A",
                                             "LY01A", "LY01B", "Nuff1", "Spinach", "Coriander",
                                             "Basil", "Dandelion")) #customise legend seperately cex=txtsize cols 1-4 stnd, concat order of places as want to display
#too busy #can show outliers
text(x=mypc[,1],y=mypc[,2], labels =ions_access, pos=2) #txt and pos to define labels on plot, 2 = top of plot, 1 bottom

#need more than 12 colors as repeating
my_pal <- scico::scico(length(unique(ions_species)), palette = "batlow")
my_pal <- scico::scico(length(unique(ions_access)), palette = "lisbon")


par(mfrow = c(1, 2))
#species all raw data
plot(mypc[,1], mypc[,2], col = ions_access,
     las=1, xlab="PC1 (26%)", ylab="PC2 (15%)",
     pch=16, cex=1.5, xlim=c(-40,40), ylim=c(-40,40)) #pcs as columns, produce xy plot, las is rotation of axis numbers, pch plot shape, ylim expand out so legend room
abline(v=0, lty=2, col="lightgrey") #draw line, lty is segmented
abline(h=0, lty=2, col="lightgrey") #0 lines dashed
legend("right", pch=16, col=unique(ions_access), cex=1, c("L. minor", "L. minuta", "S. polyrhiza", "L. turionifera", "Spinach", "Basil", "Coriander", "Dandelion")) #customise legend seperately cex=txtsize cols 1-4 stnd, concat order of places as want to display
#to show if coloring/grouping properly
text(x=mypc[,1],y=mypc[,2], labels =ions_access, pos=2) #txt and pos to define labels on plot, 2 = top of plot, 1 bottom
#dev.new()

#species 6 left
#species all raw data
plot(mypc[,1], mypc[,2], col = ions_species,
     las=1, xlab="PC1 (26%)", ylab="PC2 (15%)",
     pch=16, cex=1.5, xlim=c(-20,20), ylim=c(-20,20)) #pcs as columns, produce xy plot, las is rotation of axis numbers, pch plot shape, ylim expand out so legend room
abline(v=0, lty=2, col="lightgrey") #draw line, lty is segmented
abline(h=0, lty=2, col="lightgrey") #0 lines dashed
legend("right", pch=16, col=unique(ions_species), cex=1, c("L. minor", "L. minuta", "S. polyrhiza", "L. turionifera", "Spinach", "Dandelion")) #customise legend seperately cex=txtsize cols 1-4 stnd, concat order of places as want to display
#to show if coloring/grouping properly
text(x=mypc[,1],y=mypc[,2], labels =ions_access, pos=2) #txt and pos to define labels on plot, 2 = top of plot, 1 bottom
#dev.new()

levels(ions_species)
#9 species
ions_species <- factor(ions_species, levels=c("L. minor", "L. japonica", "L. turionifera", "L. minuta", 
                                              "S. polyrhiza", "Basil", "Coriander", "Spinach",
                                              "Dandelion"))
#8 species
#ions_species <- factor(ions_species, levels=c("L. minor", "L. minuta", "L. turionifera",
#                                "S. polyrhiza", "Basil", "Coriander", "Spinach",
#                                "Dandelion"))

#7 species
ions_species <- factor(ions_species, levels=c("L. minor", "L. japonica", "L. turionifera", "L. minuta",
                                              "S. polyrhiza", "Spinach",
                                              "Dandelion"))
#6 species
#ions_species <- factor(ions_species, levels=c("L. minor", "L. minuta", "L. turionifera",
#                                              "S. polyrhiza", "Spinach","Dandelion"))

#just dw species now 5
ions_species <- factor(ions_species, levels=c("L. minor", "L. japonica", "L. turionifera", "L. minuta",
                                              "S. polyrhiza"))

#just dw species
ions_species <- factor(ions_species, levels=c("L. minor", "L. minuta", "L. turionifera",
                                              "S. polyrhiza"))

#8species col scheme
#not been used in ggplots
my_cols <- c("red", "blue", "green", "purple", "gold", "black", "orange", "pink")

#try ggplot for pca
#install.packages("ggforce") #used to plot all elipses <3 points
library(ggforce)
PC1<-ions_acc.pca$x[,1]
PC2<-ions_acc.pca$x[,2]
ggplot(ions_all, 
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

#try change cols
pca <- ggplot(ions_all, aes(x = PC1, y = PC2, color = ions_species)) +
  geom_point(size=5) +
  stat_ellipse() +
  theme_classic() +
  scale_color_manual(values = c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00", "lightblue", "gold", "#FF0535", "grey")) +
  geom_hline(yintercept=0, color="gray", linetype="dashed") + 
  geom_vline(xintercept=0, color="gray", linetype="dashed") +
  scale_x_continuous(limits = c(-50, 50), name="PC1 (31%)") +
  scale_y_continuous(limits = c(-50, 50), name="PC2 (22%)") +
  theme(legend.title=element_blank()) +
  theme(legend.text = element_text(face = "italic"))
ggsave("pca_flavour+herbs_col_newcol+Ljp.tiff", dpi = 300, width = 17, height = 15 , units = "cm") 

pca

#7 sp when removed basil and coriander
pcax <- ggplot(ions_all, aes(x = PC1, y = PC2, color = ions_species)) +
  geom_point(size=5) +
 stat_ellipse() +
  theme_classic() +
  scale_color_manual(values = c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00", "#FF0535", "grey")) +
  geom_hline(yintercept=0, color="gray", linetype="dashed") + 
  geom_vline(xintercept=0, color="gray", linetype="dashed") +
  scale_x_continuous(limits = c(-50, 50), name="PC1 (31%)") +
  scale_y_continuous(limits = c(-50, 50), name="PC2 (22%)") +
  theme(legend.title=element_blank()) +
  theme(legend.text = element_text(face = "italic"))
ggsave("pca_flavour+herbs_col+Ljp.tiff", dpi = 300, width = 17, height = 15 , units = "cm") 
pcax

#6 species
pca1 <- ggplot(ions_all, aes(x = PC1, y = PC2, color = ions_species)) +
  geom_point(size=5) +
  stat_ellipse() +
  theme_classic() +
  scale_color_manual(values = c("red", "blue", "green", "purple", "orange", "pink")) +
  geom_hline(yintercept=0, color="gray", linetype="dashed") + 
  geom_vline(xintercept=0, color="gray", linetype="dashed") +
  scale_x_continuous(limits = c(-45, 45), name="PC1 (31%)") +
  scale_y_continuous(limits = c(-45, 45), name="PC2 (22%)") +
  theme(legend.title=element_blank()) +
  theme(legend.text = element_text(face = "italic"))
ggsave("pca_flavour+da+sp_col_newcol.tiff", dpi = 300, width = 17, height = 15 , units = "cm") 

pca1

#dw species #now 5
pca2 <- ggplot(ions_all, aes(x = PC1, y = PC2, color = ions_species)) +
  geom_point(size=5) +
  stat_ellipse() +
  theme_classic() +
  scale_color_manual(values = c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00")) +
  geom_hline(yintercept=0, color="gray", linetype="dashed") + 
  geom_vline(xintercept=0, color="gray", linetype="dashed") +
  scale_x_continuous(limits = c(-20, 20), name="PC1 (31%)") +
  scale_y_continuous(limits = c(-20, 20), name="PC2 (22%)") +
  theme(legend.title=element_blank()) +
  theme(legend.text = element_text(face = "italic"))
ggsave("pca_flavour+dw_only_col_newcol+Ljp.tiff", dpi = 300, width = 17, height = 15 , units = "cm") 

pca2

#dw species
pca3 <- ggplot(ions_all, aes(x = PC1, y = PC2, color = ions_species)) +
  geom_point(size=5) +
  stat_ellipse() +
  theme_classic() +
  scale_color_manual(values = c("red", "blue", "green", "purple")) +
  geom_hline(yintercept=0, color="gray", linetype="dashed") + 
  geom_vline(xintercept=0, color="gray", linetype="dashed") +
  scale_x_continuous(limits = c(-20, 20), name="PC1 (31%)") +
  scale_y_continuous(limits = c(-20, 20), name="PC2 (22%)") +
  theme(legend.title=element_blank()) +
  theme(legend.text = element_text(face = "italic"))
ggsave("pca_flavour+dw_only_redcomps.tiff", dpi = 300, width = 17, height = 15 , units = "cm") 

pca3
library(gridExtra)
#cant run them together as use same scripts, only last 1 saved
grid.arrange(pca, pcax, pca2, ncol=2)
#pcax?
  
#functional groups
#NOT DONE YET?
#species no basil or coriander
plot(mypc[,1], mypc[,2], col = ions_species,
     las=1, xlab="PC1 (26%)", ylab="PC2 (15%)",
     pch=16, cex=1.5, xlim=c(-20,40), ylim=c(-20,20)) #pcs as columns, produce xy plot, las is rotation of axis numbers, pch plot shape, ylim expand out so legend room
abline(v=0, lty=2, col="lightgrey") #draw line, lty is segmented
abline(h=0, lty=2, col="lightgrey") #0 lines dashed
legend("right", pch=16, col=unique(ions_species), cex=1, c("L. minor", "L. minuta", "S. polyrhiza", "L. turionifera", "Spinach", "Dandelion")) #customise legend seperately cex=txtsize cols 1-4 stnd, concat order of places as want to display
#to show if coloring/grouping properly
text(x=mypc[,1],y=mypc[,2], labels =ions_access, pos=2) #txt and pos to define labels on plot, 2 = top of plot, 1 bottom
#dev.new()

#species no basil or coriander or dandelion
plot(mypc[,1], mypc[,2], col = ions_species,
     las=1, xlab="PC1 (26%)", ylab="PC2 (15%)",
     pch=16, cex=1.5, xlim=c(-20,20), ylim=c(-20,20)) #pcs as columns, produce xy plot, las is rotation of axis numbers, pch plot shape, ylim expand out so legend room
abline(v=0, lty=2, col="lightgrey") #draw line, lty is segmented
abline(h=0, lty=2, col="lightgrey") #0 lines dashed
legend("right", pch=16, col=unique(ions_species), cex=1, c("L. minor", "L. minuta", "S. polyrhiza", "L. turionifera", "Spinach")) #customise legend seperately cex=txtsize cols 1-4 stnd, concat order of places as want to display
#to show if coloring/grouping properly
text(x=mypc[,1],y=mypc[,2], labels =ions_access, pos=2) #txt and pos to define labels on plot, 2 = top of plot, 1 bottom
#dev.new()

#run pca with summary data for ppb only so only 25 points
#work out ppb data
#gr <- read.csv("Flavour_NEW2_mixedRT.csv")
setwd("C:\\Users\\kelli\\OneDrive - The University of Nottingham\\Nottingham\\Year 2\\UK Lemna - Bradford - Hull\\Flavour analysis\\Kellie - flavour\\Duckweed flavour analysis")
#read in ppb data
grppb <- read.csv("Flavour_NEW2_mixedRT_PPB.csv")
#add ljp
grppb <- read.csv("Flavour_NEW2_mixedRT_PPB+Ljp.csv")
library(dplyr)

#do group means to see which not in duckweeds
#one by one
#test
library(plyr)
mu <- ddply(grppb, "Species", summarise, grp.mean=mean(Acetaldehyde), grp.sd=sd(Acetaldehyde))

#loop this function so do for each column?
library(dplyr)
sum <- grppb %>% group_by(Species) %>% 
  summarise_all(.funs = c(mean="mean"))

#write.csv(sum, "Flav_speciesmeans_tocutoutnondwcompounds+Ljp.csv")
#write.csv(sum, "Flav_speciesmeans_tocutoutnondwcompounds.csv")

#these groupings tell which <1ppb in duckweed, manually remove
#read back in with reduced compounds in dataset
setwd("C:\\Users\\kelli\\OneDrive - The University of Nottingham\\Nottingham\\Year 2\\UK Lemna - Bradford - Hull\\Flavour analysis\\Kellie - flavour\\Duckweed flavour analysis")
#read in ppb data
grppb <- read.csv("Flavour_NEW2_mixedRT_PPB_low_in_dw.csv")
grppb <- read.csv("Flavour_NEW2_mixedRT_PPB_low_in_dw+Ljp.csv")

#colnames(gr)[which(names(gr) == "Ethyl.butyl.ketone")] <- "X3.Heptanone"
#names(gr)
#Summary <- gr %>%
#  summarise(X.3.Heptanone_avg = mean(X3.Heptanone), stdev = sd(X3.Heptanone), N_n= n())
#Summary
#
#Summary$X.3.Heptanone_avg

#formula
#ppb = (x / 45667985 * 22)) where x is column number
#do in excel, cant do it
setwd("C:\\Users\\kelli\\OneDrive - The University of Nottingham\\Nottingham\\Year 2\\UK Lemna - Bradford - Hull\\Flavour analysis\\Kellie - flavour\\Duckweed flavour analysis\\test")
setwd("C:\\Users\\kelli\\OneDrive - The University of Nottingham\\Nottingham\\Year 2\\UK Lemna - Bradford - Hull\\Flavour analysis\\Kellie - flavour\\Duckweed flavour analysis\\test2")
#dws only
setwd("C:\\Users\\kelli\\OneDrive - The University of Nottingham\\Nottingham\\Year 2\\UK Lemna - Bradford - Hull\\Flavour analysis\\Kellie - flavour\\Duckweed flavour analysis\\test3")
#dw + Ljp
setwd("C:\\Users\\kelli\\OneDrive - The University of Nottingham\\Nottingham\\Year 2\\UK Lemna - Bradford - Hull\\Flavour analysis\\Kellie - flavour\\Duckweed flavour analysis\\test4Ljp")

grp <- grppb[,(6:65)] #just numeric
grp <- grppb2[,(6:65)] #just numeric

#just show on screen, not good if many
for (i in 1:ncol(grp)) {
  boxplot(grp[, i] ~ grppb2$Species)
}

#print to file and flick through
for (i in 1:ncol(grp)) {
  png(file = paste(names(grp)[i], ".png", sep=""))
  boxplot(grp[, i] ~ grppb2$Species, ylab = names(grp)[i], 
          xlab = "", las=2)
  dev.off()
}

#dev.new()
#many compounds are basil specific need to reduce these manually
#.gamma..Muurolene,
#Bergamotene..alpha...cis.. etc
#remove coriander and basil as nothing like duckweed
library(dplyr)
grppb2 <- filter(grppb, Species != "Coriander") %>% droplevels()
grppb2 <- filter(grppb2, Species != "Basil") %>% droplevels()
grppb2 <- filter(grppb2, Species != "Dandelion") %>% droplevels()
grppb2 <- filter(grppb2, Species != "Spinach") %>% droplevels()

#104 obs just 4 dw species
species <- grppb2[,5]
unique(species)
#by indexing leaves species names in graphs
#grppb <- grppb[-(109:116),]
#re run above

#loop kruskal wallis for columns
results <- list()
for(i in names(grppb2[,6:97])){  
  results[[i]] <- kruskal.test(formula(paste(i, "~ Species")), data = grppb2)
}

results <- list()
for(i in names(grppb2[,6:65])){  
  results[[i]] <- kruskal.test(formula(paste(i, "~ Species")), data = grppb2)
}
r1 <- unlist(results)
r1 <- as.data.frame(r1)
write.csv(r1, "KWtest_Ljp.csv")

#manually check file for sig compounds
  
#manual check
kruskal.test(Acetaldehyde ~ Species, data = grppb2)
kruskal.test(X2.Thiapropane ~ Species, data = grppb2)
kruskal.test(Propanal..2.methyl. ~ Species, data = grppb2)
#matches results

#takes sig compounds from df as a list
kw <- read.csv("KWtest_sigpvals.csv")
compounds <- print(kw$Compound)
sapply(compounds, levels)
compounds <- as.data.frame(compounds)
print(compounds)
#remove .p.value from names
compounds$compounds <-gsub("p.value","",as.character(compounds$compounds))

#manual look for sig compounds 6 species = 64
#sig compounds just dw: 4 species = 26

#find in grppb2 and do boxplots and further stats

#26 compounds sig diff by species see KW test results
#still same number with Ljp, any diff comps?

#look at ones most prevelant in dw
grppb2[do.call(order, grppb2),]
#min max of each column
apply(grppb2,2,min)
apply(grppb2,2,max) #most useful
#dont know how to arrange cols by max

max(grppb2$Acetaldehyde)

#make raw data with just significant 26 compounds in
setwd("C:\\Users\\kelli\\OneDrive - The University of Nottingham\\Nottingham\\Year 2\\UK Lemna - Bradford - Hull\\Flavour analysis\\Kellie - flavour\\Duckweed flavour analysis")
sig <- read.csv("Flavour_NEW2_mixedRT_PPB_species_justsig.csv")
sig <- read.csv("Flavour_NEW2_mixedRT_PPB_species_justsig+Ljp.csv")

#remove herbs
library(dplyr)
sig <- filter(sig, Species != "Coriander") %>% droplevels()
sig <- filter(sig, Species != "Basil") %>% droplevels()
sig <- filter(sig, Species != "Dandelion") %>% droplevels()
sig <- filter(sig, Species != "Spinach") %>% droplevels()

#summarise by species and SD
#avgs per species
library(dplyr)
colMeans(grppb2[6:97])

dw <- read.csv("Flavour_NEW2_mixedRT_PPB_low_in_dw+Ljp_rename.csv")
result1 <- dw %>% 
  group_by(Species) %>%
  summarise_all("mean")

write.csv(result1, "Flavour_dwcompounds_speciesavgs_Ljp.csv")

acc <- dw %>% 
  group_by(Accession) %>%
  summarise_all("mean")

write.csv(acc, "Flavour_dwcompounds_accessionavgs_Ljp.csv")

#this was used to add to combined script for accession correlations

result <- sig %>% 
  group_by(Species) %>%
  summarise_all("mean")

write.csv(result, "Flavour_31compounds_speciesavgs_Ljp.csv")

result_sig <- grppb2 %>% 
  group_by(Species) %>%
  summarise_all("mean")

write.csv(result_sig, "Flavour_64compounds_speciesavgs_Ljp.csv")

#remove herbs
library(dplyr)
result_sig <- filter(result_sig, Species != "Coriander") %>% droplevels()
result_sig <- filter(result_sig, Species != "Basil") %>% droplevels()
result_sig <- filter(result_sig, Species != "Dandelion") %>% droplevels()
result_sig <- filter(result_sig, Species != "Spinach") %>% droplevels()


#these not working
#sum <- grppb2 %>% group_by(Accession) %>% select(-Component.Name, -Rep, -Sample.no, -Species) %>% summarise_all(mean)
#sum
#sds per species
#sum3 <- grppb2 %>% group_by(Species) %>% select(-Component.Name, -Rep, -Sample.no, -Accession) %>% summarise_all(sd)
#sum3

#plot in new heatmap
library(RColorBrewer)
#scale numeric data
ions_all <- grppb2[6:97]
#reduced no as removed basil compounds
ions_all <- grppb2[6:65]

ions_all <- result[6:97]
ions_all <- result_sig[6:31]

ions_all <- result[6:31]

ions_sc <- scale(ions_all,center=T,scale=T)
ions_sc <- as.matrix(ions_sc)

#assign row names from df
rownames(ions_sc) <- grppb2$Accession
rownames(ions_sc) <- grppb2$Species

#current for Ljp sig compounds
rownames(ions_sc) <- result_sig$Species

heatmap(ions_sc)
heatmap(ions_sc, Colv = NA, Rowv = NA, scale="column")
#change color scheme yellow blue
colfunc <- colorRampPalette(c("blue", "yellow"))
heatmap(ions_sc,col=colfunc(11),scale="row")

#tiff  normally shows full dataset but not here
#tto make pretty, rename, add sclae bar?
par(mar=c(7,4,4,2)+0.1)
#tiff('Heatmap_species_flavours_1.tiff', units="in", width=14, height=12, res=300, compression = 'lzw')
#tiff('Heatmap_allflav_avgs_species_newsp.tiff', units="in", width=12, height=12, res=300, compression = 'lzw')
tiff('Heatmap_flavsig_avgs_species_newsp_Ljp.tiff', units="in", width=15, height=12, res=300, compression = 'lzw')
par(mfrow = c(3, 3),  mar=c(5,4.5,4,2))
colfunc <- colorRampPalette(c("blue", "yellow"))
heatmap(ions_sc,col=colfunc(15),scale="column", cexCol=1.2,margins=c(12,8))
#heatmap(ions_sc,col=colfunc(15),scale="column",cexCol=0.9,margins=c(12,8))
dev.off()

#boxplots OLD
#compounds higher in l. minuta
#furan-3-methyl, 1-hexanol, furan-alpha-ethyl, hexen-1-ol, thiapropane,
#pentanedione, linalyl anthranilate

#NEW
#linalyl anthranilate

#try statistical test
pairwise.wilcox.test(sig$Furan..3.methyl., sig$Species,
                     p.adjust.method = "BH")
#n.s
pairwise.wilcox.test(sig$X1.Hexanol, sig$Species,
                     p.adjust.method = "BH")
# l minu sig to l turion
pairwise.wilcox.test(sig$Furan..alpha...ethyl.., sig$Species,
                     p.adjust.method = "BH")
#l minor diff to all, l minu diff to l turion, s poly diff to l turion
pairwise.wilcox.test(sig$X3.Hexen.1.ol...Z.., sig$Species,
                     p.adjust.method = "BH")
#l turion diff to all
pairwise.wilcox.test(sig$X2.Thiapropane, sig$Species,
                     p.adjust.method = "BH")
#l mino diff to l minu
pairwise.wilcox.test(sig$X2.3.Pentanedione, sig$Species,
                     p.adjust.method = "BH")
#s poly diff to all
pairwise.wilcox.test(sig$Linalyl.anthranilate, sig$Species,
                     p.adjust.method = "BH")
#n.s

#cut out Furan..3.methyl, linalyl.anthranilate

#plots 1
par(mfrow = c(4, 6))
boxplot(X1.Hexanol~Species,data=sig, col = c("red", "blue", "darkgreen", "purple"),
        ylab = "1-Hexanol", xlab =  "", las=2,
        ylim=c(0,300),
        cex.lab=1.5, cex.axis=1.5)#***
text(1, 280, "ab", cex=1.5)
text(2, 280, "a", cex=1.5)
text(3, 280, "b", cex=1.5)
text(4, 280, "ab", cex=1.5)
boxplot(Furan..alpha...ethyl..~Species,data=sig, col = c("red", "blue", "darkgreen", "purple"),
        ylab = "2-Ethylfuran", xlab =  "", las=2,
        ylim=c(0,15),
        cex.lab=1.5, cex.axis=1.5)#***
text(1, 14.5, "b", cex=1.5)
text(2, 14.5, "a", cex=1.5)
text(3, 14.5, "c", cex=1.5)
text(4, 14.5, "ab", cex=1.5)
boxplot(X3.Hexen.1.ol...Z..~Species,data=sig, col = c("red", "blue", "darkgreen", "purple"),
        ylab = "cis-3-Hexen-1-ol", xlab =  "", las=2,
        ylim=c(0,220),
        cex.lab=1.5, cex.axis=1.5)#***
text(1, 219, "a", cex=1.5)
text(2, 219, "a", cex=1.5)
text(3, 219, "b", cex=1.5)
text(4, 219, "a", cex=1.5)
boxplot(X2.Thiapropane~Species,data=sig, col = c("red", "blue", "darkgreen", "purple"),
        ylab = "Dimethyl sulfide", xlab =  "", las=2,
        ylim=c(0,80),
        cex.lab=1.5, cex.axis=1.5)#***
text(1, 79, "b", cex=1.5)
text(2, 79, "a", cex=1.5)
text(3, 79, "ab", cex=1.5)
text(4, 79, "ab", cex=1.5)
boxplot(X2.3.Pentanedione~Species,data=sig, col = c("red", "blue", "darkgreen", "purple"),
        ylab = "2,3-Pentanedione", xlab =  "", las=2,
        ylim=c(0,40),
        cex.lab=1.5, cex.axis=1.5)#***
text(1, 39.5, "a", cex=1.5)
text(2, 39.5, "a", cex=1.5)
text(3, 39.5, "a", cex=1.5)
text(4, 39.5, "b", cex=1.5)

#look higher in S.poly
#Hexen-1-ol (E), 2-hexenal, 5-Isopropyl-2-methylbicyclo[3.1.0]hexan-2-ol,
#Tetradec-2-enal <trans->, cis,cis-7,10,-Hexadecadienal
pairwise.wilcox.test(sig$X2.Hexen.1.ol...E.., sig$Species,
                     p.adjust.method = "BH")
#l mino diff to l minu, s poly diff to l mino and l turion
pairwise.wilcox.test(sig$X2.Hexenal, sig$Species,
                     p.adjust.method = "BH")
#l turion diff to all
#s poly diff to l mino and l turionn
pairwise.wilcox.test(sig$X5.Isopropyl.2.methylbicyclo.3.1.0.hexan.2.ol..., sig$Species,
                     p.adjust.method = "BH")
#s poly diff to all, l mino and l minu diff
pairwise.wilcox.test(sig$Tetradec.2.enal..trans.., sig$Species,
                     p.adjust.method = "BH")
#l mino diff to l minu, l minu diff to l turion and s poly
pairwise.wilcox.test(sig$cis.cis.7.10..Hexadecadienal, sig$Species,
                     p.adjust.method = "BH")
# l minu diff to l mino, l minu diff to l turion and s poly

#plots 2
boxplot(X2.Hexen.1.ol...E..~Species,data=sig, col = c("red", "blue", "darkgreen", "purple"),
        ylab = "trans-2-Hexen-1-ol", xlab =  "", las=2,
        ylim=c(0,80),
        cex.lab=1.5, cex.axis=1.5)#***
text(1, 79, "b", cex=1.5)
text(2, 79, "a", cex=1.5)
text(3, 79, "b", cex=1.5)
text(4, 79, "a", cex=1.5)
boxplot(X2.Hexenal~Species,data=sig, col = c("red", "blue", "darkgreen", "purple"),
        ylab = "trans-2-Hexenal", xlab =  "", las=2,
        ylim=c(0,500),
        cex.lab=1.5, cex.axis=1.5)#***
text(1, 490, "b", cex=1.5)
text(2, 490, "ab", cex=1.5)
text(3, 490, "c", cex=1.5)
text(4, 490, "a", cex=1.5)
boxplot(X5.Isopropyl.2.methylbicyclo.3.1.0.hexan.2.ol...~Species,data=sig, col = c("red", "blue", "darkgreen", "purple"),
        ylab = "Sabinene hydrate", xlab =  "", las=2,
        ylim=c(0,3),
        cex.lab=1.5, cex.axis=1.5)#***
text(1, 2.9, "a", cex=1.5)
text(2, 2.9, "b", cex=1.5)
text(3, 2.9, "ab", cex=1.5)
text(4, 2.9, "a", cex=1.5)
boxplot(Tetradec.2.enal..trans..~Species,data=sig, col = c("red", "blue", "darkgreen", "purple"),
        ylab = "trans-Tetradec-2-enal", xlab =  "", las=2,
        ylim=c(0,10),
        cex.lab=1.5, cex.axis=1.5)#***
text(1, 9.8, "a", cex=1.5)
text(2, 9.8, "b", cex=1.5)
text(3, 9.8, "ab", cex=1.5)
text(4, 9.8, "a", cex=1.5)
boxplot(cis.cis.7.10..Hexadecadienal~Species,data=sig, col = c("red", "blue", "darkgreen", "purple"),
        ylab = "cis,cis-7,10,-Hexadecadienal", xlab =  "", las=2,
        ylim=c(0,20),
        cex.lab=1.5, cex.axis=1.5)#***
text(1, 19.8, "ab", cex=1.5)
text(2, 19.8, "b", cex=1.5)
text(3, 19.8, "a", cex=1.5)
text(4, 19.8, "a", cex=1.5)

#look higher in l turion
#Pentene <3-hydroxy->
pairwise.wilcox.test(sig$Pentene..3.hydroxy.., sig$Species,
                     p.adjust.method = "BH")
#n.s CUT FROM ANALYSIS

#rem l minor
#X1.octen.3.ol, Pent-2-enal, Heptanal, 5-Ethylcyclopent-1-enecarboxaldehyde,
#Pent-2-zenol, Bicyclo[3.1.1]heptane, 6,6-dimethyl-2-methylene-, (1S)-,
#Hept-5-en-2-one <6-methyl->, Tetradec-(11Z)-enal,
#Tetradecanal, pentadecanal, pyrrole, tridecanal, Gerylacetone
pairwise.wilcox.test(sig$X1.Octen.3.ol, sig$Species,
                     p.adjust.method = "BH")
# l mino diff to l minu, s poly diff to l mino and l turion
pairwise.wilcox.test(sig$Pent..2E..enal, sig$Species,
                     p.adjust.method = "BH")
# l mino diff to all
pairwise.wilcox.test(sig$X5.Ethylcyclopent.1.enecarboxaldehyde, sig$Species,
                     p.adjust.method = "BH")
# l mino diff to ljp, l tu diff to l mino l jp, l minu diff to l mino l jp, s poly diff to l mino l jp l turi

pairwise.wilcox.test(sig$Pent..2Z..enol, sig$Species,
                     p.adjust.method = "BH")
#ns cut from analysis

pairwise.wilcox.test(sig$Heptanal..n.., sig$Species,
                     p.adjust.method = "BH")
#l mino diff to all
pairwise.wilcox.test(sig$Hept.5.en.2.one..6.methyl.., sig$Species,
                     p.adjust.method = "BH")
# l mino diff to l minu
pairwise.wilcox.test(sig$Tetradec..11Z..enal, sig$Species,
                     p.adjust.method = "BH")
#l mino diff to l minu, l minu diff to l turion
pairwise.wilcox.test(sig$Tetradecanal, sig$Species,
                     p.adjust.method = "BH")
#l mino diff to l minu, l minu diff to l turion and s poly
pairwise.wilcox.test(sig$Pentadecanal., sig$Species,
                     p.adjust.method = "BH")
#l mino diff to l minu, l minu diff to l turion
pairwise.wilcox.test(sig$Pyrrole, sig$Species,
                     p.adjust.method = "BH")
#l turion diff to all, s poly diff to all
pairwise.wilcox.test(sig$Tridecanal..n.., sig$Species,
                     p.adjust.method = "BH")
#l mino diff to l inu, l ,minu diff to l turion
pairwise.wilcox.test(sig$trans.Geranylacetone, sig$Species,
                     p.adjust.method = "BH")
#l mino diff to l minu

boxplot(X1.Octen.3.ol~Species,data=sig, col = c("red", "blue", "darkgreen", "purple"),
        ylab = "1-Octen-3-ol", xlab =  "", las=2,
        ylim=c(0,120),
        cex.lab=1.5, cex.axis=1.5)#***
text(1, 119.6, "a", cex=1.5)
text(2, 119.6, "b", cex=1.5)
text(3, 119.6, "ab", cex=1.5)
text(4, 119.6, "b", cex=1.5)
boxplot(Pent..2E..enal~Species,data=sig, col = c("red", "blue", "darkgreen", "purple"),
        ylab = "(E)-pent-2-enal", xlab =  "", las=2,
        ylim=c(0,200),
        cex.lab=1.5, cex.axis=1.5)#***
text(1, 199.6, "a", cex=1.5)
text(2, 199.6, "b", cex=1.5)
text(3, 199.6, "b", cex=1.5)
text(4, 199.6, "b", cex=1.5)
boxplot(X5.Ethylcyclopent.1.enecarboxaldehyde~Species,data=sig, col = c("red", "blue", "darkgreen", "purple"),
        ylab = "5-Ethylcyclopent-1-enecarboxaldehyde", xlab =  "", las=2,
        ylim=c(0,130),
        cex.lab=1.5, cex.axis=1.5)#***
text(1, 128.6, "a", cex=1.5)
text(2, 128.6, "b", cex=1.5)
text(3, 128.6, "b", cex=1.5)
text(4, 128.6, "b", cex=1.5)
boxplot(Heptanal..n..~Species,data=sig, col = c("red", "blue", "darkgreen", "purple"),
        ylab = "Heptanal", xlab =  "", las=2,
        ylim=c(0,40),
        cex.lab=1.5, cex.axis=1.5)#***
text(1, 38.6, "a", cex=1.5)
text(2, 38.6, "b", cex=1.5)
text(3, 38.6, "b", cex=1.5)
text(4, 38.6, "b", cex=1.5)
boxplot(Hept.5.en.2.one..6.methyl..~Species,data=sig, col = c("red", "blue", "darkgreen", "purple"),
        ylab = "Sulcatone", xlab =  "", las=2,
        ylim=c(0,18),
        cex.lab=1.5, cex.axis=1.5)#***
text(1, 17.6, "a", cex=1.5)
text(2, 17.6, "b", cex=1.5)
text(3, 17.6, "ab", cex=1.5)
text(4, 17.6, "ab", cex=1.5)
boxplot(Tetradec..11Z..enal~Species,data=sig, col = c("red", "blue", "darkgreen", "purple"),
        ylab = "Z-11-Tetradecenal", xlab =  "", las=2,
        ylim=c(0,25),
        cex.lab=1.5, cex.axis=1.5)#***
text(1, 24.6, "a", cex=1.5)
text(2, 24.6, "b", cex=1.5)
text(3, 24.6, "a", cex=1.5)
text(4, 24.6, "ab", cex=1.5)
boxplot(Tetradecanal~Species,data=sig, col = c("red", "blue", "darkgreen", "purple"),
        ylab = "Tetradecanal", xlab =  "", las=2,
        ylim=c(0,120),
        cex.lab=1.5, cex.axis=1.5)#***
text(1, 118.6, "a", cex=1.5)
text(2, 118.6, "b", cex=1.5)
text(3, 118.6, "a", cex=1.5)
text(4, 118.6, "a", cex=1.5)
boxplot(Pentadecanal.~Species,data=sig, col = c("red", "blue", "darkgreen", "purple"),
        ylab = "Pentadecanal", xlab =  "", las=2,
        ylim=c(0,280),
        cex.lab=1.5, cex.axis=1.5)#***
text(1, 278.6, "a", cex=1.5)
text(2, 278.6, "b", cex=1.5)
text(3, 278.6, "a", cex=1.5)
text(4, 278.6, "ab", cex=1.5)
boxplot(Pyrrole~Species,data=sig, col = c("red", "blue", "darkgreen", "purple"),
        ylab = "Pyrrole", xlab =  "", las=2,
        ylim=c(0,10),
        cex.lab=1.5, cex.axis=1.5)#***
text(1, 9.6, "b", cex=1.5)
text(2, 9.6, "b", cex=1.5)
text(3, 9.6, "a", cex=1.5)
text(4, 9.6, "c", cex=1.5)
boxplot(Tridecanal..n..~Species,data=sig, col = c("red", "blue", "darkgreen", "purple"),
        ylab = "Tridecanal", xlab =  "", las=2,
        ylim=c(0,75),
        cex.lab=1.5, cex.axis=1.5)#***
text(1, 73.6, "a", cex=1.5)
text(2, 73.6, "b", cex=1.5)
text(3, 73.6, "a", cex=1.5)
text(4, 73.6, "ab", cex=1.5)
boxplot(trans.Geranylacetone~Species,data=sig, col = c("red", "blue", "darkgreen", "purple"),
        ylab = "trans-Geranylacetone", xlab =  "", las=2,
        ylim=c(0,18),
        cex.lab=1.5, cex.axis=1.5)#***
text(1, 17.6, "a", cex=1.5)
text(2, 17.6, "b", cex=1.5)
text(3, 17.6, "ab", cex=1.5)
text(4, 17.6, "ab", cex=1.5)

#DO FOR ALL IN DUCKWEED RENAMED AND DO LETTERS ON TABLE, SHOW MORE DATA
#dw <- read.csv("Flavour_NEW2_mixedRT_PPB_low_in_dw_rename.csv")
dw <- read.csv("Flavour_NEW2_mixedRT_PPB_low_in_dw+Ljp_rename.csv")
library(dplyr)
dw <- filter(dw, Species != "Coriander") %>% droplevels()
dw <- filter(dw, Species != "Basil") %>% droplevels()
dw <- filter(dw, Species != "Dandelion") %>% droplevels()
dw <- filter(dw, Species != "Spinach") %>% droplevels()

pairwise.wilcox.test(dw$Acetaldehyde, dw$Species,
                     p.adjust.method = "BH") #ns
pairwise.wilcox.test(dw$Dimethyl.sulfide, dw$Species,
                     p.adjust.method = "BH")
pairwise.wilcox.test(dw$X2.methyl.1.propanal, dw$Species,
                     p.adjust.method = "BH")
pairwise.wilcox.test(dw$Methyl.acetate, dw$Species,
                     p.adjust.method = "BH")
pairwise.wilcox.test(dw$X3.methyl.furan, dw$Species,
                     p.adjust.method = "BH")
pairwise.wilcox.test(dw$Methyl.ethyl.ketone, dw$Species,
                     p.adjust.method = "BH")
pairwise.wilcox.test(dw$Butanal..3.methyl., dw$Species,
                     p.adjust.method = "BH") #L mino L minu
pairwise.wilcox.test(dw$Ethanol, dw$Species,
                     p.adjust.method = "BH")
pairwise.wilcox.test(dw$X2.ethyl.furan, dw$Species,
                     p.adjust.method = "BH") #L minu diff L mino, L turion diff all, S poly diff L jp L tu
pairwise.wilcox.test(dw$X3.pentanone, dw$Species,
                     p.adjust.method = "BH")
pairwise.wilcox.test(dw$X1.penten.3.one, dw$Species,
                     p.adjust.method = "BH")
pairwise.wilcox.test(dw$Toluene, dw$Species,
                     p.adjust.method = "BH")
pairwise.wilcox.test(dw$X2.butenal, dw$Species,
                     p.adjust.method = "BH")
pairwise.wilcox.test(dw$X2.3.pentanedione, dw$Species,
                     p.adjust.method = "BH") #spoly all
pairwise.wilcox.test(dw$X1.cyclopropylpropan.1.one, dw$Species,
                     p.adjust.method = "BH")
pairwise.wilcox.test(dw$Hexanal, dw$Species,
                     p.adjust.method = "BH")
pairwise.wilcox.test(dw$X1.undecanol, dw$Species,
                     p.adjust.method = "BH")
pairwise.wilcox.test(dw$trans.2.Pentenal, dw$Species,
                     p.adjust.method = "BH") # L minu diff Ljp, Lmino, L tu diff L jp, mino, S poly diff to L jp L mino
pairwise.wilcox.test(dw$X1.penten.3.ol, dw$Species,
                     p.adjust.method = "BH")
pairwise.wilcox.test(dw$Heptanal, dw$Species,
                     p.adjust.method = "BH") # L jp diff L mino, L minu diff L jp, L mino, L tu diff L jp L mino, S poly diff L mino
pairwise.wilcox.test(dw$p.Xylene, dw$Species,
                     p.adjust.method = "BH")
pairwise.wilcox.test(dw$X2.hexenal, dw$Species,
                     p.adjust.method = "BH") # L tu diff to L jp, mino, minu, S poly diff L jp L mino
pairwise.wilcox.test(dw$X2.pentylfuran, dw$Species,
                     p.adjust.method = "BH") # L mino diff L jp, L minu diff L mino, S poly diff L mino
pairwise.wilcox.test(dw$X3.octanone, dw$Species,
                     p.adjust.method = "BH")
pairwise.wilcox.test(dw$Styrene, dw$Species,
                     p.adjust.method = "BH")
pairwise.wilcox.test(dw$Cumene, dw$Species,
                     p.adjust.method = "BH")

pairwise.wilcox.test(dw$cis.2.pentenol, dw$Species,
                     p.adjust.method = "BH")
pairwise.wilcox.test(dw$Vinyl.hexanoate, dw$Species,
                     p.adjust.method = "BH")
pairwise.wilcox.test(dw$Cistus.cyclohexanone, dw$Species,
                     p.adjust.method = "BH")
pairwise.wilcox.test(dw$trans.2.heptenal, dw$Species,
                     p.adjust.method = "BH")
pairwise.wilcox.test(dw$Methyl.heptenone, dw$Species,
                     p.adjust.method = "BH") # L minu from L mino, S poly from L mino
pairwise.wilcox.test(dw$Mesitylene, dw$Species,
                     p.adjust.method = "BH")
pairwise.wilcox.test(dw$Hexanol, dw$Species,
                     p.adjust.method = "BH")
pairwise.wilcox.test(dw$Benzene..2.propenyl., dw$Species,
                     p.adjust.method = "BH")
pairwise.wilcox.test(dw$cis.3.hexen.1.ol, dw$Species,
                     p.adjust.method = "BH") #L mino L jp, L tu diff L jp L mino L minu
pairwise.wilcox.test(dw$trans.2.hexen.1.ol, dw$Species,
                     p.adjust.method = "BH") # L minu from L jp, mino, L tu from L minu, S poly from L jp, mino, tu
pairwise.wilcox.test(dw$L.fenchone, dw$Species,
                     p.adjust.method = "BH")
pairwise.wilcox.test(dw$Benzene..1.3.bis.1.1.dimethylethyl.., dw$Species,
                     p.adjust.method = "BH")
pairwise.wilcox.test(dw$Ethylcyclopent.1.carboxaldehydee, dw$Species,
                     p.adjust.method = "BH")
pairwise.wilcox.test(dw$X1.octen.3.ol, dw$Species,
                     p.adjust.method = "BH") # L minu from L mino, S poly from L jp L mino L tu
pairwise.wilcox.test(dw$Decanal, dw$Species,
                     p.adjust.method = "BH")
pairwise.wilcox.test(dw$Pyrrole, dw$Species,
                     p.adjust.method = "BH") # L  tu diff L jp, mino, minu, S poly diff L jp, mino, tu
pairwise.wilcox.test(dw$D.camphor, dw$Species,
                     p.adjust.method = "BH")
pairwise.wilcox.test(dw$Benzaldehyde, dw$Species,
                     p.adjust.method = "BH")
pairwise.wilcox.test(dw$Linalyl.anthranilate, dw$Species,
                     p.adjust.method = "BH") # L mino from L jp, L minu L jp, L tu from L minu, S poly from L mino L minu
pairwise.wilcox.test(dw$Beta.cyclocitral, dw$Species,
                     p.adjust.method = "BH")
pairwise.wilcox.test(dw$trans.2.decenal, dw$Species,
                     p.adjust.method = "BH") # S poly L mino
pairwise.wilcox.test(dw$Naphthalene, dw$Species,
                     p.adjust.method = "BH")
pairwise.wilcox.test(dw$Tridecanal, dw$Species,
                     p.adjust.method = "BH") # L minu from L jp, L mino, L tu from L minu
pairwise.wilcox.test(dw$trans.geranylacetone, dw$Species,
                     p.adjust.method = "BH") #L minu from L jp, L mino
pairwise.wilcox.test(dw$X2.methyl.naphthalene, dw$Species,
                     p.adjust.method = "BH")
pairwise.wilcox.test(dw$Tetradecanal, dw$Species,
                     p.adjust.method = "BH") # L minu from L jp, L mino, L tu from Lmu, S poly from L minu
pairwise.wilcox.test(dw$trans.beta.ionone, dw$Species,
                     p.adjust.method = "BH")
pairwise.wilcox.test(dw$Pentadecanal, dw$Species,
                     p.adjust.method = "BH") # L minu from L jp, L mino, L tu from L mu
pairwise.wilcox.test(dw$cis.11.tetradecenal, dw$Species,
                     p.adjust.method = "BH") # L minu from L jp, L mo, L tu from L mu
pairwise.wilcox.test(dw$X2.tetradecenal, dw$Species,
                     p.adjust.method = "BH") # L mu from L mo L jp, L tu from L mu, S piky from L mu
pairwise.wilcox.test(dw$X2.ethyl.3.methyl.maleimide, dw$Species,
                     p.adjust.method = "BH")
pairwise.wilcox.test(dw$cis.cis.7.10..Hexadecadienal, dw$Species,
                     p.adjust.method = "BH") # L mu from L jp L mu, L tu from L mu, S poly from L mu

sig <- dw

#with Ljp 21 sig diff species
# Trans 2 decanal # S poly L mino

#order of species
sig$Species <- factor(sig$Species, levels = c("L. minor", "L. japonica", "L. turionifera", "L. minuta", "S. polyrhiza"))
cols <- c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00")
#replot
#new names as in excel sheet
#plots 1
#par("mar")
#remake boxplots using ggplot
library(ggplot2)
one <- ggplot(sig, aes(x=X2.pentylfuran, y=Species, group=Species)) + 
  geom_boxplot(aes(fill=Species))
las=1
one

#save pdf via console
par(mar=c(1,1,1,1))
#saving as tiff not working
#tiff('Flavour_species_boxplotssig_Ljp_22comp.tiff', units="px", width=1000, height=600, res=300, compression = 'lzw')
pdf('Flavour_species_boxplotssig_Ljp_22comp_testpt1.pdf', width=800, height=800)
par(mfrow = c(3,4),  mar=c(5,4.5,4,2))
#par(mar=c(1,1,1,1))
boxplot(X2.pentylfuran~Species,data=sig, col = c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        ylab = "2-pentylfuran", xlab =  "", las=2,
        ylim=c(0,100),
        cex.lab=1.1, cex.axis=1.1)#***
text(1, 100, "a", cex=1.1)
text(2, 100, "b", cex=1.1)
text(3, 100, "ab", cex=1.1)
text(4, 100, "b", cex=1.1)
text(5, 100, "b", cex=1.1)
# Pentyl furan new L mino diff L jp, L minu diff L mino, S poly diff L mino
stripchart(X2.pentylfuran ~ Species, sig,  at =c(1,2,3,4,5),           # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.5,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)
boxplot(X2.ethyl.furan~Species,data=sig, col = c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        ylab = "2-Ethylfuran", xlab =  "", las=2,
        ylim=c(0,15),
        cex.lab=1.1, cex.axis=1.1)#***
text(1, 14.5, "b", cex=1.1)
text(2, 14.5, "ab", cex=1.1)
text(3, 14.5, "c", cex=1.1)
text(4, 14.5, "a", cex=1.1)
text(5, 14.5, "ab", cex=1.1)
stripchart(X2.ethyl.furan ~ Species, sig,  at =c(1,2,3,4,5),           # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.5,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)
#L minu diff L mino, L turion diff all, S poly diff L jp L tu
boxplot(cis.3.hexen.1.ol~Species,data=sig, col = c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        ylab = "cis-3-Hexen-1-ol", xlab =  "", las=2,
        ylim=c(0,220),
        cex.lab=1.1, cex.axis=1.1)#***
text(1, 219, "a", cex=1.1)
text(2, 219, "b", cex=1.1)
text(3, 219, "c", cex=1.1)
text(4, 219, "ab", cex=1.1)
text(5, 219, "ab", cex=1.1)
#L mino L jp, L tu diff L jp L mino L minu
stripchart(cis.3.hexen.1.ol ~ Species, sig,  at =c(1,2,3,4,5),           # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.5,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)
boxplot(Butanal..3.methyl.~Species,data=sig, col = c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        ylab = "Butanal-3-methyl", xlab =  "", las=2,
        ylim=c(0,150),
        cex.lab=1.1, cex.axis=1.1)#***
text(1, 145, "b", cex=1.1)
text(2, 145, "ab", cex=1.1)
text(3, 145, "ab", cex=1.1)
text(4, 145, "a", cex=1.1)
text(5, 145, "ab", cex=1.1)
stripchart(Butanal..3.methyl. ~ Species, sig,  at =c(1,2,3,4,5),           # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.5,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)
boxplot(X2.3.pentanedione~Species,data=sig, col = c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        ylab = "2,3-Pentanedione", xlab =  "", las=2,
        ylim=c(0,40),
        cex.lab=1.1, cex.axis=1.1)#***
text(1, 39.5, "a", cex=1.1)
text(2, 39.5, "a", cex=1.1)
text(3, 39.5, "a", cex=1.1)
text(4, 39.5, "a", cex=1.1)
text(5, 39.5, "b", cex=1.1)
stripchart(X2.3.pentanedione ~ Species, sig,  at =c(1,2,3,4,5),           # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.5,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)
#spoly all
boxplot(trans.2.hexen.1.ol~Species,data=sig, col = c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        ylab = "trans-2-Hexen-1-ol", xlab =  "", las=2,
        ylim=c(0,80),
        cex.lab=1.1, cex.axis=1.1)#***
text(1, 79, "bc", cex=1.1)
text(2, 79, "c", cex=1.1)
text(3, 79, "c", cex=1.1)
text(4, 79, "ab", cex=1.1)
text(5, 79, "a", cex=1.1)
# L minu from L jp, mino, L tu from L minu, S poly from L jp, mino, tu
stripchart(trans.2.hexen.1.ol ~ Species, sig,  at =c(1,2,3,4,5),           # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.5,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)
boxplot(X2.hexenal~Species,data=sig, col = c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        ylab = "2-hexenal", xlab =  "", las=2,
        ylim=c(0,500),
        cex.lab=1.1, cex.axis=1.1)#***
text(1, 490, "b", cex=1.1)
text(2, 490, "b", cex=1.1)
text(3, 490, "c", cex=1.1)
text(4, 490, "b", cex=1.1)
text(5, 490, "a", cex=1.1)
# L tu diff to L jp, mino, minu, S poly diff L jp L mino
stripchart(X2.hexenal ~ Species, sig,  at =c(1,2,3,4,5),           # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.5,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)
boxplot(Linalyl.anthranilate~Species,data=sig, col = c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        ylab = "Linalyl anthranilate", xlab =  "", las=2,
        ylim=c(0,12),
        cex.lab=1.1, cex.axis=1.1)#***
text(1, 11.7, "a", cex=1.1)
text(2, 11.7, "b", cex=1.1)
text(3, 11.7, "b", cex=1.1)
text(4, 11.7, "a", cex=1.1)
text(5, 11.7, "b", cex=1.1)
# L mino from L jp, L minu L jp, L tu from L minu, S poly from L mino L minu
stripchart(Linalyl.anthranilate ~ Species, sig,  at =c(1,2,3,4,5),           # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.5,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)
boxplot(X2.tetradecenal~Species,data=sig, col = c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        ylab = "2-tetradecenal", xlab =  "", las=2,
        ylim=c(0,8),
        cex.lab=1.1, cex.axis=1.1)#***
text(1, 7.8, "a", cex=1.1)
text(2, 7.8, "a", cex=1.1)
text(3, 7.8, "a", cex=1.1)
text(4, 7.8, "b", cex=1.1)
text(5, 7.8, "a", cex=1.1)
# L mu from L mo L jp, L tu from L mu, S piky from L mu
stripchart(X2.tetradecenal ~ Species, sig,  at =c(1,2,3,4,5),           # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.5,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)
boxplot(cis.cis.7.10..Hexadecadienal~Species,data=sig, col = c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        ylab = "cis,cis-7,10,-Hexadecadienal", xlab =  "", las=2,
        ylim=c(0,20),
        cex.lab=1, cex.axis=1.1)#***
text(1, 19.8, "a", cex=1.1)
text(2, 19.8, "a", cex=1.1)
text(3, 19.8, "a", cex=1.1)
text(4, 19.8, "b", cex=1.1)
text(5, 19.8, "a", cex=1.1)
# L mu from L jp L mu, L tu from L mu, S poly from L mu
stripchart(cis.cis.7.10..Hexadecadienal ~ Species, sig,  at =c(1,2,3,4,5),           # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.5,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)
boxplot(X1.octen.3.ol~Species,data=sig, col = c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        ylab = "1-octen-3-ol", xlab =  "", las=2,
        ylim=c(0,120),
        cex.lab=1.1, cex.axis=1.1)#***
text(1, 119.6, "a", cex=1.1)
text(2, 119.6, "ab", cex=1.1)
text(3, 119.6, "ab", cex=1.1)
text(4, 119.6, "b", cex=1.1)
text(5, 119.6, "c", cex=1.1)
# L minu from L mino, S poly from L jp L mino L tu
stripchart(X1.octen.3.ol ~ Species, sig,  at =c(1,2,3,4,5),           # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.5,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)
boxplot(trans.2.Pentenal~Species,data=sig, col = c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        ylab = "trans-2-pentenal", xlab =  "", las=2,
        ylim=c(0,200),
        cex.lab=1.1, cex.axis=1.1)#***
text(1, 199.6, "a", cex=1.1)
text(2, 199.6, "a", cex=1.1)
text(3, 199.6, "b", cex=1.1)
text(4, 199.6, "b", cex=1.1)
text(5, 199.6, "b", cex=1.1)
#L minu diff Ljp, Lmino, L tu diff L jp, mino, S poly diff to L jp L mino
stripchart(trans.2.Pentenal ~ Species, sig,  at =c(1,2,3,4,5),           # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.5,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)
#install.packages("grDevices")
#library(savePlot)
#saveplot("test1.pdf", type = "pdf") 
#library(ggplot2)
#ggsave("Flav_boxplots_pt1.pdf")
dev.off()
#pdf('Flavour_species_boxplotssig_Ljp_22comp_pt2.pdf', width=1400, height=800)
par(mfrow = c(3,4),  mar=c(5,4.5,4,2))
boxplot(trans.2.decenal~Species,data=sig, col = c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        ylab = "Trans-2-decenal", xlab =  "", las=2,
        ylim=c(0,20),
        cex.lab=0.5, cex.axis=1.1)#***
text(1, 20, "a", cex=1.1)
text(2, 20, "ab", cex=1.1)
text(3, 20, "ab", cex=1.1)
text(4, 20, "ab", cex=1.1)
text(5, 20, "bc", cex=1.1)
# S poly L mino
stripchart(trans.2.decenal ~ Species, sig,  at =c(1,2,3,4,5),           # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.5,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)
boxplot(Heptanal~Species,data=sig, col = c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        ylab = "Heptanal", xlab =  "", las=2,
        ylim=c(0,40),
        cex.lab=1.1, cex.axis=1.1)#***
text(1, 38.6, "a", cex=1.1)
text(2, 38.6, "b", cex=1.1)
text(3, 38.6, "c", cex=1.1)
text(4, 38.6, "c", cex=1.1)
text(5, 38.6, "bc", cex=1.1)
# L jp diff L mino, L minu diff L jp, L mino, L tu diff L jp L mino, S poly diff L mino
stripchart(Heptanal ~ Species, sig,  at =c(1,2,3,4,5),           # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.5,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)
boxplot(Methyl.heptenone~Species,data=sig, col = c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        ylab = "Methyl heptenone", xlab =  "", xaxt="n", las=2,
        ylim=c(0,18),
        cex.lab=1.1, cex.axis=1.1)#***
text(1, 17.6, "a", cex=1.1)
text(2, 17.6, "ab", cex=1.1)
text(3, 17.6, "ab", cex=1.1)
text(4, 17.6, "b", cex=1.1)
text(5, 17.6, "b", cex=1.1)
# L minu from L mino, S poly from L mino
stripchart(Methyl.heptenone ~ Species, sig,  at =c(1,2,3,4,5),           # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.5,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)
boxplot(Ethylcyclopent.1.carboxaldehyde~Species,data=sig, col = c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        ylab = "5-Ethylcyclopentene-1-carboxaldehyde", xlab =  "", xaxt="n", las=2,
        ylim=c(0,130),
        cex.lab=1.1, cex.axis=1.1)#***
text(1, 128.6, "a", cex=1.1)
text(2, 128.6, "b", cex=1.1)
text(3, 128.6, "c", cex=1.1)
text(4, 128.6, "c", cex=1.1)
text(5, 128.6, "c", cex=1.1)
# l mino diff to ljp, l tu diff to l mino l jp, l minu diff to l mino l jp, s poly diff to l mino l jp l turi
stripchart(Ethylcyclopent.1.carboxaldehyde ~ Species, sig,  at =c(1,2,3,4,5),           # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.5,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)
boxplot(cis.11.tetradecenal~Species,data=sig, col = c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        ylab = "cis-11-tetradecenal", xlab =  "", xaxt="n", las=2,
        ylim=c(0,25),
        cex.lab=1.1, cex.axis=1.1)#***
text(1, 24.6, "a", cex=1.1)
text(2, 24.6, "a", cex=1.1)
text(3, 24.6, "a", cex=1.1)
text(4, 24.6, "b", cex=1.1)
text(5, 24.6, "ab", cex=1.1)
# L minu from L jp, L mo, L tu from L mu
stripchart(cis.11.tetradecenal ~ Species, sig,  at =c(1,2,3,4,5),           # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.5,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)
boxplot(Tetradecanal~Species,data=sig, col = c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        ylab = "Tetradecanal", xlab =  "", xaxt="n", las=2,
        ylim=c(0,120),
        cex.lab=1.1, cex.axis=1.1)#***
text(1, 118.6, "a", cex=1.1)
text(2, 118.6, "a", cex=1.1)
text(3, 118.6, "a", cex=1.1)
text(4, 118.6, "b", cex=1.1)
text(5, 118.6, "a", cex=1.1)
# L minu from L jp, L mino, L tu from Lmu, S poly from L minu
stripchart(Tetradecanal ~ Species, sig,  at =c(1,2,3,4,5),           # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.5,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)
boxplot(Pentadecanal~Species,data=sig, col = c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        ylab = "Pentadecanal", xlab =  "", xaxt="n", las=2,
        ylim=c(0,280),
        cex.lab=1.1, cex.axis=1.1)#***
text(1, 278.6, "a", cex=1.1)
text(2, 278.6, "a", cex=1.1)
text(3, 278.6, "a", cex=1.1)
text(4, 278.6, "b", cex=1.1)
text(5, 278.6, "a", cex=1.1)
stripchart(Pentadecanal ~ Species, sig,  at =c(1,2,3,4,5),           # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.5,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)
boxplot(Pyrrole~Species,data=sig, col = c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        ylab = "Pyrrole", xlab =  "", xaxt="n", las=2,
        ylim=c(0,10),
        cex.lab=1.1, cex.axis=1.1)#***
text(1, 9.6, "b", cex=1.1)
text(2, 9.6, "b", cex=1.1)
text(3, 9.6, "a", cex=1.1)
text(4, 9.6, "bc", cex=1.1)
text(5, 9.6, "c", cex=1.1)
# L  tu diff L jp, mino, minu, S poly diff L jp, mino, tu
stripchart(Pyrrole ~ Species, sig,  at =c(1,2,3,4,5),           # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.5,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)
boxplot(Tridecanal~Species,data=sig, col = c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        ylab = "Tridecanal", xlab =  "", xaxt="n", las=2,
        ylim=c(0,75),
        cex.lab=1.1, cex.axis=1.1)#***
text(1, 73.6, "a", cex=1.1)
text(2, 73.6, "a", cex=1.1)
text(3, 73.6, "a", cex=1.1)
text(4, 73.6, "b", cex=1.1)
text(5, 73.6, "ab", cex=1.1)
# L minu from L jp, L mino, L tu from L minu
stripchart(Tridecanal ~ Species, sig,  at =c(1,2,3,4,5),           # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.5,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)
boxplot(trans.geranylacetone~Species,data=sig, col = c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        ylab = "trans-Geranylacetone", xlab =  "", xaxt="n", las=2,
        ylim=c(0,18),
        cex.lab=1.1, cex.axis=1.1)#***
text(1, 17.6, "a", cex=1.1)
text(2, 17.6, "a", cex=1.1)
text(3, 17.6, "ab", cex=1.1)
text(4, 17.6, "b", cex=1.1)
text(5, 17.6, "ab", cex=1.1)
#L minu from L jp, L mino
stripchart(trans.geranylacetone ~ Species, sig,  at =c(1,2,3,4,5),           # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.5,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)
ggsave("Flav_boxplots_pt2.pdf")
#dev.off()
#dev.new()

