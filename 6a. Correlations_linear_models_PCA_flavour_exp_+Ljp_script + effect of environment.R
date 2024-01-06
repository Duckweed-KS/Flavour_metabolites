#combined flavour experiment
setwd("C:\\Users\\kelli\\OneDrive - The University of Nottingham\\Nottingham\\Year 2\\UK Lemna - Bradford - Hull\\Flavour analysis\\HS")

forpca <- read.csv("") 
#dat <- read.csv("Comb_hs+photos+biomass+lightlevs+flavour+met_summ.csv")

#removed flavour as old data, replaced hs summaries, species updated
#dat <- read.csv("Comb_hs+photos+biomass+lightlevs_summ_spchange.csv")
#correlates light levels, photos and fieldspec with measurements, mets

#for mets
dat <- read.csv("Comb_hs+photos+biomass+lightlevs_summ_spchange+Ljp.csv")

#for  flavs
#correlates light levels, photos and fieldspec with measurements, flavs
dat <- read.csv("Comb_hs+photos+biomass+lightlevs_summ_spchange+Ljp_flavs.csv")

names(dat) 

#corr plot

str(dat)
sapply(dat, class)

#non automated approach
names(dat)
comb <-dat
#correlations fieldspec with real
corr <- cor.test(comb$Green_area_mean, comb$X.Area_mean,
                 method = "pearson"
)
corr$estimate #pos
corr$p.value #sig
#correlations fieldspec with real
corr <- cor.test(comb$Green_area_mean, comb$NDVI_mean,
                 method = "pearson"
)
corr$estimate #pos
corr$p.value #sig
corr <- cor.test(comb$Green_area_mean, comb$GM_mean,
                 method = "pearson"
)
corr$estimate #pos
corr$p.value #sig
corr <- cor.test(comb$Green_area_mean, comb$GI_mean,
                 method = "pearson"
)
corr$estimate #pos
corr$p.value #sig

#all greenness correlates with green area

#growth biomass correlartes with green area
corr <- cor.test(comb$Green_area_mean, comb$FW.g.norm,
                 method = "pearson"
)
corr$estimate #pos
corr$p.value #sig
corr <- cor.test(comb$Green_area_mean, comb$g.6.months_norm,
                 method = "pearson"
)
corr$estimate #pos
corr$p.value #sig
corr <- cor.test(comb$Green_area_mean, comb$FDW..g.,
                 method = "pearson"
)
corr$estimate #weak pos
corr$p.value #sig
corr <- cor.test(comb$Green_area_mean, comb$FDW.g.6.months,
                 method = "pearson"
)
corr$estimate #weak pos
corr$p.value #sig

#good correlation with fw (6 month better) but fdw weaker

#percent area with growth biomass
#growth biomass correlartes with %area
corr <- cor.test(comb$X.Area_mean, comb$FW.g.norm,
                 method = "pearson"
)
corr$estimate #pos
corr$p.value #sig
corr <- cor.test(comb$X.Area_mean, comb$g.6.months_norm,
                 method = "pearson"
)
corr$estimate #pos
corr$p.value #sig
corr <- cor.test(comb$X.Area_mean, comb$FDW..g.,
                 method = "pearson"
)
corr$estimate #weak pos
corr$p.value #sig
corr <- cor.test(comb$X.Area_mean, comb$FDW.g.6.months,
                 method = "pearson"
)
corr$estimate #weak pos
corr$p.value #sig

#greenness with fieldspec param
comb$Species <- factor(comb$Species, levels = c("L. minor", "L. japonica", "L. turionifera", "L. minuta", "S. polyrhiza"))


#best 2 to plot
plot(comb$Green_area_mean ~ comb$NDVI_mean, 
     pch = 19, col = factor(comb$Species))

plot(comb$X.Area_mean ~ comb$FW.g.norm,
     pch = 19, col = factor(comb$Species))

library(ggplot2)
library(ggpubr)
library(ggrepel)

comb$Species <- factor(dat$Species, levels = c("L. minor", "L. japonica", "L. turionifera", "L. minuta", "S. polyrhiza"))


#do this for fieldspec datasets
#create scatterplot with points colored by group
one <- ggplot(comb, aes(Green_area_mean, NDVI_mean)) +
  geom_point(aes(color=Species), size=3) +
  geom_smooth(method = "lm", col = "black") +
  xlab("NDVI") +
  ylab("Green area") +
  scale_color_manual(values=c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"))  +
  theme(legend.text=element_text(size=12, face="italic"))+
theme_bw() +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#do this for measured data
comb$name <- (comb$Accession)
#create scatterplot with points colored by group
two <- ggplot(comb, aes(X.Area_mean, FW.g.norm)) +
  geom_point(aes(color=Species), size=3) +
  geom_smooth(method = "lm", col = "black") +
  stat_cor(method = "pearson", label.x = 0, label.y = 400) +
  xlab("% Coverage area") +
  ylab("Fresh weight (g)") +
  scale_color_manual(values=c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"))  +
  geom_text_repel(data=subset(comb, X.Area_mean > 70 | FW.g.norm > 350),
            aes(X.Area_mean,FW.g.norm,label="KS12")) +
  geom_text_repel(data=subset(comb, X.Area_mean < 10 | FW.g.norm < 100),
                  aes(X.Area_mean,FW.g.norm,label=name)) +
  theme(legend.text=element_text(size=12, face="italic"))+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
two

#create scatterplot with points colored by group
three <- ggplot(comb, aes(Green_area_mean, X.Area_mean)) +
  geom_point(aes(color=Species), size=3) +
  geom_smooth(method = "lm", col = "black") +
  stat_cor(method = "pearson", label.x = 100, label.y = 68) +
  xlab("Green area") +
  ylab("% Coverage area") +
  scale_color_manual(values=c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"))  +
  geom_text_repel(data=subset(comb,Green_area_mean > 150 | X.Area_mean > 70),
                  aes(Green_area_mean,X.Area_mean,label="KS12")) +
  geom_text_repel(data=subset(comb, Green_area_mean < 110 | X.Area_mean < 10),
                  aes(Green_area_mean,X.Area_mean,label=name)) +
  theme(legend.text=element_text(size=12, face="italic"))+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
three

threet <- ggplot(comb, aes(GM_mean, Green_area_mean)) +
  geom_point(aes(color=Species), size=3) +
  geom_smooth(method = "lm", col = "black") +
  stat_cor(method = "pearson", label.x = 0, label.y = 200) +
  xlab("Green model") +
  ylab("Green area") +
  scale_color_manual(values=c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"))  +
  geom_text_repel(data=subset(comb, Green_area_mean > 160),
                  aes(GM_mean,Green_area_mean,label=name)) +
  geom_text_repel(data=subset(comb, GM_mean < 0.5 | Green_area_mean < 100),
                  aes(GM_mean,Green_area_mean,label=name)) +
  theme(legend.text=element_text(size=12, face="italic"))+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
threet

#changed from PRI to OVI
fourt <- ggplot(comb, aes(OVI_mean, FW.g.norm)) +
  geom_point(aes(color=Species), size=3) +
  geom_smooth(method = "lm", col = "black") +
  stat_cor(method = "pearson", label.x = 0.95, label.y = 10) +
  xlab("OVI mean") +
  ylab("Total fresh weight (g)") +
  scale_color_manual(values=c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"))  +
  geom_text_repel(data=subset(comb, OVI_mean < -0.1 | FW.g.norm > 350),
                aes(OVI_mean,FW.g.norm,label=name)) +
  geom_text_repel(data=subset(comb,OVI_mean < -0.3 | FW.g.norm < 100),
                  aes(OVI_mean,FW.g.norm,label=name)) +
  theme(legend.text=element_text(size=12, face="italic"))+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
fourt

library(gridExtra)
pdf('Scatterplots_GreenareaNDVI+CovareaFW_+Ljp.pdf', width=22, height=13)
#tiff('Scatterplots_GreenareaNDVI+CovareaFW_+Ljp.tiff', units="in", width=22, height=13, res=300, compression = 'lzw')
grid.arrange(one, two, ncol=2, nrow = 2)

pdf('Scatterplots_CovareaFW+Greenareacov_+Ljp.pdf', width=22, height=13)
#tiff('Scatterplots_CovareaFW+Greenareacov_+Ljp.tiff', units="in", width=22, height=13, res=300, compression = 'lzw')
grid.arrange(two, three, ncol=3, nrow = 2)
dev.off()
#dev.new()

pdf('Scatterplots_GMgreenarea_OVI_FW+Ljp1.pdf', width=22, height=13)
#tiff('Scatterplots_GMgreenarea_PRI_FW+Ljp.tiff', units="in", width=22, height=13, res=300, compression = 'lzw')
grid.arrange(threet, fourt, ncol=3, nrow = 2)
dev.off()
#dev.new()

#try with raw reps so more points on graph, but low correls
#no of raw reps for NDVI = 4
#no of raw reps for green area = 4 or 40
#no of raw reps for % cov area = 4 or 40
#no of raw reps for fw = 4

#read in Biomass
bio <- read.csv("Biomass.csv")
#read in Green_area_%cov_summary
gacov <- read.csv("Green_area_%cov_summary.csv")
#read in HS_param_summ_accrep
ndvi <- read.csv("HS_param_summ_accrep.csv")
#include light readings, metabolites, flavour in join
li <- read.csv("Light_spectro_spchange.csv") #100 obs
me <- read.csv("Metabolite_all_spnew.csv") #100 obs added ks06b in
fl <- read.csv("Flavour_new_minforjoin.csv") #made to be 100 obs

#attach
join <- cbind(bio, gacov, ndvi)
join1 <- cbind(bio, gacov, ndvi, li, me, fl)
#might need to use join function instead
#load in ordered by accession_rep so all right to start

#arrange by accession_rep
library(dplyr)
j1 <- join1[order(join1$Accession_rep),]

#numeric only?
numeric_cols <- lapply(j1, is.numeric) %>% unlist() %>% unname()
j2 <- j1[, numeric_cols] #jus tget numeric
numeric_cols #check which numeric

#remove pointless cols
#remove anything with _n or stdev
j3 <- j2[, -grep("\\_n$", colnames(j2))]
j3 <- j3[, -grep("\\_stdev$", colnames(j3))]
#remove cols by index so just meaningful data included
j4 <- j3[-c(1,2,11,16,35,36,37,58,59,104)]

#fails at tryptamine - remove
j4 <- j4[-c(77)]

j4$Species <- join1$Species

join <- j4

#do i really need to analyse everything all together?
#just do quick linear model output rather than all graphs?
setwd("C:\\Users\\kelli\\OneDrive - The University of Nottingham\\Nottingham\\Year 2\\UK Lemna - Bradford - Hull\\Flavour analysis\\HS\\lm test massive")
lms <- expand.grid(1:186, 1:186)
lms_names <- expand.grid(names(j4)[1:186], names(j4)[1:186])
out <- vector(mode = "list", length = nrow(lms))
for(i in 1:nrow(lms)){
  lms_col_2 <- lms$Var2[i]
  lms_col_1 <- lms$Var1[i]
  plot_name <- paste0(stringr::str_pad(i, width = 3, pad = "0"), " ", lms_names$Var2[i], " vs. ", lms_names$Var1[i], ".png")
  png(plot_name, width = 500, height = 500, type = "cairo")
  plot(j4[, lms_col_1], j4[, lms_col_2], main = paste0(lms_names$Var2[i], " vs. ", lms_names$Var1[i]), type = "n")
  text(j4[, lms_col_1], j4[, lms_col_2], row.names(j4), cex = 0.8)
  abline(lm(j4[, lms_col_2] ~ j4[, lms_col_1]))
  #plots per column and output as png with correspinding names
  dev.off()
  out[[i]] <- data.frame(r.squared = summary(lm(j4[, lms_col_2] ~ j4[, lms_col_1]))$r.squared,
                         pos_neg = ifelse(summary(lm(j4[, lms_col_2] ~ j4[, lms_col_1]))$coef[2, 1] > 0, "+", "-"),
                         p.value = summary(lm(j4[, lms_col_2] ~ j4[, lms_col_1]))$coefficients[2, 4])
}
#select r squared and p val from table and make df output
#use estimate col as gradient to make + or - col

#too many to output?

#write to files
library(dplyr)
(all_output <- bind_cols(lms_names, do.call(rbind, out)))
all_output <- all_output[all_output$r.squared != 1, ]
all_output <- all_output[all_output$p.value < 0.05, ]
all_output$r.squared <- ifelse(all_output$pos_neg == "+", all_output$r.squared, -all_output$r.squared)

write.csv(all_output, "Lm_100obs_correlations.csv")

#analyse in more targetted aproach?
corr <- cor.test(join$GM_mean, join$NDVI_mean,
                 method = "pearson"
)
corr$estimate #0.94
corr$p.value #sig
corr <- cor.test(join$GI_mean, join$NDVI_mean,
                 method = "pearson"
)
corr$estimate #0.94
corr$p.value #sig
corr <- cor.test(join$GI_mean, join$GM_mean,
                 method = "pearson"
)
corr$estimate #0.94
corr$p.value #sig

#to edit after
corr <- cor.test(join$Green_area_mean, join$NDVI_mean,
                 method = "pearson"
)
corr$estimate #weak
corr$p.value #ns
corr <- cor.test(join$X.Area_mean, join$Tot_FW,
                 method = "pearson"
)
corr$estimate #weak
corr$p.value #sig

plot(join$Green_area_mean ~ join$NDVI_mean, 
     pch = 19, col = factor(join$Species))
plot(join$X.Area_mean ~ join$Tot_FW, 
     pch = 19, col = factor(join$Species))

library(ggplot2)
library(ggpubr)
library(ggrepel)

four <- ggplot(join, aes(PRI_mean, Tot_FW)) +
  geom_point(aes(color=Species), size=3) +
  geom_smooth(method = "lm", col = "black") +
  stat_cor(method = "pearson", label.x = 0, label.y = 200) +
  xlab("PRI mean") +
  ylab("Total fresh weight (g)") +
  scale_color_manual(values=c('red', 'cyan', 'green', 'purple'))  +
  #geom_text_repel(data=subset(comb,Green_area_mean > 150 | X.Area_mean > 70),
  #                aes(Green_area_mean,X.Area_mean,label="KS12")) +
  #geom_text_repel(data=subset(comb, Green_area_mean < 110 | X.Area_mean < 10),
  #                aes(Green_area_mean,X.Area_mean,label=name)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
four

five <- ggplot(join, aes(GM_mean, Green_area_mean)) +
  geom_point(aes(color=Species), size=3) +
  geom_smooth(method = "lm", col = "black") +
  stat_cor(method = "pearson", label.x = 0, label.y = 200) +
  xlab("Green model") +
  ylab("Green area") +
  scale_color_manual(values=c('red', 'cyan', 'green', 'purple'))  +
  #geom_text_repel(data=subset(comb,Green_area_mean > 150 | X.Area_mean > 70),
  #                aes(Green_area_mean,X.Area_mean,label="KS12")) +
  #geom_text_repel(data=subset(comb, Green_area_mean < 110 | X.Area_mean < 10),
  #                aes(Green_area_mean,X.Area_mean,label=name)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
five

#all greenness with cov and greenness photo about 0.5 r2s with raw data
#PRI with both just 0.1 ns
library(gridExtra)
tiff('Scatterplots_raw_GMGreaanarea_PRIwithFW.tiff', units="in", width=22, height=13, res=300, compression = 'lzw')
grid.arrange(five, four, ncol=3, nrow = 2)
dev.off()
#dev.new()



names(join)

names(dat)

gr <- dat

gr$Dec_LambdaP_mean <- as.numeric(gr$Dec_LambdaP_mean)
gr$Feb2_LambdaP_mean <- as.numeric(gr$Feb2_LambdaP_mean)

library(dplyr)

numeric_cols <- lapply(dat, is.numeric) %>% unlist() %>% unname()
dat2 <- dat[, numeric_cols] #jus tget numeric
numeric_cols #check which numeric
library(corrplot)
#?corrplot
cor(dat2)

#this line not needed for flav data as complete
dat2 <- dat2[complete.cases(dat2), ] #just taken KS06A out

#remove cyanadin-3-gluc
colnames(dat2)
dat2 <- dat2[,-52]

#create matrix

M<-cor(dat2,use='complete.obs')
#M<-cor(dat2)
M[M < 0.5 & M > -0.5] <- 0

#M <- read.csv("Comb_mat_spnew_rmcyandine.csv")
as.matrix(M)

col_classes <- unname(sapply(M, class)) #do I need to do this to define col_classes to then be able to do something else with this? where is 'class' coming from?
Msig <- which(col_classes == 0.5:-0.5)
Msig <- which(col_classes == 0)
replace(Msig, NA)

library(dplyr)
library(corrplot)

#PLOTS NOT WORKING, JUST USING SAVED EXCEL MATRIX

### Computing p-value of correlations. mat : is a matrix of data.
#This is creating a function to compute p-values, you just need to run this.
# ... : further arguments to pass to the native R cor.test function
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)

    n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

# Create a matrix of the p-value of the correlation:
p.mat <- cor.mtest(dat2)

#slightly better plot
plotcolour<- colorRampPalette(c("darkred", "gray88", "darkcyan"))(20) 

corrplot(M, type="upper", order="hclust",method = "circle",
         p.mat = p.mat, sig.level = 0.05, insig = "blank",tl.cex = 0.5, tl.col = 'black',col=plotcolour)

#to run nicer plot
plotcolour<- colorRampPalette(c("darkred", "gray88", "darkcyan"))(20)
#This is where I've created a colour scheme of my choosing. 
#tiff('Comb_corrplot.tiff', units="in", width=14, height=12, res=300, compression = 'lzw')
tiff('Comb_corrplot_flavs+ARI.tiff', units="in", width=14, height=12, res=300, compression = 'lzw')
corrplot(M, col=plotcolour, type="upper", order="hclust", # Add coefficient of correlation
         tl.col="black", tl.cex = 0.4, tl.offset = 0.5,
         #Text label color and rotation
         # Combine with significance
         p.mat = p.mat,
         sig.level = c(.001, .01, .05), pch.cex = 0.7, #add signficance stars
         insig = "label_sig", pch.col = "black",
         # hide correlation coefficient on the principal diagonal
         diag=FALSE )
dev.off()
#tl.cex changes txt size and pch.cex in plot size

#SAVE MATRIX AS EXCEL FILE
write.csv(M, "Comb_mat_spnew_rmcyandine.csv")

#automate linear models
#using dat2 raw for accessions 25 without flav
setwd("C:\\Users\\kelli\\OneDrive - The University of Nottingham\\Nottingham\\Year 2\\UK Lemna - Bradford - Hull\\Flavour analysis\\HS\\lm test")
lms <- expand.grid(1:77, 1:77)
lms_names <- expand.grid(names(dat2)[1:77], names(dat2)[1:77])
out <- vector(mode = "list", length = nrow(lms))
for(i in 1:nrow(lms)){
  lms_col_2 <- lms$Var2[i]
  lms_col_1 <- lms$Var1[i]
  plot_name <- paste0(stringr::str_pad(i, width = 3, pad = "0"), " ", lms_names$Var2[i], " vs. ", lms_names$Var1[i], ".png")
  png(plot_name, width = 500, height = 500, type = "cairo")
  plot(dat2[, lms_col_1], dat2[, lms_col_2], main = paste0(lms_names$Var2[i], " vs. ", lms_names$Var1[i]), type = "n")
  text(dat2[, lms_col_1], dat2[, lms_col_2], row.names(dat2), cex = 0.8)
  abline(lm(dat2[, lms_col_2] ~ dat2[, lms_col_1]))
  #plots per column and output as png with correspinding names
  dev.off()
  out[[i]] <- data.frame(r.squared = summary(lm(dat2[, lms_col_2] ~ dat2[, lms_col_1]))$r.squared,
                         pos_neg = ifelse(summary(lm(dat2[, lms_col_2] ~ dat2[, lms_col_1]))$coef[2, 1] > 0, "+", "-"),
                         p.value = summary(lm(dat2[, lms_col_2] ~ dat2[, lms_col_1]))$coefficients[2, 4])
}
#select r squared and p val from table and make df output
#use estimate col as gradient to make + or - col

#write to files
library(dplyr)
(all_output <- bind_cols(lms_names, do.call(rbind, out)))
all_output <- all_output[all_output$r.squared != 1, ]
all_output <- all_output[all_output$p.value < 0.05, ]
all_output$r.squared <- ifelse(all_output$pos_neg == "+", all_output$r.squared, -all_output$r.squared)

write.csv(all_output, "Lm_correlations.csv")

#pca for corr sig subset
#pca
names(dat)
ions_all <- dat #5 includes species in col 13
# turn accession into factor and classify data to include
sapply(ions_all, class)

#ions_all$EnvLight <- as.factor(ions_all$EnvLight)
options(ggrepel.max.overlaps = Inf)
library(FactoMineR)
library(factoextra)

#plot all vars against each other on pca
#can either remove all observations for KS6A or
#do imputation
ions_all <- read.csv("Comb_hs+photos+biomass+lightlevs_summ_spchange+Ljp_impNAs.csv")
ions_all <- read.csv("Comb_hs+photos+biomass+lightlevs_summ_spchange+Ljp_forreducedpca.csv")

#remove cyanidin column all 0s
ions_acc <- ions_all[ ,c(3:42)] #reduced no light env
ions_acc <- dat[ ,c(3:51,53:81)] #mets
ions_acc <- ions_all[ ,c(3:97)] #flavs 
ions_acc.pca <- PCA(ions_acc, quali.sup=94)
ions_acc.pca <- PCA(ions_acc, quali.sup=78) #cant exceed max col number
ions_acc.pca <- PCA(ions_acc, quali.sup=40) #cant exceed max col number

print(ions_acc.pca)
ions_acc.pca$eig #24 and 18% 1 and 2
#ions_acc.pca$quali.sup
head(ions_acc.pca) #values for each row in pc
print(summary(ions_acc.pca)) #shows pc contirbutions of eigens

#remove columns with 0s in just cyanidin column for met

#do this first to get pc columns, then visualise variable
#ions_acc.pca<- prcomp(ions_all[,c(3:152)], center=T,scale=T) #list can include mix mat and df
ions_acc.pca<- prcomp(ions_acc, center=T,scale=T) #list can include mix mat and df
#tutorial
#tutorial
str(ions_acc.pca)
library(vegan) #used for ellipses

mypc <- ions_acc.pca$x #define new variable x = pcs just need to plot these on xy graph

#FOR ACCESSION POINTS ON LANDSCAPE
png("Comp_flavour_allavgd.png", type="cairo", width = 650, height = 650) #save as png

plot(mypc[,1], mypc[,2], col = reordered_groups,
     las=1, xlab="PC1 (24%)", ylab="PC2 (18%)",
     pch=16, cex=1.5, ylim=c(-12,12), xlim=c(-12,12)) #pcs as columns, produce xy plot, las is rotation of axis numbers, pch plot shape, ylim expand out so legend room
abline(v=0, lty=2, col="lightgrey") #draw line, lty is segmented
abline(h=0, lty=2, col="lightgrey") #0 lines dashed
legend("topright", pch=16, col=reordered_groups, cex=0.6, c("L. minor", "L. japonica", "L. turionifera", "L. minuta", "S. polyrhiza")) #customise legend seperately cex=txtsize cols 1-4 stnd, concat order of places as want to display
text(x=mypc[,1],y=mypc[,2], labels =dat$Accession, pos=2) #txt and pos to define labels on plot, 2 = top of plot, 1 bottom

dev.off() #needs to be off to save it
#dev.new()

#change order of species
dat$Species <- factor(dat$Species, levels = c("L. minor", "L. japonica", "L. turionifera", "L. minuta", "S. polyrhiza"))
reordered_groups <- c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00")

#basic plot
biplot(ions_acc.pca)

fviz_pca_biplot(ions_acc.pca, repel=TRUE, pointsize=6, pointshape=21, col.var="black", arrowsize=0.6, labelsize=5, col.ind=dat$Species, palette=c("green2", "gold", "green2", "gold", "gold"), addEllipses=TRUE, ellipse.type="confidence")

#FOR VARIABLE COS 2 GRAPH
#plot cos 2 as bar graph #high = good representation on pc
fviz_cos2(ions_acc.pca, choice = "var", axes = 1:2, cex.axis=0.05)

fviz_pca_var(ions_acc.pca, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
)

#do as ggplot
pdf('PCA_variables_mets.pdf', width=14, height=12)
plot1 <- fviz_pca_var(ions_acc.pca, col.var = "cos2",
                      gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                      repel = TRUE)
plot1 + xlab("PC1 (36%)") + ylab("PC2 (18%)") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
                     plot.title = element_blank())
dev.off()
#remove title and grids

#plots for seemingly related hs with flavs
library(ggplot2)
one <- ggscatter(dat2, x = "NDWI1_mean", y = "trans.2.heptenal", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "", ylab = "")
#0.51, p =0.009
two <- ggscatter(dat2, x = "NDWI1_mean", y = "Tetradecanal", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "", ylab = "")
#0.51, p =0.009
thr <- ggscatter(dat2, x = "NDWI1_mean", y = "Pentadecanal", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "", ylab = "")
#0.58, p =0.002

fou <- ggscatter(dat2, x = "NDWI_mean", y = "X3.Heptanone", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "", ylab = "")
#-0.55, p =0.004
fou
fiv <- ggscatter(dat2, x = "PRI_mean", y = "Ethanol", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "", ylab = "")
#-0.63, p =0.0006

six <- ggscatter(dat2, x = "Green_area_mean", y = "Ethanol", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "", ylab = "")
#0.52, p =0.007

sev <- ggscatter(dat2, x = "OVI_mean", y = "Ethanol", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "", ylab = "")
#0.52, p =0.008
eig <- ggscatter(dat2, x = "OVI_mean", y = "X2.3.pentanedione", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "", ylab = "")
#0.56, p =0.003

nin <- ggscatter(dat2, x = "NDVI_mean", y = "Benzaldehyde", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "", ylab = "")
#-0.64, p =0.0005
ten <- ggscatter(dat2, x = "GI_mean", y = "Benzaldehyde", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "", ylab = "")
#-0.55, p =0.004

#flav with light levels
twe <- ggscatter(dat2, x = "PRI_mean", y = "Feb2_PPFD_mean", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "", ylab = "")
#-0.5, p =0.001
thi <- ggscatter(dat2, x = "Pyrrole", y = "Feb2_PPFD_mean", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "", ylab = "")
#-0.36, p =0.07
fourt <- ggscatter(dat2, x = "Pyrrole", y = "Feb2_UV_mean", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Pyrrole", ylab = "UV")
#-0.5, p =0.01
fif <- ggscatter(dat2, x = "Decanal", y = "Feb2_G_mean", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Decanal", ylab = "G")
#-0.55, p =0.004
sixt <- ggscatter(dat2, x = "Decanal", y = "Feb2_R_mean", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Decanal", ylab = "R")
#-0.51, p =0.01
sixt
library(gridExtra)
grid.arrange(one, two, thr, fou, fiv, six, sev, eig, nin, ten,
             twe,thi,fourt,fif,sixt)

#Decanal and pyrrole anything to do with light spectrum?
#coriander high red light decreases decanal, no diff in tomato
#uv induces pyrrole, absorbs uv, precusor for chl
pdf('Flavbylight_scatters.pdf', width=14, height=12)
grid.arrange(fourt, fif, sixt, ncol=4, nrow=3)
dev.off()

#now correlate mets with light and hs meas
library(ggplot2)
one <- ggscatter(dat2, x = "GM_mean", y = "Glucose", 
                 add = "reg.line", conf.int = TRUE, 
                 cor.coef = TRUE, cor.method = "pearson",
                 xlab = "", ylab = "")
one
#-0.53, p =0.007
two <- ggscatter(dat2, x = "PRI_mean", y = "FW.g.norm", 
                 add = "reg.line", conf.int = TRUE, 
                 cor.coef = TRUE, cor.method = "pearson",
                 xlab = "", ylab = "")
two
#-0.67, p =0.0002
thr <- ggscatter(dat2, x = "PRI_mean", y = "g.6.months_norm", 
                 add = "reg.line", conf.int = TRUE, 
                 cor.coef = TRUE, cor.method = "pearson",
                 xlab = "", ylab = "")
#-0.67, p =0.0002
thr
fou <- ggscatter(dat2, x = "PRI_mean", y = "Apigenine.7.O.Glc", 
                 add = "reg.line", conf.int = TRUE, 
                 cor.coef = TRUE, cor.method = "pearson",
                 xlab = "", ylab = "")
#-0.56, p =0.004 doesnt look very good
fou
fiv <- ggscatter(dat2, x = "PRI_mean", y = "putative.Chlorogenic.acid.Isomer", 
                 add = "reg.line", conf.int = TRUE, 
                 cor.coef = TRUE, cor.method = "pearson",
                 xlab = "", ylab = "")
#-0.52, p =0.009 doesnt look very good
fiv
six <- ggscatter(dat2, x = "PRI_mean", y = "Cyanidine.3.Glc.1", 
                 add = "reg.line", conf.int = TRUE, 
                 cor.coef = TRUE, cor.method = "pearson",
                 xlab = "", ylab = "")
#-0.51, p =0.01 doesnt look very good
six
sev <- ggscatter(dat2, x = "PRI_mean", y = "Luteoline.7.O.Glc", 
                 add = "reg.line", conf.int = TRUE, 
                 cor.coef = TRUE, cor.method = "pearson",
                 xlab = "", ylab = "")
#-0.51, p =0.01 doesnt look very good
sev
eig <- ggscatter(dat2, x = "PRI_mean", y = "Cyanidin.Mal.Glc", 
                 add = "reg.line", conf.int = TRUE, 
                 cor.coef = TRUE, cor.method = "pearson",
                 xlab = "", ylab = "")
#-0.51, p =0.01 doesnt look very good
eig

#not working well with sec mets as either have a lot or pretty much 0

nin <- ggscatter(dat2, x = "GM_mean", y = "L.Methionine", 
                 add = "reg.line", conf.int = TRUE, 
                 cor.coef = TRUE, cor.method = "pearson",
                 xlab = "", ylab = "")
#0.55, p =0.05
nin
ten <- ggscatter(dat2, x = "GM_mean", y = "L.Histidine", 
                 add = "reg.line", conf.int = TRUE, 
                 cor.coef = TRUE, cor.method = "pearson",
                 xlab = "", ylab = "")
#0.57, p =0.003
ten
elv <- ggscatter(dat2, x = "GM_mean", y = "L.Tyrosine", 
                 add = "reg.line", conf.int = TRUE, 
                 cor.coef = TRUE, cor.method = "pearson",
                 xlab = "", ylab = "")
#0.62, p =0.001 (2)
elv
twe <- ggscatter(dat2, x = "GM_mean", y = "L.Proline", 
                 add = "reg.line", conf.int = TRUE, 
                 cor.coef = TRUE, cor.method = "pearson",
                 xlab = "", ylab = "")
#0.59, p =0.002
twe
thi <- ggscatter(dat2, x = "GM_mean", y = "L.Leucine_2", 
                 add = "reg.line", conf.int = TRUE, 
                 cor.coef = TRUE, cor.method = "pearson",
                 xlab = "", ylab = "")
#0.59, p =0.002
thi
fourt <- ggscatter(dat2, x = "GM_mean", y = "L.Phenylalanine", 
                 add = "reg.line", conf.int = TRUE, 
                 cor.coef = TRUE, cor.method = "pearson",
                 xlab = "", ylab = "")
#0.65, p =0.0006 (1)
fourt
fif <- ggscatter(dat2, x = "GM_mean", y = "L.Isoleucine_1", 
                   add = "reg.line", conf.int = TRUE, 
                   cor.coef = TRUE, cor.method = "pearson",
                   xlab = "", ylab = "")
#0.59, p =0.002
fif
sixt <- ggscatter(dat2, x = "GM_mean", y = "Glycine", 
                 add = "reg.line", conf.int = TRUE, 
                 cor.coef = TRUE, cor.method = "pearson",
                 xlab = "", ylab = "")
#0.58, p =0.003
sixt
sevt <- ggscatter(dat2, x = "Fructose", y = "Glycine", 
                  add = "reg.line", conf.int = TRUE, 
                  cor.coef = TRUE, cor.method = "pearson",
                  xlab = "", ylab = "")
#-0.51, p =0.011
sevt
eigt <- ggscatter(dat2, x = "Glucose", y = "L.Serine", 
                  add = "reg.line", conf.int = TRUE, 
                  cor.coef = TRUE, cor.method = "pearson",
                  xlab = "", ylab = "")
#-0.53, p =0.007
eigt
nint <- ggscatter(dat2, x = "Starch.mg...g", y = "Tyramine", 
                  add = "reg.line", conf.int = TRUE, 
                  cor.coef = TRUE, cor.method = "pearson",
                  xlab = "", ylab = "")
#-0.51, p =0.01
nint
twent <- ggscatter(dat2, x = "Fructose", y = "L.Glutamine", 
                  add = "reg.line", conf.int = TRUE, 
                  cor.coef = TRUE, cor.method = "pearson",
                  xlab = "", ylab = "")
#-0.52, p =0.009
twent
twentone <- ggscatter(dat2, x = "Fructose", y = "L.Asparagine", 
                   add = "reg.line", conf.int = TRUE, 
                   cor.coef = TRUE, cor.method = "pearson",
                   xlab = "", ylab = "")
#-0.51, p =0.01
twentone
twenttwo <- ggscatter(dat2, x = "Fructose", y = "L.Histidine", 
                      add = "reg.line", conf.int = TRUE, 
                      cor.coef = TRUE, cor.method = "pearson",
                      xlab = "", ylab = "")
#-0.53, p =0.007
twenttwo
twentthree <- ggscatter(dat2, x = "Fructose", y = "L.Proline", 
                      add = "reg.line", conf.int = TRUE, 
                      cor.coef = TRUE, cor.method = "pearson",
                      xlab = "", ylab = "")
#-0.55, p =0.005
twentthree
twentfour <- ggscatter(dat2, x = "Fructose", y = "Glycine", 
                      add = "reg.line", conf.int = TRUE, 
                      cor.coef = TRUE, cor.method = "pearson",
                      xlab = "", ylab = "")
#-0.51, p =0.01
twentfour
twentfive <- ggscatter(dat2, x = "Fructose", y = "L.Alanine", 
                      add = "reg.line", conf.int = TRUE, 
                      cor.coef = TRUE, cor.method = "pearson",
                      xlab = "", ylab = "")
#-0.54, p =0.007
twentfive
twentsix <- ggscatter(dat2, x = "Fructose", y = "L.Serine", 
                        add = "reg.line", conf.int = TRUE, 
                        cor.coef = TRUE, cor.method = "pearson",
                        xlab = "", ylab = "")
#-0.62, p =0.001
twentsix
twentsev <- ggscatter(dat2, x = "Fructose", y = "L.Threonine", 
                       add = "reg.line", conf.int = TRUE, 
                       cor.coef = TRUE, cor.method = "pearson",
                       xlab = "", ylab = "")
#-0.54, p =0.006
twentsev
twenteig <- ggscatter(dat2, x = "Fructose", y = "Apigenin.7.O.Glc", 
                      add = "reg.line", conf.int = TRUE, 
                      cor.coef = TRUE, cor.method = "pearson",
                      xlab = "", ylab = "")
#0.66, p =0.0004
twenteig

#put sig correlatiosn into own correlation matrix
sig <- read.csv("Comb_justsig.csv")
sig <- read.csv("Comb_justsig_HS.csv") #just ones correlating with hs
sig <- read.csv("Comb_justsig_HS_remhept.csv") #rem ndwi and heptanone
sig <- read.csv("Comb_justsig_ENV.csv") #correlate light env
library(dplyr)

numeric_cols <- lapply(sig, is.numeric) %>% unlist() %>% unname()
sig2 <- sig[, numeric_cols] #jus tget numeric
numeric_cols #check which numeric
library(corrplot)
#?corrplot
cor(sig2)

#this line not needed for flav data as complete
sig2 <- sig2[complete.cases(sig2), ] #just taken KS06A out

#remove cyanadin-3-gluc
colnames(sig2)
sig2 <- sig2[,-62]
sig2 <- sig2[,-49]

#create matrix

M<-cor(sig2,use='complete.obs')
#M<-cor(dat2)
M[M < 0.5 & M > -0.5] <- 0

#M <- read.csv("Comb_mat_spnew_rmcyandine.csv")
as.matrix(M)

col_classes <- unname(sapply(M, class)) #do I need to do this to define col_classes to then be able to do something else with this? where is 'class' coming from?
Msig <- which(col_classes == 0.5:-0.5)
Msig <- which(col_classes == 0)
replace(Msig, NA)

library(dplyr)
library(corrplot)

#PLOTS NOT WORKING, JUST USING SAVED EXCEL MATRIX

### Computing p-value of correlations. mat : is a matrix of data.
#This is creating a function to compute p-values, you just need to run this.
# ... : further arguments to pass to the native R cor.test function
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

# Create a matrix of the p-value of the correlation:
p.mat <- cor.mtest(sig2)

#slightly better plot
plotcolour<- colorRampPalette(c("darkred", "gray88", "darkcyan"))(20) 

corrplot(M, type="upper", order="hclust",method = "circle",
         p.mat = p.mat, sig.level = 0.05, insig = "blank",tl.cex = 0.5, tl.col = 'black',col=plotcolour)

#to run nicer plot
plotcolour<- colorRampPalette(c("darkred", "gray88", "darkcyan"))(20)
#This is where I've created a colour scheme of my choosing. 
pdf('Comb_corrplot_SIG_HS_remhept.pdf', width=14, height=12)
tiff('Comb_corrplot_SIG_HS_remhept.tiff', units="in", width=14, height=12, res=300, compression = 'lzw')
corrplot(M, col=plotcolour, type="upper", order="hclust", # Add coefficient of correlation
         tl.col="black", tl.cex = 1, tl.offset = 0.5,
         #Text label color and rotation
         # Combine with significance
         p.mat = p.mat,
         sig.level = c(.001, .01, .05), pch.cex = 0.7, #add signficance stars
         insig = "label_sig", pch.col = "black",
         # hide correlation coefficient on the principal diagonal
         diag=FALSE )
dev.off()

#looked sig before but now not? ndwi1 with pentadecanal, tetradecanal, trans-2-heptanal
#PRI and ethanol
one <- ggscatter(sig2, x = "PRI_mean", y = "Ethanol", 
                 add = "reg.line", conf.int = TRUE, 
                 cor.coef = TRUE, cor.method = "pearson",
                 xlab = "", ylab = "")
one
#0.51, p =0.009

#include environmental data and see if effect measurements with stats/boxplots
env <- read.csv("Comb_hs+photos+biomass+lightlevs_summ_spchange+Ljp_impNAs_renamed+flavours+env.csv")

#loop kruskal wallis for columns
results <- list()
for(i in names(env[,2:143])){  
  results[[i]] <- kruskal.test(formula(paste(i, "~ Light.habitat")), data = env)
}

r1 <- unlist(results)
r1 <- as.data.frame(r1)
write.csv(r1, "KWtest_Ljp_all_envlight.csv")

#manually check file for sig compounds
kruskal.test(NDVI_mean ~ Light.habitat, data = env)
#same as in excel print out so look for p val < 0.05

#sig variables
#Water.pH
#Env_Light
#Longitude
#Benzene..1.3.bis.1.1.dimethylethyl

#higher ph, light, benzene in dhl
#higher longitude in dll

kruskal.test(NDVI_mean ~ Light.habitat, data = env)

#boxplots for sig
boxplot(Water.pH~Light.habitat,data=env, col = c("red", "blue"),
        ylab = "Water.pH", xlab =  "", las=2,
        ylim=c(5,10),
        cex.lab=1.5, cex.axis=1.5)#** 0.02

boxplot(Env_Light~Light.habitat,data=env, col = c("red", "blue"),
        ylab = "Env_light", xlab =  "", las=2,
        ylim=c(0,2000),
        cex.lab=1.5, cex.axis=1.5)#*** 0.0002

boxplot(Longitude~Light.habitat,data=env, col = c("red", "blue"),
        ylab = "Longitude", xlab =  "", las=2,
        ylim=c(-6,2),
        cex.lab=1.5, cex.axis=1.5)#** 0.02

boxplot(Benzene..1.3.bis.1.1.dimethylethyl..~Light.habitat,data=env, col = c("red", "blue"),
        ylab = "Benzene compound", xlab =  "", las=2,
        ylim=c(0,22),
        cex.lab=1.5, cex.axis=1.5)#* 0.04

#loop kruskal wallis for columns
results <- list()
for(i in names(env[,2:144])){  
  results[[i]] <- kruskal.test(formula(paste(i, "~ Location")), data = env)
}

r1 <- unlist(results)
r1 <- as.data.frame(r1)
write.csv(r1, "KWtest_Ljp_all_loco.csv")

#sites in the south had more light, higher a, w temps, benzaldehyde
#L. japonica from sites in the north had more 1-penten-3-one
#L. japonica and S. polyrhiza from south had more benzaldehyde
#inbalanced comparison only 1 Ljp south, only 1 Sp north

boxplot(X1.penten.3.one~Location,data=env, col = c("red", "blue"),
        ylab = "X1.penten.3.one", xlab =  "", las=2,
        ylim=c(0,500),
        cex.lab=1.5, cex.axis=1.5)#* 0.03

boxplot(Benzaldehyde~Location,data=env, col = c("red", "blue"),
        ylab = "Benzaldehyde", xlab =  "", las=2,
        ylim=c(0,50),
        cex.lab=1.5, cex.axis=1.5)#* 0.04

#check with other tests as KW supposed to be for 3 variables
pairwise.wilcox.test(env$X1.penten.3.one, env$Location,
                     p.adjust.method = "BH")
pairwise.wilcox.test(env$Benzaldehyde, env$Location,
                     p.adjust.method = "BH")
pairwise.wilcox.test(env$Benzene..1.3.bis.1.1.dimethylethyl.., env$Light.habitat,
                     p.adjust.method = "BH")
