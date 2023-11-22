setwd("C:\\Users\\kelli\\OneDrive - The University of Nottingham\\Nottingham\\Year 2\\UK Lemna - Bradford - Hull\\Flavour analysis\\HS")

hs <- read.csv("Glasshouse_temp.csv")

#install.packages("xts")
library(xts)

#avgs per day
day <- aggregate(Temperature ~ Date, hs, mean)

sample.xts <- as.xts(hs)

#lauras
##TINYTAG DATA ANALYSIS

#LOAD PACKAGES 
library(tidyverse) #Comprises ggplot2 dplyr, tidyr. Contains the pipe %>%
library(car) #levene's test
library(rcompanion) #Normality histogram with overlaid line
library(ggstatsplot) #Used for the ggbetweenstats - allows labelling of boxplots for outlier detection
library(gridExtra) #Making multiplots
library(ggsignif) #Add P values or sig. stars to your plots
library(corrplot) #correlation matrix 
library(Hmisc) #correlation matrix
library(ggpubr) #correlations pearsons and spearman



attach(hs)

data <- hs

#SPLIT BY TOD-Time of day
Day <- filter(data, TOD == "Day")
Night <- filter(data, TOD == "Night")

#ADD IN - WORK OUT DAT
#AVERAGES FOR DAY AND NIGHT
summaryT <- group_by(Day, DAT)
summaryT <- summarise(summaryT,n = n(), avg = mean(T), sd = sd(T), se=sd(T)/sqrt(n()))
#write.csv(summaryT, "SummaryT.csv")

summaryT2 <- group_by(Night, DAT)
summaryT2 <- summarise(summaryT2,n = n(), avg = mean(T), sd = sd(T), se=sd(T)/sqrt(n()))
#write.csv(summaryT2, "SummaryT2.csv")

##LINE GRAPHS
##DAYTIME TEMP
p1 <- ggplot(data = summaryT, aes(x = DAT, y = avg,group=1))+
  geom_path()+
  geom_point(color = "blue")+
  geom_line(color = "blue")+
  geom_hline (yintercept = mean(Day$T), linetype = 2, color= "blue")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+ #Remove the standard background and grid lines by setting to element_blank
  theme(axis.line.x = element_line(color = "black", size = 0.5),axis.line.y = element_line(color = "black", size = 0.5))+ 
  theme(axis.title.y = element_text(size = 14, angle = 90, vjust = 6, hjust = 0.5))+ theme(axis.title.x = element_text(size = 14, angle = 0, vjust = -3, hjust = 0.5))+
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  theme(axis.text.x = element_text(color = "black", size = 12, angle = 0), axis.text.y = element_text(color = "black",size = 12, angle = 0))+ #Define format for x and y axis
  theme(legend.position = c(0.5,0.9))+ # Select coordinates to position the legend
  theme(legend.direction = "horizontal")+ #Set the legend direction to horizontal 
  theme(legend.text = element_text(size = 12))+ #Edit the text size of the legend
  theme(panel.border = element_rect(color = "black",fill = NA, size = 0.5))+ #Add a border around the graph
  xlab ("Experimental days")+ 
  ylab(expression("Daytime temperature ("*~degree*C*")")) + #Add a title with subscript font
  labs(tag = substitute(paste(bold('A'))))+ #Label the graph alphabetically in bold
  geom_errorbar(aes(ymin = avg-se, ymax = avg + se), width=.2, position = position_dodge(.9))+ #errorbars
  scale_y_continuous(expand = c(0,0), breaks = seq(15, 25, by = 3), limits = c(15,25))+
  scale_x_continuous(expand = c(0,0), breaks = seq(0, 120, by = 10), limits = c(0,120))

#NIGHT-TIME TEMP
p2 <- ggplot(data = summaryT2, aes(x = DAT, y = avg,group = 1))+
  geom_path()+
  geom_point()+
  geom_hline(yintercept = mean(Night$T), linetype = 2)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+ #Remove the standard background and grid lines by setting to element_blank
  theme(axis.line.x = element_line(color = "black", size = 0.5),axis.line.y = element_line(color = "black", size = 0.5))+ 
  theme(axis.title.y = element_text(size = 14, angle = 90, vjust = 6, hjust = 0.5))+ theme(axis.title.x = element_text(size = 14, angle = 0, vjust = -3, hjust = 0.5))+
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  theme(axis.text.x = element_text(color = "black", size = 12, angle = 0), axis.text.y = element_text(color = "black",size = 12, angle = 0))+ #Define format for x and y axis
  theme(legend.position = c(0.5,0.9))+ # Select coordinates to position the legend
  theme(legend.direction = "horizontal")+ #Set the legend direction to horizontal 
  theme(legend.text = element_text(size = 12))+ #Edit the text size of the legend
  theme(panel.border = element_rect(color = "black",fill = NA, size = 0.5))+ #Add a border around the graph
  xlab ("Experimental days")+ 
  ylab(expression("Night-time temperature ("*~degree*C*")")) + #Add a title with subscript font
  labs(tag = substitute(paste(bold('B'))))+ #Label the graph alphabetically in bold
  geom_errorbar(aes(ymin = avg-se, ymax = avg + se), width=.2, position = position_dodge(.9))+ #errorbars
  scale_y_continuous(expand = c(0,0), breaks = seq(15, 25, by = 3), limits = c(15,25))+
  scale_x_continuous(expand = c(0,0), breaks = seq(0, 120, by = 10), limits = c(0,120))

grid.arrange(p1, p2, ncol = 2)
