library(ggplot2)
library(ncar)
library(tidyverse)
library(DescTools)
library(ggridges)

#clean environment of old items
rm(list = ls())

#identify all files in result folder
Result_Folder <- "D:/Martin_KI/RunningProjects/CASK/ImmunoFluorescence/Puncta_quantification/Homer_Syn_Quantification/ResultFiles/" #Change this path to the folder with your result .csv files from ImageJ macro
Files <- list.files(Result_Folder, pattern = "\\.csv$")

#read sample to cell line code
sample_code <- read.table("D:/Martin_KI/RunningProjects/CASK/ImmunoFluorescence/Puncta_quantification/Puncta_quantification/Sample_code.txt", sep = "\t", header = TRUE)

#Set up empty data frame to contain data
summary_Homer <- data.frame()
summary_Syn <- data.frame()
summary_Nuclei <- data.frame()

#for loop to populate data frame with individual puncti from each image one-by-one
for (i in 1:length(Files)){
  this_file <- Files[i]
  print(this_file)
  this_cell <- unlist(strsplit(this_file,"_"))
  this_sample <- this_cell[1]
  this_cellLine <- sample_code[as.integer(this_sample),2]
  this_time <- sample_code[as.integer(this_sample),4]
  this_replicate <- this_cell[3]
  Punctae <- read.table(paste("D:/Martin_KI/RunningProjects/CASK/ImmunoFluorescence/Puncta_quantification/Homer_Syn_Quantification/ResultFiles/", this_file, sep = ""), sep = ",", header = TRUE)
  this_cell_puncta <- data.frame(CellLine = this_cellLine, Time = this_time, Sample = this_sample, Replicate = this_replicate, puncti = Punctae$Area, Circularity = Punctae$Circ., IntDen = Punctae$IntDen)
  if (startsWith(this_cell[length(this_cell)], "H")){
    summary_Homer<- rbind(summary_Homer,this_cell_puncta)
  }
  else if (startsWith(this_cell[length(this_cell)], "S")){
    summary_Syn<- rbind(summary_Syn,this_cell_puncta)
  }
  else if (startsWith(this_cell[length(this_cell)], "N")){
    summary_Nuclei<- rbind(summary_Nuclei,this_cell_puncta)
  }
  else {
    next()
  }
}

#Set boundary (area) for nuclei
lower_Nuclei <- 1000
summary_Nuclei_interest <- summary_Nuclei[(summary_Nuclei$puncti > lower_Nuclei),]

summary_Nuclei_interest$CellLine <- factor(summary_Nuclei_interest$CellLine, c("CTRL9II", "AF22", "ASD17AII", "AFA0I"))
p <- ggplot(data = summary_Nuclei_interest, aes(x = CellLine, y = puncti, color = Replicate, fill = Sample)) +
  geom_boxplot() +
  theme_classic() +
  ggtitle("Nuclei_Size")#+
  #ylim(0,40000) #+
#coord_cartesian(ylim = c(0,10))
p


p <- ggplot(summary_Nuclei_interest, aes(x = CellLine, y = puncti))+
  geom_boxplot() +
  theme_classic() +
  coord_cartesian(ylim = c(0,8000))
p

#compute summary statistics for Nuclei
summary_stats_Nuclei <- summary_Nuclei_interest %>%
  group_by(Replicate, Sample, CellLine) %>%
  summarise(average = mean(puncti), medianNuc = median(puncti), sd = sd(puncti), count = length(puncti), sum_of_area = sum(puncti))
summary_stats_Nuclei

summary_stats_Nuclei$CellLine <- factor(summary_stats_Nuclei$CellLine, c("CTRL9II", "AF22", "ASD17AII", "AFA0I"))
p <- ggplot(summary_stats_Nuclei, aes(x = CellLine, y = average))+
  geom_boxplot() +
  theme_classic() +
  ylim(0,9000) +
  ggtitle("NucleiSize")
p

p <-  ggplot(data = summary_stats_Nuclei, aes(x = CellLine, y = count, color = CellLine)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(shape= 16, position = position_jitter(0.2), size = 2)+
  theme_classic() +
  ggtitle("Nuclei_Count")#+
#ylim(0,20)
p


summary_stats_Nuclei_NoReplicate <- summary_Nuclei_interest %>%
  group_by(Sample, CellLine) %>%
  summarise(average = mean(puncti), medianNuc = median(puncti), sd = sd(puncti), count = length(puncti), sum_of_area = sum(puncti))
res.Nucleiaov <- aov(average ~ CellLine, data = summary_stats_Nuclei_NoReplicate)
summary(res.Nucleiaov)
TukeyHSD(res.Nucleiaov)

#identify all files in result folder
Synapse_Folder <- "D:/Martin_KI/RunningProjects/CASK/ImmunoFluorescence/Puncta_quantification/Homer_Syn_Quantification/ResultFiles/Synapse" #Change this path to the folder with your result .csv files from ImageJ macro
Synapse_Files <- list.files(Synapse_Folder, pattern = "SynapseCounterResults\\.csv$")

#Set up empty data frame to contain data
summary_Synapses <- data.frame()

#for loop to populate data frame with individual puncti from each image one-by-one
for (i in 1:length(Synapse_Files)){
  this_file <- Synapse_Files[i]
  print(this_file)
  this_cell <- unlist(strsplit(this_file,"_"))
  this_sample <- this_cell[1]
  this_cellLine <- sample_code[as.integer(this_sample),2]
  this_time <- sample_code[as.integer(this_sample),4]
  this_replicate <- this_cell[3]
  Synapse <- read.table(paste("D:/Martin_KI/RunningProjects/CASK/ImmunoFluorescence/Puncta_quantification/Homer_Syn_Quantification/ResultFiles/Synapse/", this_file, sep = ""), sep = ",", header = TRUE)
  this_cell_puncta <- data.frame(CellLine = this_cellLine, Time = this_time, Sample = this_sample, Replicate = this_replicate, PreSynapseNo = Synapse$Presyn..N, PreSynapseSize = Synapse$Presyn..mean.size, PostSynapseNo = Synapse$Postsyn..N, PostSynapseSize = Synapse$Postsyn..mean.size, ColoNo = Synapse$Coloc..N, ColoSize = Synapse$Coloc..mean.size)
  summary_Synapses<- rbind(summary_Synapses,this_cell_puncta)
}

summary_Synapses <- merge(summary_Synapses, summary_stats_Nuclei, by = c("Sample", "Replicate"))
summary_Synapses <- summary_Synapses %>% 
  rename(
    CellLine = CellLine.x
  )

summary_Synapses <- summary_Synapses %>%
  mutate(PostSynapsePerNuclei = PostSynapseNo/count)
summary_Synapses <- summary_Synapses %>%
  mutate(PreSynapsePerNuclei = PreSynapseNo/count)
summary_Synapses <- summary_Synapses %>%
  mutate(ColocSynapse = ColoNo/count)

#this part needs to be changed dynamically to test:
#PreSynapseNo, PreSynapseSize, PostSynapseNo, PostSynapseSize, ColoNo or ColoSize
#PostSynapsePerNuclei, PreSynapsePerNuclei or ColocSynapsePerNuclei
testthis <- "PostSynapsePerNuclei"

summary_Synapses$Status <- ifelse(summary_Synapses$CellLine == "ASD17AII", "ASD", ifelse(summary_Synapses$CellLine == "AFA0I", "MICPCH", "Control"))

res.Synapseaov <- aov(get(testthis) ~ CellLine, data = summary_Synapses)
TukeyHSD(res.Synapseaov)

summary_stats_Synapses <- summary_Synapses %>%
  group_by(Sample, Status) %>%
  summarise(average = mean(get(testthis)), sd = sd(get(testthis)))
#summary_stats_Synapses
res.Synapseaov <- aov(average ~ Status, data = summary_stats_Synapses)
TukeyHSD(res.Synapseaov)

summary_stats_Synapses <- summary_Synapses %>%
  group_by(Sample, CellLine) %>%
  summarise(average = mean(get(testthis)), sd = sd(get(testthis)))
#summary_stats_Synapses
res.Synapseaov <- aov(average ~ CellLine, data = summary_stats_Synapses)
reportStatistics <- TukeyHSD(res.Synapseaov)
reportStatistics <- as.data.frame(reportStatistics[1:1])

#write.csv(reportStatistics, paste("Homer1_Syn_", testthis, ".csv", sep=""))

summary_stats_Synapses$CellLine <- factor(summary_stats_Synapses$CellLine, c("CTRL9II", "AF22", "ASD17AII", "AFA0I"))
p <- ggplot(data = summary_stats_Synapses, aes(x = CellLine, y = average, color =  CellLine)) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic() +
  geom_jitter()+
  ggtitle(testthis)
p

#ggsave(paste(testthis,"Homer_Syn.svg", sep = ""), height = 5, width = 5, dpi = 600)
#ggsave(paste(testthis,"Homer_Syn.svg", sep = ""), height = 5, width = 5, dpi = 600)
#ggsave(paste(testthis,"Homer_Syn.svg", sep = ""), height = 5, width = 5, dpi = 600)
