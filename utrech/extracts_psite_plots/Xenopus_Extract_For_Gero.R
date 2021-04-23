#install packages

# install.packages("tidyverse")
# install.packages("gplots")
# install.packages("devtools")
library(devtools)
# install_github("jokergoo/ComplexHeatmap")
# install.packages("heatmap3")
devtools::install_github("yutannihilation/gghighlight")

#Load packages#

library(ggpubr)
library(tidyr)
library(gplots)
library(dplyr)
library(tibble)
library(ggplot2)
library(ComplexHeatmap)
library(heatmap3)
library(reshape2)
library(gghighlight)
library(ggsci)

#Load data (beware that had to remove first useless three rows, otherwise everything was loaded as character)

getwd()

pSites_ANOVASig_Zscored <- read.csv('Phosphosites_XL_Extract_ANOVA0.05_MedianNorm_Avgs_Zsocred_CorrectGroups.csv', stringsAsFactors = F)
str(pSites_ANOVASig_Zscored)
summary(pSites_ANOVASig_Zscored)
head(pSites_ANOVASig_Zscored)

#Creating unique useful identifier

pSites_ANOVASig_Zscored_Unique_ID <- unite(pSites_ANOVASig_Zscored, 
                                           'Protein_Amino.acid_Site_NumberPhospho', 
                                           c('Protein', 'Amino.acid', 'Position', 'Multiplicity'),
                                           sep ='|', remove = FALSE)


#Df with only intensities of time point of interest
#Individual (per bio replicate)

pSites_ANOVASig_Zscored_Unique_ID_BioRep <-
  pSites_ANOVASig_Zscored_Unique_ID[, c(1:21)]

#Means

pSites_ANOVASig_Zscored_Unique_ID_Means <-
  pSites_ANOVASig_Zscored_Unique_ID[, c(21, 47:53)]

#Medians

pSites_ANOVASig_Zscored_Unique_ID_Medians <-
  pSites_ANOVASig_Zscored_Unique_ID[, c(21, 55:61)]

#Preparing data
#reordering columns

col_order_BioRep <- c('Interphase_000_A', 'Interphase_000_B', 'Interphase_000_C', 'Interphase_030_A', 'Interphase_030_B_R', 
                      'Interphase_030_C', 'Interphase_060_A', 'Interphase_060_B', 'Interphase_060_C', 'Interphase_090_A', 
                      'Interphase_090_B', 'Interphase_090_C', 'Interphase_120_A_R', 'Interphase_120_B_R', 'Interphase_120_C', 
                      'Mitosis_CyB_A', 'Mitosis_CyB_C_R', 'Mitosis_CyB_D', 'Mitosis_CSF', 'Mitosis_CSF_TR_R','Protein_Amino.acid_Site_NumberPhospho')
col_order_Means <- c('Mean.Interphase_000', "Mean.Interphase_030", 'Mean.Interphase_060', 'Mean.Interphase_090',
                     'Mean.Interphase_120', 'Mean.Mitosis_CyB','Mean.Mitosis_CSF',
                     'Protein_Amino.acid_Site_NumberPhospho')
col_order_Medians <- c('Median.Interphase_000', "Median.Interphase_030", 'Median.Interphase_060', 'Median.Interphase_090',
                       'Median.Interphase_120','Median.Mitosis_CyB', 'Median.Mitosis_CSF',
                       'Protein_Amino.acid_Site_NumberPhospho')

pSites_ANOVASig_Zscored_Unique_ID_BioRep_Ord <- pSites_ANOVASig_Zscored_Unique_ID_BioRep[, col_order_BioRep]

pSites_ANOVASig_Zscored_Unique_ID_Means_Ord <- pSites_ANOVASig_Zscored_Unique_ID_Means[, col_order_Means]

pSites_ANOVASig_Zscored_Unique_ID_Medians_Ord <- pSites_ANOVASig_Zscored_Unique_ID_Medians[, col_order_Medians]

#Transforming to numeric matrix

pSites_ANOVASig_Zscored_Unique_ID_BioRep_Ord_Matrix <-
  as.matrix(sapply(pSites_ANOVASig_Zscored_Unique_ID_BioRep_Ord, as.numeric))

pSites_ANOVASig_Zscored_Unique_ID_Means_Ord_Matrix <-
  as.matrix(sapply(pSites_ANOVASig_Zscored_Unique_ID_Means_Ord, as.numeric))

pSites_ANOVASig_Zscored_Unique_ID_Medians_Ord_Matrix <-
  as.matrix(sapply(pSites_ANOVASig_Zscored_Unique_ID_Medians_Ord, as.numeric))

#Adding rownames

rownames(pSites_ANOVASig_Zscored_Unique_ID_BioRep_Ord_Matrix) <- pSites_ANOVASig_Zscored_Unique_ID_BioRep_Ord[, 21]
rownames(pSites_ANOVASig_Zscored_Unique_ID_Means_Ord_Matrix) <- pSites_ANOVASig_Zscored_Unique_ID_Means_Ord[, 8]
rownames(pSites_ANOVASig_Zscored_Unique_ID_Medians_Ord_Matrix) <- pSites_ANOVASig_Zscored_Unique_ID_Medians_Ord[, 8]

#removing useless column

pSites_ANOVASig_Zscored_Unique_ID_BioRep_Ord_Matrix_ready <- pSites_ANOVASig_Zscored_Unique_ID_BioRep_Ord_Matrix[, -c(21)]
pSites_ANOVASig_Zscored_Unique_ID_Means_Ord_Matrix_ready <- pSites_ANOVASig_Zscored_Unique_ID_Means_Ord_Matrix[, -c(8)]
pSites_ANOVASig_Zscored_Unique_ID_Medians_Ord_Matrix_ready <- pSites_ANOVASig_Zscored_Unique_ID_Medians_Ord_Matrix[, -c(8)]

#Plotting with Complex Heatmap
library(circlize)
col_fun = colorRamp2(c(-2.5, 0, 2.5), c("cyan", "black", "yellow"))
summary(pSites_ANOVASig_Zscored_Unique_ID_Means_Ord_Matrix_ready)

set.seed(5346)
HM = Heatmap(pSites_ANOVASig_Zscored_Unique_ID_Means_Ord_Matrix_ready, 
             column_order = colnames(pSites_ANOVASig_Zscored_Unique_ID_Means_Ord_Matrix_ready),
             show_row_names = FALSE, row_km = 6, row_km_repeats = 100, clustering_distance_rows = 'spearman', 
             clustering_method_rows = 'average', name = 'Z-Score', col = col_fun, row_dend_width = unit(3, 'cm'))


#Extracting row clusters
#From https://github.com/jokergoo/ComplexHeatmap/issues/136

HM = draw(HM)           #assign Heat Map to object and draw it
r.dend <- row_dend(HM)  #Extract row dendrogram
rcl.list <- row_order(HM)  #Extract clusters (output is a list)
lapply(rcl.list, function(x) length(x))  #check/confirm size clusters

# loop to extract genes for each cluster.
for (i in 1:length(row_order(HM))){
  if (i == 1) {
    clu <- t(t(row.names(pSites_ANOVASig_Zscored_Unique_ID_Means_Ord_Matrix_ready[row_order(HM)[[i]],])))
    out <- cbind(clu, paste("cluster", i, sep=""))
    colnames(out) <- c("ID", "Cluster")
  } else {
    clu <- t(t(row.names(pSites_ANOVASig_Zscored_Unique_ID_Means_Ord_Matrix_ready[row_order(HM)[[i]],])))
    clu <- cbind(clu, paste("cluster", i, sep=""))
    out <- rbind(out, clu)
  }
}

out

#Matching cluster to identifier in scaled matrix, in order to plot heatmap for different clusters
#number of clusters corresponds to the number from top to bottom of the heatmap (logical order) not the onw shown on the heatmap
#transforming matrices to data frames in order to use merge

pSites_ANOVASig_Zscored_Unique_ID_Means_Ord_Matrix_ready_df <-
  as.data.frame(pSites_ANOVASig_Zscored_Unique_ID_Means_Ord_Matrix_ready)
out_df <- as.data.frame(out)


#rownames to column to use merge

pSites_ANOVASig_Zscored_Unique_ID_Means_Ord_Matrix_ready_df_NoRowNames <- 
  pSites_ANOVASig_Zscored_Unique_ID_Means_Ord_Matrix_ready_df %>% rownames_to_column('ID')

#merging

pSites_with_clusters_fromHM_Scaled <- merge(x = pSites_ANOVASig_Zscored_Unique_ID_Means_Ord_Matrix_ready_df_NoRowNames,
                                            y = out_df, by.x = 'ID', by.y = 'ID')


#Filtering based on cluster to get data frame for each cluster

Cluster1 <- pSites_with_clusters_fromHM_Scaled[pSites_with_clusters_fromHM_Scaled$Cluster == 'cluster1', ]

Cluster2 <- pSites_with_clusters_fromHM_Scaled[pSites_with_clusters_fromHM_Scaled$Cluster == 'cluster2', ]

Cluster3 <- pSites_with_clusters_fromHM_Scaled[pSites_with_clusters_fromHM_Scaled$Cluster == 'cluster3', ]

Cluster4 <- pSites_with_clusters_fromHM_Scaled[pSites_with_clusters_fromHM_Scaled$Cluster == 'cluster4', ]

Cluster5 <- pSites_with_clusters_fromHM_Scaled[pSites_with_clusters_fromHM_Scaled$Cluster == 'cluster5', ]

Cluster6 <- pSites_with_clusters_fromHM_Scaled[pSites_with_clusters_fromHM_Scaled$Cluster == 'cluster6', ]

#writing tables

write.table(Cluster1, file = "Cluster1.txt", sep="\t", quote=F, row.names=FALSE)
write.table(Cluster2, file = "Cluster2.txt", sep="\t", quote=F, row.names=FALSE)
write.table(Cluster3, file = "Cluster3.txt", sep="\t", quote=F, row.names=FALSE)
write.table(Cluster4, file = "Cluster4.txt", sep="\t", quote=F, row.names=FALSE)
write.table(Cluster5, file = "Cluster5.txt", sep="\t", quote=F, row.names=FALSE)
write.table(Cluster6, file = "Cluster6.txt", sep="\t", quote=F, row.names=FALSE)

#transforming to matrix

Cluster1_matrix <- data.matrix(Cluster1)
rownames(Cluster1_matrix) <- Cluster1[, 1]
Cluster1_matrix <- Cluster1_matrix[, -c(1, 9)]

Cluster2_matrix <- data.matrix(Cluster2)
rownames(Cluster2_matrix) <- Cluster2[, 1]
Cluster2_matrix <- Cluster2_matrix[, -c(1, 9)]

Cluster3_matrix <- data.matrix(Cluster3)
rownames(Cluster3_matrix) <- Cluster3[, 1]
Cluster3_matrix <- Cluster3_matrix[, -c(1, 9)]

Cluster4_matrix <- data.matrix(Cluster4)
rownames(Cluster4_matrix) <- Cluster4_matrix[, 1]
Cluster4_matrix <- Cluster4_matrix[, -c(1, 9)]

Cluster5_matrix <- data.matrix(Cluster5)
rownames(Cluster5_matrix) <- Cluster5_matrix[, 1]
Cluster5_matrix <- Cluster5_matrix[, -c(1, 9)]

Cluster6_matrix <- data.matrix(Cluster6)
rownames(Cluster6_matrix) <- Cluster6_matrix[, 1]
Cluster6_matrix <- Cluster6_matrix[, -c(1, 9)]



#Loading data before Z-scoring in order to plot trends of normalized intensities per bio replicate with ggplot

pSites_ANOVASig_Normalized <- read.csv('Phosphosites_XL_Extract_ANOVASig_Norm_NoOA.csv', stringsAsFactors = F)

#Creating unique useful identifier

pSites_ANOVASig_Normalized_Unique_ID <- unite(pSites_ANOVASig_Normalized, 
                                              'Protein_Amino.acid_Site_NumberPhospho', 
                                              c('Protein', 'Amino.acid', 'Position', 'Multiplicity'),
                                              sep ='|', remove = FALSE)


#Df with only columns of interest


pSites_ANOVASig_Normalized_Unique_ID_BioRep <-
  pSites_ANOVASig_Normalized_Unique_ID[, c(1:22, 26, 38, 40:42)]

#Creating object to organize columnes

col_order_BioRep_Normalized <- c('Interphase_000_A', 'Interphase_000_B', 'Interphase_000_C', 'Interphase_030_A', 'Interphase_030_B_R', 
                                 'Interphase_030_C', 'Interphase_060_A', 'Interphase_060_B', 'Interphase_060_C', 'Interphase_090_A', 
                                 'Interphase_090_B', 'Interphase_090_C', 'Interphase_120_A_R', 'Interphase_120_B_R', 'Interphase_120_C', 
                                 'Mitosis_CyB_A', 'Mitosis_CyB_C_R', 'Mitosis_CyB_D', 'Mitosis_CSF', 'Mitosis_CSF_TR_R','Protein_Amino.acid_Site_NumberPhospho',
                                 'Amino.acid', 'Multiplicity', 'Proteins', 'Leading.proteins', 'Protein', 'Sequence.window')

#reordering columns 

pSites_ANOVASig_Normalized_Unique_ID_BioRep_Ord <- pSites_ANOVASig_Normalized_Unique_ID_BioRep[, col_order_BioRep_Normalized]

#merging

pSites_with_clusters_fromHM_Normalized <- merge(x = pSites_ANOVASig_Normalized_Unique_ID_BioRep_Ord,
                                                y = out_df, by.x = 'Protein_Amino.acid_Site_NumberPhospho', by.y = 'ID')

#Tidying up the data

#Changing colnames 

pSites_with_clusters_fromHM_Normalized_tidy <- plyr::rename(pSites_with_clusters_fromHM_Normalized, 
                                                            c('Interphase_120_A_R' = 'Interphase_120_A', 'Interphase_120_B_R' = 'Interphase_120_B',
                                                              'Interphase_030_B_R' = 'Interphase_030_B', 'Mitosis_CyB_C_R' = 'Mitosis_CyB_C',
                                                              'Mitosis_CSF_TR_R' = 'Mitosis_CSF_B', 'Mitosis_CSF' = 'Mitosis_CSF_A', 'Protein_Amino.acid_Site_NumberPhospho' = 'ID'))
#checking properties of the data frame

summary(pSites_with_clusters_fromHM_Normalized_tidy)
str(pSites_with_clusters_fromHM_Normalized_tidy)

#changing several columns to factors

pSites_with_clusters_fromHM_Normalized_tidy$Amino.acid <- as.factor(pSites_with_clusters_fromHM_Normalized_tidy[, 'Amino.acid'])
pSites_with_clusters_fromHM_Normalized_tidy$Multiplicity <- as.factor(pSites_with_clusters_fromHM_Normalized_tidy[, 'Multiplicity'])
pSites_with_clusters_fromHM_Normalized_tidy$Proteins <- as.factor(pSites_with_clusters_fromHM_Normalized_tidy[, 'Proteins'])
pSites_with_clusters_fromHM_Normalized_tidy$Leading.proteins <- as.factor(pSites_with_clusters_fromHM_Normalized_tidy[, 'Leading.proteins'])
pSites_with_clusters_fromHM_Normalized_tidy$Protein <- as.factor(pSites_with_clusters_fromHM_Normalized_tidy[, 'Protein'])
pSites_with_clusters_fromHM_Normalized_tidy$Cluster <- as.factor(pSites_with_clusters_fromHM_Normalized_tidy[, 'Cluster'])

#Writing table to get sequences and clusters for motif analysis

write.table(pSites_with_clusters_fromHM_Normalized_tidy,
            file = "pSites_with_clusters_fromHM_Normalized_tidy.txt", sep="\t", 
            quote=F, row.names=FALSE)


#From wide format to long format

pSites_with_clusters_fromHM_Normalized_tidy_v1 <- melt(pSites_with_clusters_fromHM_Normalized_tidy)

#Separate the variable into 'Condition', 'Time/Treatment', 'Replicate'

pSites_with_clusters_fromHM_Normalized_tidy_v2 <- pSites_with_clusters_fromHM_Normalized_tidy_v1 %>%
  separate(variable, c('Condition', 'Time_Or_Treatment', 'Replicate'), '_')

str(pSites_with_clusters_fromHM_Normalized_tidy_v2)
summary(pSites_with_clusters_fromHM_Normalized_tidy_v2)

#changing several columns to factors

pSites_with_clusters_fromHM_Normalized_tidy_v2$Condition <- as.factor(pSites_with_clusters_fromHM_Normalized_tidy_v2[, 'Condition'])
pSites_with_clusters_fromHM_Normalized_tidy_v2$Time_Or_Treatment <- as.factor(pSites_with_clusters_fromHM_Normalized_tidy_v2[, 'Time_Or_Treatment'])
pSites_with_clusters_fromHM_Normalized_tidy_v2$Replicate <- as.factor(pSites_with_clusters_fromHM_Normalized_tidy_v2[, 'Replicate'])

#Transforming all treatments to numbers to use geom_smooth to plot overall trends
#Naming mitotic condition as 150

pSites_with_clusters_fromHM_Normalized_tidy_TxAsNumeric <- plyr::rename(pSites_with_clusters_fromHM_Normalized, 
                                                                        c('Interphase_120_A_R' = 'Interphase_120_A', 'Interphase_120_B_R' = 'Interphase_120_B',
                                                                          'Interphase_030_B_R' = 'Interphase_030_B', 'Mitosis_CyB_C_R' = 'Mitosis_150_CyBC',
                                                                          'Mitosis_CyB_A' = 'Mitosis_150_CyBA', 'Mitosis_CSF' = 'Mitosis_150_CSFA', 'Mitosis_CSF_TR_R' = 'Mitosis_150_CSFB',
                                                                          'Mitosis_CyB_D' = 'Mitosis_150_CyBD', 
                                                                          'Protein_Amino.acid_Site_NumberPhospho' = 'ID'))

#changing several columns to factors

pSites_with_clusters_fromHM_Normalized_tidy_TxAsNumeric$Amino.acid <- as.factor(pSites_with_clusters_fromHM_Normalized_tidy_TxAsNumeric[, 'Amino.acid'])
pSites_with_clusters_fromHM_Normalized_tidy_TxAsNumeric$Multiplicity <- as.factor(pSites_with_clusters_fromHM_Normalized_tidy_TxAsNumeric[, 'Multiplicity'])
pSites_with_clusters_fromHM_Normalized_tidy_TxAsNumeric$Proteins <- as.factor(pSites_with_clusters_fromHM_Normalized_tidy_TxAsNumeric[, 'Proteins'])
pSites_with_clusters_fromHM_Normalized_tidy_TxAsNumeric$Leading.proteins <- as.factor(pSites_with_clusters_fromHM_Normalized_tidy_TxAsNumeric[, 'Leading.proteins'])
pSites_with_clusters_fromHM_Normalized_tidy_TxAsNumeric$Protein <- as.factor(pSites_with_clusters_fromHM_Normalized_tidy_TxAsNumeric[, 'Protein'])
pSites_with_clusters_fromHM_Normalized_tidy_TxAsNumeric$Cluster <- as.factor(pSites_with_clusters_fromHM_Normalized_tidy_TxAsNumeric[, 'Cluster'])

#From wide format to long format

pSites_with_clusters_fromHM_Normalized_tidy_TxAsNumeric_v1 <- melt(pSites_with_clusters_fromHM_Normalized_tidy_TxAsNumeric)

#Separate the variable into 'Condition', 'Time/Treatment', 'Replicate'

pSites_with_clusters_fromHM_Normalized_tidy_TxAsNumeric_v2 <- pSites_with_clusters_fromHM_Normalized_tidy_TxAsNumeric_v1 %>%
  separate(variable, c('Condition', 'Time_Or_Treatment', 'Replicate'), '_')

str(pSites_with_clusters_fromHM_Normalized_tidy_TxAsNumeric_v2)


#Chaging time column to numeric

pSites_with_clusters_fromHM_Normalized_tidy_TxAsNumeric_v2$Time_Or_Treatment <- as.numeric(as.character(pSites_with_clusters_fromHM_Normalized_tidy_TxAsNumeric_v2$Time_Or_Treatment))
pSites_with_clusters_fromHM_Normalized_tidy_TxAsNumeric_v2$Condition <- as.factor(pSites_with_clusters_fromHM_Normalized_tidy_TxAsNumeric_v2[, 'Condition'])
pSites_with_clusters_fromHM_Normalized_tidy_TxAsNumeric_v2$Replicate <- as.factor(pSites_with_clusters_fromHM_Normalized_tidy_TxAsNumeric_v2[, 'Replicate'])

#PLotting
#mcm4|S|31|___1
#mcm4|S|31|___2
#mcm4|T|23|___2
#mcm4|S|87|___3
#mcm4|S|88|___3
#mcm4|T|94|___3
#nup98|T|912|___1
#nup98|S|693|___1
#nup98|S|693|___2

y <- ggplot(data = subset(pSites_with_clusters_fromHM_Normalized_tidy_TxAsNumeric_v2, ID %in% c('mcm4|S|31|___1', 'mcm4|S|31|___2')), 
            aes(x = Time_Or_Treatment, y = value, color = ID)) +
  geom_point(size = 2, shape = 21, stroke = 1.5)+
  geom_smooth(se = TRUE, alpha = 0.2, span = 0.75)+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), 
        axis.line = element_line(colour = "black", size = 0.25), axis.text.x = element_text(size = 14, color = 'black'), axis.ticks = element_line(size = 0.25), 
        axis.text.y = element_text(size = 14, color = 'black'), axis.title = element_text(size = 16, color = 'black'), plot.title = element_text(size = 18), 
        legend.title = element_text(size = 9, color = 'black'), legend.text = element_text(size = 9, color = 'black'), legend.position = 'top') +  
  ylab('Normalized Intensity') +
  xlab('Time (min)') +
  ggtitle('mcm4 S31 AND S31|T23') 

y + scale_color_tron()

#loading hits per motif

load('Hits_Aurora_All.RData')
load('Hits_Cdc7_All.RData')
load('Hits_CdkMinimal_All.RData')
load('Hits_CdkFull_All.RData')
Hits_Aurora_All_List
Hits_Cdc7_All_List
Hits_Cdk_Minimal_All_List
Hits_Cdk_Full_All_List

#Filtering for only IDs with Aurora motif
pSites_with_clusters_fromHM_Normalized_tidy_TxAsNumeric_v2_Aurora_Only <- 
  subset(pSites_with_clusters_fromHM_Normalized_tidy_TxAsNumeric_v2, ID %in% Hits_Aurora_All_List)

#Filtering for only IDs with Cdc7 motif
pSites_with_clusters_fromHM_Normalized_tidy_TxAsNumeric_v2_Cdc7_Only <- 
  subset(pSites_with_clusters_fromHM_Normalized_tidy_TxAsNumeric_v2, ID %in% Hits_Cdc7_All_List)

#Filtering for only IDs with cdk minimal motif
pSites_with_clusters_fromHM_Normalized_tidy_TxAsNumeric_v2_Cdk_Minimal_Only <- 
  subset(pSites_with_clusters_fromHM_Normalized_tidy_TxAsNumeric_v2, ID %in% Hits_Cdk_Minimal_All_List)

#Filtering for only IDs with cdk full motif
pSites_with_clusters_fromHM_Normalized_tidy_TxAsNumeric_v2_Cdk_Full_Only <- 
  subset(pSites_with_clusters_fromHM_Normalized_tidy_TxAsNumeric_v2, ID %in% Hits_Cdk_Full_All_List)

#Trying to generate a column with motif identifier



pSites_with_clusters_fromHM_Normalized_tidy_TxAsNumeric_v2_Aurora_Only$Motif <- 'aurora'
pSites_with_clusters_fromHM_Normalized_tidy_TxAsNumeric_v2_Cdc7_Only$Motif <- 'cdc7'
pSites_with_clusters_fromHM_Normalized_tidy_TxAsNumeric_v2_Cdk_Minimal_Only$Motif <- 'cdk_min'
pSites_with_clusters_fromHM_Normalized_tidy_TxAsNumeric_v2_Cdk_Full_Only$Motif <- 'cdk_full'

pSites_with_clusters_fromHM_Normalized_tidy_TxAsNumeric_v2_WithMotif <- rbind(pSites_with_clusters_fromHM_Normalized_tidy_TxAsNumeric_v2_Aurora_Only, 
                                                                              pSites_with_clusters_fromHM_Normalized_tidy_TxAsNumeric_v2_Cdc7_Only, 
                                                                              pSites_with_clusters_fromHM_Normalized_tidy_TxAsNumeric_v2_Cdk_Minimal_Only, 
                                                                              pSites_with_clusters_fromHM_Normalized_tidy_TxAsNumeric_v2_Cdk_Full_Only)


#Plotting Per kinase motif



s <- ggplot(data = subset(pSites_with_clusters_fromHM_Normalized_tidy_TxAsNumeric_v2_Cdk_Minimal_Only, Cluster %in% c('cluster4', 'cluster5', 'cluster6')), 
            aes(x = Time_Or_Treatment, y = value, color = Cluster)) +
  #geom_jitter(size = 2, shape = 21, stroke = 1.3)+
  geom_smooth(se = TRUE, alpha = 0.2, span = 0.65, method = 'loess')+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), 
        axis.line = element_line(colour = "black", size = 0.25), axis.text.x = element_text(size = 14, color = 'black'), axis.ticks = element_line(size = 0.25), 
        axis.text.y = element_text(size = 14, color = 'black'), axis.title = element_text(size = 16, color = 'black'), plot.title = element_text(size = 18), 
        legend.title = element_text(size = 9, color = 'black'), legend.text = element_text(size = 9, color = 'black'), legend.position = 'top') +  
  ylab('Normalized Intensity') +
  xlab('Time (min)') +
  ggtitle('Cdk minimal motif') 
s
s + scale_color_tron()

z <- ggplot(data = subset(pSites_with_clusters_fromHM_Normalized_tidy_TxAsNumeric_v2_WithMotif, Cluster %in% c('cluster4', 'cluster5', 'cluster6')), 
            aes(x = Time_Or_Treatment, y = value, color = Motif)) +
  #geom_jitter(size = 2, shape = 21, stroke = 1.3)+
  geom_smooth(se = TRUE, alpha = 0.2, span = 0.65, method = 'loess')+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), 
        axis.line = element_line(colour = "black", size = 0.25), axis.text.x = element_text(size = 14, color = 'black'), axis.ticks = element_line(size = 0.25), 
        axis.text.y = element_text(size = 14, color = 'black'), axis.title = element_text(size = 16, color = 'black'), plot.title = element_text(size = 18), 
        legend.title = element_text(size = 9, color = 'black'), legend.text = element_text(size = 9, color = 'black'), legend.position = 'top') +  
  ylab('Normalized Intensity') +
  xlab('Time (min)') +
  ggtitle('Interphase phosphorylations 
          (clusters 4, 5, 6)') 
z
z + scale_color_tron()

# Danger: Gero's code ahead

z1 <- ggplot(data = subset(pSites_with_clusters_fromHM_Normalized_tidy_TxAsNumeric_v2_WithMotif, Cluster %in% c('cluster4', 'cluster5', 'cluster6') & Time_Or_Treatment<=120), 
             aes(x = Time_Or_Treatment, y = value, color = Motif)) +
  #geom_jitter(size = 2, shape = 21, stroke = 1.3)+
  geom_smooth(se = TRUE, alpha = 0.2, span = 0.65, method = 'loess')+
  theme_bw() + coord_cartesian(ylim=c(-1.1,1.1), xlim=c(0,130))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), 
        axis.line = element_line(colour = "black", size = 0.25),axis.ticks = element_line(size = 0.25),
        axis.text.x = element_text(size = 14, color = 'black'),axis.text.y = element_text(size = 14, color = 'black'),
        axis.title.y = element_text(size = 16, color = 'black'), axis.title.x = element_text(size = 16, color = 'black'), plot.title = element_text(size = 18), 
        legend.title = element_text(size = 9, color = 'black'), legend.text = element_text(size = 9, color = 'black'), legend.position = 'top') +  
  ylab('Normalized Intensity') +
  xlab('Time (min)') +
  # ggtitle('Interphase phosphorylations clusters 4, 5, 6)') + 
  scale_x_continuous(breaks=c(0,30,60,90,120))
z1 <- z1 + scale_color_tron()


z2 <- ggplot(data = subset(pSites_with_clusters_fromHM_Normalized_tidy_TxAsNumeric_v2_WithMotif, Cluster %in% c('cluster4', 'cluster5', 'cluster6') & Time_Or_Treatment>120),
             aes(x = Time_Or_Treatment, y = value, color = Motif)) +
  geom_point(stat = "summary",fun = "mean", size = 15,shape="-") +
  theme_bw() + coord_cartesian(ylim=c(-1.1,1.1), xlim=c(130,170))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), 
        axis.line = element_line(colour = "black", size = 0.25),axis.line.y = element_blank(),axis.ticks = element_blank(),
        axis.text.x = element_text(size = 14, color = 'black'),axis.text.y = element_blank(),
        axis.title.y = element_blank(), axis.title.x = element_blank(), plot.title = element_text(size = 18), 
        legend.title = element_text(size = 9, color = 'black'), legend.text = element_text(size = 9, color = 'black'), legend.position = 'top') +  
  scale_x_continuous(breaks=c(150),labels = c("Mitosis"))
z2 <- z2 + scale_color_tron()

ggarrange(z1,z2,ncol = 2,nrow = 1,align = "h",widths = c(13,4),common.legend = T)
