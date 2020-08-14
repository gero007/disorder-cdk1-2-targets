#install packages

# install.packages("tidyverse")
# install.packages("gplots")
# install.packages("devtools")
# library(devtools)
# install_github("jokergoo/ComplexHeatmap")
# install.packages("heatmap3")

#Load packages#


library(tidyr)
library(gplots)
library(dplyr)
library(tibble)
library(ggplot2)
library(ComplexHeatmap)
library(heatmap3)
library(reshape2)

#Load data (beware that had to remove first useless three rows, otherwise everything was loaded as character)

pSites_Median_Norm_Full_Matrix <- read.csv('utrech/Phosphoproteomics_Shotgun_Clustering_Gero/All_pSites_ANOVA_Statistics.csv', stringsAsFactors = F)
str(pSites_Median_Norm_Full_Matrix)
summary(pSites_Median_Norm_Full_Matrix)
head(pSites_Median_Norm_Full_Matrix)

#Creating unique useful identifier

pSites_Median_Norm_Full_Matrix_Unique_ID <- unite(pSites_Median_Norm_Full_Matrix, 
                                                  'Protein_Amino.acid_Site_NumberPhospho', 
                                                  c('Proteins', 'Amino.acid', 'Positions.within.proteins', 'Multiplicity'),
                                                  sep ='|', remove = FALSE)

#filtering for ANOVA +

pSites_Median_Norm_Full_Matrix_Unique_ID_ANOVA_sig <- 
  pSites_Median_Norm_Full_Matrix_Unique_ID[pSites_Median_Norm_Full_Matrix_Unique_ID$ANOVA.Significant == '+', ] 

#Df with only mean intensities of time point of interest

pSites_Median_Norm_Full_Matrix_Unique_ID_ANOVA_sig_Mean <- 
  pSites_Median_Norm_Full_Matrix_Unique_ID_ANOVA_sig[, c(86:103)]

summary(pSites_Median_Norm_Full_Matrix_Unique_ID_ANOVA_sig_Mean)
str(pSites_Median_Norm_Full_Matrix_Unique_ID_ANOVA_sig_Mean)
typeof(pSites_Median_Norm_Full_Matrix_Unique_ID_ANOVA_sig_Mean)

#Preparing data
#Transforming to matrix

pSites_Median_Norm_Full_Matrix_Unique_ID_ANOVA_sig_Mean_Matrix <-
  as.matrix(sapply(pSites_Median_Norm_Full_Matrix_Unique_ID_ANOVA_sig_Mean, as.numeric))

#Changing column names

colnames(pSites_Median_Norm_Full_Matrix_Unique_ID_ANOVA_sig_Mean_Matrix) <- c('0h00', '0h15','0h30','0h45','1h00','1h15','1h25',
                                                                              '1h45','1h50','1h55','2h05','2h20','2h30','2h50',
                                                                              '3h00','3h10','3h20','3h35')

#Changing rownames

rownames(pSites_Median_Norm_Full_Matrix_Unique_ID_ANOVA_sig_Mean_Matrix) <- pSites_Median_Norm_Full_Matrix_Unique_ID_ANOVA_sig[, 68]


#Plotting with Complex Heatmap
#Scaling is needed, not done automatically by Complex HeatMap

pSites_Median_Norm_Full_Matrix_Unique_ID_ANOVA_sig_Mean_Matrix_Scaled <-
  t(scale(t(pSites_Median_Norm_Full_Matrix_Unique_ID_ANOVA_sig_Mean_Matrix)))

#Plotting using pearson, average, seems to give the best clustering thus far

Heatmap(pSites_Median_Norm_Full_Matrix_Unique_ID_ANOVA_sig_Mean_Matrix_Scaled, 
        column_order = colnames(pSites_Median_Norm_Full_Matrix_Unique_ID_ANOVA_sig_Mean_Matrix),
        show_row_names = FALSE, clustering_distance_rows = 'pearson', clustering_method_rows = 'average')

#spearman, average, splits the oscillating cluster

Heatmap(pSites_Median_Norm_Full_Matrix_Unique_ID_ANOVA_sig_Mean_Matrix_Scaled, 
        column_order = colnames(pSites_Median_Norm_Full_Matrix_Unique_ID_ANOVA_sig_Mean_Matrix),
        show_row_names = FALSE, clustering_distance_rows = 'spearman', clustering_method_rows = 'average')

#trying k-means; so far this is the best option
#it does a previous round of splitting using k-means. This allows to better group the features and additionally 
#highlight the patterns
#When doing k-means splitting, spearman seems to give the best clustering (average, complete, single look all fine). Pearson also
#works fine

HM = Heatmap(pSites_Median_Norm_Full_Matrix_Unique_ID_ANOVA_sig_Mean_Matrix_Scaled, 
             column_order = colnames(pSites_Median_Norm_Full_Matrix_Unique_ID_ANOVA_sig_Mean_Matrix),
             show_row_names = FALSE, row_km = 4, row_km_repeats = 100, clustering_distance_rows = 'spearman', 
             clustering_method_rows = 'average', name = 'Z-Score')

#Extracting row clusters
#From https://github.com/jokergoo/ComplexHeatmap/issues/136

HM = draw(HM)           #assign Heat Map to object and draw it
r.dend <- row_dend(HM)  #Extract row dendrogram
rcl.list <- row_order(HM)  #Extract clusters (output is a list)
lapply(rcl.list, function(x) length(x))  #check/confirm size clusters

# loop to extract genes for each cluster.
for (i in 1:length(row_order(HM))){
  if (i == 1) {
    clu <- t(t(row.names(pSites_Median_Norm_Full_Matrix_Unique_ID_ANOVA_sig_Mean_Matrix_Scaled[row_order(HM)[[i]],])))
    out <- cbind(clu, paste("cluster", i, sep=""))
    colnames(out) <- c("ID", "Cluster")
  } else {
    clu <- t(t(row.names(pSites_Median_Norm_Full_Matrix_Unique_ID_ANOVA_sig_Mean_Matrix_Scaled[row_order(HM)[[i]],])))
    clu <- cbind(clu, paste("cluster", i, sep=""))
    out <- rbind(out, clu)
  }
}

out