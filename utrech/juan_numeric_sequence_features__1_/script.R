options(stringsAsFactors = F)

library(tidyverse)
library(ggplot2)
library(magrittr)
library(heatmap3)
library(viridisLite)
library(caret)

## Read data xenopus

xenl.features <- read_delim("calculator/Xenopus_Database_MPI.FASTA.features", delim="\t")
xenl.features %<>% select(-net_charge_P) %>% modify_if(is.double, function(x) {(x - mean(x)) / sd(x)} )
xenl.features %<>% mutate(TRG_ER_FFAT_1=0, LIG_GLEBS_BUB3_1=0) # correction for an omission in python code - fixed there now

## Read data yeast
elife <- read.delim("input/Clusterplot_evolsig.cdt", stringsAsFactors=FALSE)
elife <- elife[2:nrow(elife), 2:ncol(elife)]
rownames(elife) <- elife$idr_name
elife %<>% select(-idr_name, -GWEIGHT, -NAME, -mean_calc_net_charge_P)
elife %<>% select(starts_with("mean_calc"))
names(elife) <- gsub("mean_calc_","",names(elife))

xenl.coordinates <- xenl.features[,"IDR"]
xenl.coordinates <- xenl.coordinates %>% separate("IDR", c("id", "protein"), sep="\\|", remove = F) %>% separate("protein", c("protein", "start","end"), sep="_")

expclusters <- read_csv("input/expression_clusters.csv")

# join 
# after https://stackoverflow.com/questions/37289405/dplyr-left-join-by-less-than-greater-than-condition
colortable <- xenl.coordinates %>% 
  left_join(expclusters, by="protein") %>%
  filter( psite - 1 >= start, psite - 1 <= end) %>%
  select(-id, -protein, -start, -end, -psite) %>%
  distinct(IDR, cluster, .keep_all = TRUE) %>%
  pivot_wider(names_from = cluster, values_from = cluster) %>%
  mutate_at(vars(-IDR), list (~ !is.na(.))) %>%
  mutate(color = if_else(Exit_Metaphase & Oscillating_Mitotic, "yellow", if_else(Exit_Metaphase, "blue","red")))

rowcolors <- xenl.features[,"IDR"] %>% left_join(colortable, by = "IDR")

### Sort the xenopus data the same way as the yeast data in a dataframe
xenl.features.df <- xenl.features[,names(elife)]
rownames(xenl.features.df) <- xenl.features$IDR

#svm <- train("svmLinear3", xenl.features.df, rowcolors$color)

pdf("output/heatmap_xenopus.pdf", width=11, height=11)
xenl.cluster <- heatmap3(xenl.features.df,
                         balanceColor = T,
                         #col=colorRampPalette(c("blue","blue","black","yellow","yellow"))(1024),
                         col=viridis(1024),
                         RowSideColors = rowcolors$color,
                         distfun = dist,
                         showColDendro = F,
                         Colv = NA,
                         keep.dendro=TRUE)

dev.off()

xenl.features %>% bind_cols(rowcolors %>% select(-IDR)) %>% write_csv("xenopus_features_colors.csv")


#####
## Check if the numbers are in the same ballpark as the yeast dataset

combined.dataset <- rbind(cbind(xenl.features.df, color="red"), cbind(elife, color="black"))

pdf("output/heatmap_combined.pdf", width=11, height=12)
heatmap3(combined.dataset %>% select(-color), balanceColor = T, col=colorRampPalette(c("blue","blue","black","yellow","yellow"))(1024), distfun = dist,showColDendro = F, Colv = NA, RowSideColors = combined.dataset$color)
dev.off()

write.xlsx(xenl.features.df, sheetName =  "Xenopus", file = "output/Xenopus_normalized.xlsx", col.names = T, row.names = T, showNA = T)
xenl.features %>% bind_cols(rowcolors %>% select(-IDR)) %>% xlsx::write.xlsx("output/xenopus_features.xlsx", sheetName = "Xenopus_colors")
write.xlsx(elife,sheetName =  "Zarin", file = "output/Zarin_normalized.xlsx")
