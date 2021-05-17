library(dplyr)
library(stringr)
library(readr)
library(tidyr)



MLO_lilaDB_nr <- read_delim("MLO_lilaDB_nr.tab", "\t", escape_double = FALSE, trim_ws = TRUE)

MLO_human <- read_delim("MLO_human.tab", "\t", escape_double = FALSE, trim_ws = TRUE)

MLO_lilaDB_expanded <- merge.data.frame(MLO_lilaDB_nr,MLO_human,all = T,by.x = "uniprot",by.y = "ACC#")

MLO_lilaDB_expanded <- MLO_lilaDB_expanded %>% unite("MLO", body:MLOs, na.rm = TRUE, remove = TRUE,sep = ",")

MLO_lilaDB_expanded <- MLO_lilaDB_expanded %>% group_by(uniprot) %>% summarise_at("MLO",function(x){paste(x,collapse = ",")})

MLO_lilaDB_expanded <- subset(MLO_lilaDB_expanded,MLO!="")

write.table(MLO_lilaDB_expanded,file = "MLO_lilaDB_nr_expanded.tab",quote = F,sep = "\t",row.names = F)

MLO_lilaDB_expanded_kinaseTargets <- merge.data.frame(MLO_lilaDB_expanded,human_data[,c("ACC#","target","target_mapk","target_aurk","target_plk","target_nek","target_dyrk")],by.x = "uniprot",by.y = "ACC#",all.x = T)

MLO_lilaDB_expanded_kinaseTargets <- MLO_lilaDB_expanded_kinaseTargets[!is.na(MLO_lilaDB_expanded_kinaseTargets$target),]
