library(dplyr)
library(stringr)
library(readr)
library(tidyr)



MLO_lilaDB_nr <- read_delim("MLO_lilaDB_nr.tab", "\t", escape_double = FALSE, trim_ws = TRUE)

MLO_lilaDB_nr <- MLO_lilaDB_nr %>% group_by(uniprot) %>% summarise_at("MLOs",function(x){paste(x,collapse = ",")})

# write.table(MLO_lilaDB_nr,file = "MLO_lilaDB_nr_reformattted.tab",quote = F,sep = "\t",row.names = F)

MLO_lilaDB_kinaseTargets <- merge.data.frame(MLO_lilaDB_nr,human_data[,c("ACC#","target","target_mapk","target_aurk","target_plk","target_nek","target_dyrk")],by.x = "uniprot",by.y = "ACC#",all.x = T)

MLO_lilaDB_kinaseTargets <- MLO_lilaDB_kinaseTargets[!is.na(MLO_lilaDB_kinaseTargets$target),]
