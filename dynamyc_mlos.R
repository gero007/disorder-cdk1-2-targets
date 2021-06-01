library(dplyr)
library(stringr)
library(readr)
library(tidyr)

dynamic_symbol_to_HsUniprot <- read_delim("utrech/ids/dynamic_symbol_to_HsUniprot.tab",
                                          "\t", escape_double = FALSE, trim_ws = TRUE)

MLO_lilaDB_nr <- read_delim("MLO_lilaDB_nr.tab", "\t", escape_double = FALSE, trim_ws = TRUE)

MLO_human <- read_delim("MLO_human.tab", "\t", escape_double = FALSE, trim_ws = TRUE)

MLO_lilaDB_expanded <- merge.data.frame(MLO_lilaDB_nr,MLO_human,all = T,by.x = "uniprot",by.y = "ACC#")

MLO_lilaDB_expanded <- MLO_lilaDB_expanded %>% unite("MLO", body:MLOs, na.rm = TRUE, remove = TRUE,sep = ",")

MLO_lilaDB_expanded <- MLO_lilaDB_expanded %>% group_by(uniprot) %>% summarise_at("MLO",function(x){paste(unique(x),collapse = ",")})

MLO_lilaDB_expanded <- subset(MLO_lilaDB_expanded,MLO!="")

MLO_lilaDB_expanded_dynamicSubset <- merge.data.frame(dynamic_symbol_to_HsUniprot,MLO_lilaDB_expanded,by.x = "Human Uniprot",by.y = "uniprot",all.x = T)


write.table(MLO_lilaDB_expanded_dynamicSubset,file = "MLO_lilaDB_expanded_dynamicSubset.tab",sep = "\t",quote = F,row.names = F)




