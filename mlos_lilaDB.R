library(dplyr)
library(stringr)
library(readr)
library(tidyr)

redundant_MLO_lilaDB <-  read_delim("MLO_lilaDB_redundant_20052021.tab", "\t", escape_double = FALSE, trim_ws = TRUE)

MLO_lilaDB_nr <- redundant_MLO_lilaDB %>% group_by(`ACC#`) %>% summarise_at("body",function(x){paste(unique(x),collapse = ",")})



dynamic_symbol_to_HsUniprot <- read_delim("utrech/ids/dynamic_symbol_to_HsUniprot.tab", "\t", escape_double = FALSE, trim_ws = TRUE)

MLO_lilaDB_nr <- MLO_lilaDB_nr %>% mutate(xenopus=case_when(
    `ACC#` %in% dynamic_symbol_to_HsUniprot$`Human Uniprot` ~ "Dynamic",
    TRUE ~ "Non Dynamic"
))

cdkTargets_uniprot <- read_lines("utrech/ids/humanCDKtargets_uniprot.txt",skip = 1)

MLO_lilaDB_nr <- MLO_lilaDB_nr %>% mutate(human=case_when(
  `ACC#` %in% cdkTargets_uniprot ~ "Cdk target",
  TRUE ~ "Non Cdk target"
))

write.table(MLO_lilaDB_nr,"MLO_lilaDB_nr_260521.tab",sep = "\t",quote = F,row.names = F)
