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

# write.table(MLO_lilaDB_nr,"MLO_lilaDB_nr_260521.tab",sep = "\t",quote = F,row.names = F)


MLO_lilaDB_nr_curated <- read_delim("MLO_lilaDB_nr_260521_curated.tab", "\t", escape_double = FALSE, trim_ws = TRUE)

humanCDKtarget_MLOS <- subset(MLO_lilaDB_nr_curated, human=="Cdk target")
humanCDKtarget_MLOS$dominant_body <- as.factor(humanCDKtarget_MLOS$dominant_body)

xenopusDynamic_MLOS <- subset(MLO_lilaDB_nr_curated, xenopus=="Dynamic" ) 
xenopusDynamic_MLOS$dominant_body <- as.factor(xenopusDynamic_MLOS$dominant_body)

# Checking for other kinases: human-data must be loaded! 

mapk_targets <- subset(human_data,target_mapk=="mapk target")$`ACC#`
humanMAPKtarget_MLOS <- subset(MLO_lilaDB_nr_curated, `ACC#` %in% mapk_targets )
humanMAPKtarget_MLOS$dominant_body <- as.factor(humanMAPKtarget_MLOS$dominant_body)

aurk_targets <- subset(human_data,target_aurk=="aurk target")$`ACC#`
humanAURKtarget_MLOS <- subset(MLO_lilaDB_nr_curated, `ACC#` %in% aurk_targets )
humanAURKtarget_MLOS$dominant_body <- as.factor(humanAURKtarget_MLOS$dominant_body)

plk_targets <- subset(human_data,target_plk=="plk target")$`ACC#`
humanPLKtarget_MLOS <- subset(MLO_lilaDB_nr_curated, `ACC#` %in% plk_targets )
humanPLKtarget_MLOS$dominant_body <- as.factor(humanPLKtarget_MLOS$dominant_body)

nek_targets <- subset(human_data,target_nek=="nek target")$`ACC#`
humanNEKtarget_MLOS <- subset(MLO_lilaDB_nr_curated, `ACC#` %in% nek_targets )
humanNEKtarget_MLOS$dominant_body <- as.factor(humanNEKtarget_MLOS$dominant_body)

dyrk_targets <- subset(human_data,target_dyrk=="dyrk target")$`ACC#`
humanDYRKtarget_MLOS <- subset(MLO_lilaDB_nr_curated, `ACC#` %in% dyrk_targets )
humanDYRKtarget_MLOS$dominant_body <- as.factor(humanDYRKtarget_MLOS$dominant_body)