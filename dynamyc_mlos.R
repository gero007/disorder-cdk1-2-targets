
unique(MLO_lilaDB_nr$uniprot)
length(unique(MLO_lilaDB_nr$uniprot))
cluster_A <- read_lines("utrech/ids/xenopusclusterA_mpi.txt")
cluster_B <- read_lines("utrech/ids/xenopusclusterB_mpi.txt")
cluster_C <- read_lines("utrech/ids/xenopusclusterC_mpi.txt")
cluster_D <- read_lines("utrech/ids/xenopusclusterD_mpi.txt")
significant_anova <- Reduce(union, list(cluster_A,cluster_B,cluster_C,cluster_D))
human_CDK1targets <- read_lines("utrech/ids/humanCDKtargets_mpi.txt")
xenopus_ANOVA_data <- read_delim("utrech/Xen_phospho_AnovaPhosphosites.txt","\t", escape_double = FALSE, trim_ws = TRUE)
subset(MLO_lilaDB_nr, uniprot %in% significant_anova)
xenopus_ANOVA_data[1:10,]
library(readr)
dynamic_symbol_to_HsUniprot <- read_delim("utrech/ids/dynamic_symbol_to_HsUniprot.tab",
"\t", escape_double = FALSE, trim_ws = TRUE)
View(dynamic_symbol_to_HsUniprot)

lilaDB_dynamicSubset <- subset(MLO_lilaDB_nr, uniprot %in% dynamic_symbol_to_HsUniprot$`Human Uniprot`)
lilaDB_dynamicSubset %>% group_by(uniprot) %>% summarise(str_c(`MLOs`, collapse = ";"))