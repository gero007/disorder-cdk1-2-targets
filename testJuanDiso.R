library(iupred)
library(readr)
library(dplyr)
library(data.table)
library(ggplot2)


all_protein_id_phospho <- read_lines("utrech/Xen_phospho_allProteins.txt")
cluster_oscillating <- read_lines("utrech/Xen_phospho_ClusterD.txt")
# cluster_exit_metaphase <- read_file("")
significant_anova <- read_lines("utrech/Xen_phospho_ANOVApos.txt")


# # Read in pre-calculated file or perform iupred
# if (!file.exists("all-predictions.Rdata")) {
#   all_predictions <-
#     iupred.predict.fastafile("Xenopus_Database_MPI.FASTA",
#                              restrict = unique(all_protein_id_phospho
#                                                $Protein),
#                              idPattern = "^[^|]+\\|([^ ]+).*")
#   save(all_predictions, file="all-predictions.Rdata")
# } else {
#   print ("Loading prepared file")
#   load("all-predictions.Rdata")
# }
all_predictions <- list()
disoPath <- "predictions/IUpred_run_xenopus/"
file.names <- dir(disoPath, pattern =".iupred")
for (disofile in file.names) {
  id<-strsplit(x = disofile,split = ".",fixed = T)[[1]][1]
  aux_table <- read_delim(paste(disoPath,disofile,sep = ""),"\t", escape_double = FALSE, col_names = FALSE,col_types = cols(X1 = col_integer(),X2 = col_character(),X3 = col_double()),comment = "#", trim_ws = TRUE)
  colnames(aux_table) <- c("POS","RES","IUPRED_SCORE")
  aux_table$IUPRED_DISO <- aux_table$IUPRED_SCORE>=0.5
  all_predictions[[id]]$`.sequence` <- paste(aux_table$RES,collapse = "")
  all_predictions[[id]]$`.positions` <- aux_table$POS
  all_predictions[[id]]$`.scores` <- aux_table$IUPRED_SCORE
  all_predictions[[id]]$`.phospho` <- NULL
  all_predictions[[id]]$`.accesion` <- id
  all_predictions[[id]]$`.threshold` <- 0.5
  all_predictions[[id]]$`.disordered` <- aux_table$IUPRED_DISO
} 

# remove unannotated proteins
all_protein_id_phospho <- all_protein_id_phospho[all_protein_id_phospho %in% names(all_predictions)]

# Introduce phosphosites (from all_proteins variable) using the apply function
apply(all_protein_id_phospho, MARGIN=1, FUN=function(x) {
  all_predictions[[x[["Protein"]]]] <<- setPhospho(all_predictions[[x[["Protein"]]]] , as.numeric(x[["p-site"]]));
  return(invisible());
})

# Get the disordered/ordered phospho counts using countPhosphoDisordered
phospho_relative_to_organization<- lapply(all_predictions, countPhosphoDisordered)

# unlist into data.table
phospho_relative_to_organization <- rbindlist(phospho_relative_to_organization, id="protein")

# Convenient structure for batch processing
phospho_region <- list(
  all            = phospho_relative_to_organization,
  anova_sig      = phospho_relative_to_organization[protein %in% significant_anova],
  oscillating    = phospho_relative_to_organization[protein %in% cluster_oscillating,],
  exit_metaphase = phospho_relative_to_organization[protein %in% cluster_exit_metaphase,],
  nonsignificant = phospho_relative_to_organization[! protein %in% significant_anova,]
)

# Calculate percentage of disordered region from all proteins
disordered_fraction <- mapply(all_predictions, FUN=function(x) {sum(x$.disordered) / stringi::stri_length(x$.sequence)})

lapply(phospho_region, function(x) {
  x[,expected_disordered := (disordered + ordered)  * disordered_fraction[protein] ];
  x[,expected_ordered := (disordered + ordered)  * (1-disordered_fraction[protein]) ];
})


long_phospho_region_list <- lapply(phospho_region, FUN=melt, id="protein")
long_phospho_region <- rbindlist(long_phospho_region_list, id="group")

wide_phospho_region <- rbindlist(phospho_region, id="group")

# Binomial test
wide_phospho_region[,disordered_prop := disordered / (disordered + ordered)]
wide_phospho_region[,expect_disordered_prop := expected_disordered / (expected_disordered + expected_ordered)]
wide_phospho_region[,n := disordered + ordered]

wide_phospho_region <- wide_phospho_region[n > 0]

wide_phospho_region[,binom := purrr::pmap(.(disordered, n, expect_disordered_prop), binom.test, alternative="t")]
wide_phospho_region[,binom_p := mapply("[[", binom, "p.value", SIMPLIFY = T)]
wide_phospho_region[,binom_q := mapply(p.adjust, binom_p)]
wide_phospho_region[,binom_sig := factor(ifelse(binom_q < 0.05, ifelse(binom_q < 0.01, "1% FDR", "5% FDR"), "n.s."), levels=c("n.s.","5% FDR", "1% FDR")) ]


ggplot(long_phospho_region, aes(x=variable, y=value, fill=variable)) + geom_boxplot() + facet_grid(. ~ group)

ggplot(wide_phospho_region) + 
  geom_point(aes(x=disordered_prop, y=expect_disordered_prop, color=binom_sig)) + 
  geom_abline(slope=1, color="red") + 
  theme_minimal()

ggplot(wide_phospho_region) +
  geom_point(aes(x=disordered, y=expected_disordered, colour = binom_sig), alpha = 0.5, shape=16) +
  geom_abline(slope=1, color="red") +
  facet_grid(. ~ group ) + 
  scale_colour_manual(values = c("gray","orange","red")) + 
  theme_minimal()




