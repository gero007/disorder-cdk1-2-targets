library(iupred)
library(readxl)
library(dplyr)
library(data.table)
library(ggplot2)


all_protein_id_phospho <- read_xlsx("Phospho_Xlaevis_Proteins_PSites.xlsx", "All_IDs", col_names = T, col_types = c("text", "numeric"), range="A1:B4584")

cluster_oscillating = read_xlsx("Phospho_Xlaevis_Proteins_PSites.xlsx", "All_Oscillating", col_names = F, col_types = "text")[["...1"]]
cluster_exit_metaphase = read_xlsx("Phospho_Xlaevis_Proteins_PSites.xlsx", "All_Exit_Metaphase", col_names = T, col_types = "text")[["Protein"]]
significant_anova = read_xlsx("Phospho_Xlaevis_Proteins_PSites.xlsx", "All_ANOVA+", col_names = T, col_types = "text")[["Protein"]]


# Read in pre-calculated file or perform iupred
if (!file.exists("all-predictions.Rdata")) {
  all_predictions <-
    iupred.predict.fastafile("Xenopus_Database_MPI.FASTA",
                             restrict = unique(all_protein_id_phospho
                                              $Protein),
                             idPattern = "^[^|]+\\|([^ ]+).*")
  save(all_predictions, file="all-predictions.Rdata")
} else {
  print ("Loading prepared file")
  load("all-predictions.Rdata")
}


# remove unannotated proteins
all_protein_id_phospho <- all_protein_id_phospho[all_protein_id_phospho$Protein %in% names(all_predictions),]

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




