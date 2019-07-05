# Purpose:
# predict miRNA-mRNA interactions;
# perform GO enrichment;
# analyze reactome pathways.

# Xixi Zhou, UNM CoP, Dept. of PS: xzhou@salud.unm.edu

# Load packages:
library(miRNAtap) # mRNA prediction for miRNA targets
library(topGO) # Gene Ontology enrichment
library(ReactomePA) # Reactome pathway analysis
library(org.Mm.eg.db) # Genome wide annotation for Mouse
library(tidyverse) # Basic data science package
library(readxl) # Read excel spredsheets
library(ggpubr) # Generate figures for publication

# Read miRNA targets from two txt files, up- and down-regulated.
mir <- read_xlsx("Data.xlsx", sheet = "miRNAs", col_names = T)

# 1. Predict mRNA targets:
#   Write gene predictions from each miRNA to a txt file under "Targets" folder.
#   Save lists of predicted genes for GO analysis.
predictions <- NULL
for (i in 1:nrow(mir)) {
  targets <- getPredictedTargets(mir$miRNA[i],
                                 # For each miRNA
                                 species = "mmu",
                                 # Mouse
                                 min_src = 2,
                                 # found at least in 2 predictions
                                 method = "geom") # geom method, by default
  # "targets" is returned as a matrix with geneID as rownames
  # Now convert to tibble, with geneID as a new column
  targets <- as_tibble(targets, rownames = "mRNA")
  targets <- add_column(targets,
                        miRNA = mir$miRNA[i],
                        # Add miRNA name as a new column
                        Regulation = mir$Regulation[i]) # Annotate up or down
  predictions <- bind_rows(predictions, targets)
}

# Filter up- and down-regulated mRNAs, store as seperated vectors
# Rank value as the vector, mRNA names as vector name
predictions_list <- list()
regulation <- levels(factor(predictions$Regulation))
for (i in regulation){
  predictions_table <- filter(predictions, Regulation==i)
  predictions_vector <- pull(predictions_table, rank_product)
  names(predictions_vector) <- pull(predictions_table, mRNA)
  predictions_list[[length(predictions_list)+1]] <- predictions_vector
}
# Clean up objects
rm(i, targets,
   predictions_table,
   predictions_vector) # Delete temperary objects

# 2. GO enrichment of predicted genes
selection <- function(x) {TRUE} # We don't impose a cut off
all_GO2_genes <- annFUN.org(
  whichOnto = "BP",
  feasibleGenes = NULL,
  mapping = "org.Mm.eg.db",
  ID = "entrez"
)
# Set up lists to store results
GO_data_list <- list()
GO_results_list <- list()
ks_results_list <- list()

for (i in predictions_list) {
  GO_data <- new(
    'topGOdata',
    ontology = 'BP',
    allGenes = i,
    annot = annFUN.GO2genes,
    GO2genes = all_GO2_genes,
    geneSel = selection,
    nodeSize = 10
  )
  GO_data_list[[length(GO_data_list)+1]] <- GO_data
  ks_results <- runTest(GO_data,
                       algorithm = "classic",
                       statistic = "ks")
  ks_results_list[[length(ks_results_list)+1]] <- ks_results
  GO_results <- GenTable(GO_data,
                        KS = ks_results,
                        orderBy = "KS",
                        topNodes = 50)
  GO_results_list[[length(GO_results_list)+1]] <- GO_results
}
# Clean up
rm(i, selection, all_GO2_genes, GO_data, GO_results, ks_results)

# 3. Use predicted genes for reactome analysis
pathways_list <- list()
for (i in predictions_list) {
  pathways <- enrichPathway(
    gene = names(i),
    organism = "mouse",
    pvalueCutoff = 0.05, # Cut off at p<0.05
    readable = T
  )
  pathways_list[[length(pathways_list)+1]] <- pathways
}
# Clean up
rm(i, pathways)
# 4. Save results and generate figures
# Save data image
save.image("miRNA_data.RData")
# Save mRNA targets
write.table(
  predictions,
  file = "Results/Predictions.txt",
  append = F,
  row.names = F,
  col.names = T,
  sep = "\t",
  dec = "."
)
for (i in 1:length(predictions_list)) {
  write.table(
    GO_results_list[[i]],
    file = paste("Results/GO_", regulation[i], ".txt", sep = ""),
    append = F,
    sep = "\t",
    dec = ".",
    row.names = F,
    col.names = T
  )
}
# # plot device
pdf("Figures.pdf", width = 18, height = 9, onefile = T)
# Summary targets in a histogram
gghistogram(predictions, x="miRNA", stat="count",
            title="Predicted mRNA Targets of miRNAs",
            fill="Regulation", palette="lancet",
            color=FALSE, xlab = "", ylab = "Count",
            legend="bottom")
# Plot GO and pathway figures
showSigOfNodes(GO_data_list[[1]],
               score(ks_results_list[[1]]),
               useInfo = 'all')
mtext(regulation[1])
dotplot(pathways_list[[1]], showCategory = 15)
cnetplot(pathways_list[[1]], categorySize = "pvalue")
showSigOfNodes(GO_data_list[[2]],
               score(ks_results_list[[2]]),
               useInfo = 'all')
mtext(regulation[2])
dotplot(pathways_list[[2]], showCategory = 15)
cnetplot(pathways_list[[2]], categorySize = "pvalue")
# Not sure why loop doesn't work for plotting
dev.off()
