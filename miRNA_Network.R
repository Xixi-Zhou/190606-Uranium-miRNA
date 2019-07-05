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
miRNAs <- read_xlsx("Data.xlsx", sheet = "miRNAs", col_names = T)

# 1. Predict mRNA targets:
predictions <- lapply(miRNAs$miRNA,
                      getPredictedTargets, # Function to predict mRNA
                      species = "mmu", # Mouse
                      min_src = 2, # at least from 2 sources
                      method = "geom") # geom method by default
# Convert matrix to tibble, add gene ID as a new column
predictions <- lapply(predictions, as_tibble, rownames = "mRNA")
# Add annotations
predictions <- mapply(add_column, # Add annotations to results
                      predictions,
                      miRNA = miRNAs$miRNA, # miRNA name
                      regulation = miRNAs$regulation) # up or down
# Combine into a single tibble
predictions <- do.call("bind_rows", predictions) 

# Filter up- and down-regulated mRNAs, store as seperated vectors
# Rank value as the vector, mRNA names as vector name
mRNA_list <- split(predictions, f = predictions$regulation)
mRNA_list <- lapply(mRNA_list,
                    function(x) {
                      vector <- pull(x, rank_product)
                      names(vector) <- pull(x, mRNA)
                      vector})

# 2. GO enrichment of predicted genes
selection <- function(x) {TRUE} # We don't impose a cut off
all_GO2_genes <- annFUN.org(
  whichOnto = "BP",
  feasibleGenes = NULL,
  mapping = "org.Mm.eg.db",
  ID = "entrez")
GO_data <- lapply(mRNA_list,
                  function(x) {
                    new('topGOdata', # Database
                        ontology = 'BP', # Biological processes
                        allGenes = x,
                        annot = annFUN.GO2genes,
                        GO2genes = all_GO2_genes,
                        geneSel = selection,
                        nodeSize = 10)})
ks_results <- lapply(GO_data, runTest,
                     algorithm = "classic",
                     statistic = "ks")
# Generate table of GO results
GO_results <- mapply(function(x, y) GenTable(x,
                                             KS = y,
                                             orderBy = "KS",
                                             topNodes = 50),
                     GO_data,
                     ks_results,
                     SIMPLIFY = F) # Force to return a list
# Clean up
rm(all_GO2_genes, selection)

# 3. Use predicted genes for reactome analysis
pathways <- lapply(mRNA_list, function(x) {
  enrichPathway(
    gene = names(x),
    organism = "mouse",
    pvalueCutoff = 0.05, # Cut off at p<0.05
    readable = T)})

# 4. Save results and generate figures
# Save data image
save.image("miRNA_data.RData")
# Save mRNA targets
write.table(predictions,
            file = "Results/Predictions.txt",
            append = F,
            row.names = F,
            col.names = T,
            sep = "\t",
            dec = ".")
# Save GO results
mapply(function(x,y)
  write.table(x,
              file = paste("Results/GO_", y, ".txt", sep = ""),
              append = F,
              sep = "\t",
              dec = ".",
              row.names = F,
              col.names = T), GO_results, names(GO_results))
# plot device
pdf("Figures.pdf", width = 18, height = 9, onefile = T)
# Summary targets in a histogram
gghistogram(predictions, x="miRNA", stat="count",
            title="Predicted mRNA Targets of miRNAs",
            fill="regulation", palette="lancet",
            color=FALSE, xlab = "", ylab = "Count",
            legend="bottom")
# Plot GO and pathway figures
mapply(function(x,y) {
  showSigOfNodes(x, score(y), useInfo = 'all')},
  GO_data, ks_results)
lapply(pathways, dotplot, showCategory = 20)
lapply(pathways, cnetplot, categorySize="pvalue")
# Wrap up
dev.off()
