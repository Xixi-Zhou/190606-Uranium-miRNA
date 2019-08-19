# Predict miRNA-mRNA interactions, perform GO enrichment, 
# & analyze reactome pathways
# Xixi Zhou, UNM CoP, Dept. of PS: xzhou@salud.unm.edu
# Load packages, 'walk' and 'map' are not available for now.
sapply(c("tidyverse", "ggpubr", "rio",
         "miRNAtap", "topGO", "ReactomePA"), 
       library, character.only = T)
# Read miRNA targets from excel spreadsheet, up- and down-regulated.
miRNAs <- import("Data.xlsx") %>% as_tibble()

# 1. Predict mRNA targets, and save -----
predictions <-
  map(miRNAs$miRNA, getPredictedTargets, # Function to predict mRNA
      species = "mmu", # Mouse
      min_src = 2, # at least from 2 sources
      method = "geom") %>%
  map(as_tibble, rownames = "mRNA") %>% # Data wrangling
  list(miRNAs$miRNA, miRNAs$regulation) %>%
  pmap(~mutate(..1, miRNA = ..2, regulation = ..3)) %>% 
  bind_rows() %T>% # !!! %T>% passes left hand site value back
  export("mRNA_predictions.xlsx") # Save results to a spreadsheet
# Filter up- and down-regulated mRNAs, store as seperated vectors
# Rank value as the vector, mRNA names as vector name
mRNA_list <- split(predictions, f = predictions$regulation) %>%
  map(~{
    vector <- pull(., rank_product)
    names(vector) <- pull(., mRNA)
    vector})

# 2. GO enrichment of predicted genes -----
GO_data <- map(mRNA_list, ~new(
  'topGOdata', # Database
  ontology = 'BP', # Biological processes
  allGenes = ., annot = annFUN.GO2genes,
  GO2genes = annFUN.org(
    whichOnto = "BP", feasibleGenes = NULL,
    mapping = "org.Mm.eg.db", ID = "entrez"),
  geneSel = function(x) TRUE, # No filters
  nodeSize = 10))
ks_results <- map(GO_data, runTest, algorithm = "classic", statistic = "ks")
# Generate table of GO results, and Save GO results
# 'export' function can write the list of tables to multiple sheets
map2(GO_data, ks_results,
     ~GenTable(.x, KS = .y, orderBy = "KS", topNodes = 50)) %>% 
  export(file = "GO_results.xlsx")

# 3. Use predicted genes for reactome analysis -----
pathways <- map(mRNA_list, ~enrichPathway(
  gene = names(.), organism = "mouse", pvalueCutoff = 0.05, readable = T))

# 4. generate figures -----
# Save data image
save.image("miRNA_data.RData")
# Plot to pdf
pdf("Figures.pdf", width = 16, height = 9, onefile = T)
# Plot histogram of mRNA prediction
gghistogram(
  predictions, x="miRNA", stat="count",
  title="Predicted mRNA Targets of miRNAs", fill="regulation",
  palette="lancet", color=F, xlab = "", ylab = "Count", legend="bottom")
# Plot GO enrichment results
map2(GO_data, ks_results, ~showSigOfNodes(.x, score(.y), useInfo = 'all'))
# Plot reactome analysis results
map(pathways, dotplot, showCategory = 20)
map(pathways, cnetplot, categorySize="pvalue")
dev.off()
