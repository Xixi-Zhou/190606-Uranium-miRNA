# Predict miRNA-mRNA interactions and analyze reactome pathways
# Xixi Zhou, UNM CoP, Dept. of PS: xzhou@salud.unm.edu
sapply(c("tidyverse", "ggpubr", "rio", "miRNAtap", "ReactomePA"), 
       library, character.only = T)

# 1. Import data, predict mRNA targets, and save. -----
# Read miRNA targets from excel spreadsheet, up- and down-regulated.
miRNAs <- import("Data.xlsx") %>% as_tibble()
# Predict mRNAs
predictions <- miRNAs$miRNA %>%
  map(getPredictedTargets, # Function to predict mRNA
      species = "mmu", # Mouse
      min_src = 2, # at least from 2 sources
      method = "geom") %>%
  map(as_tibble, rownames = "mRNA") %>% # Data wrangling
  list(miRNAs$miRNA, miRNAs$regulation) %>%
  pmap(~ mutate(..1, miRNA = ..2, regulation = ..3)) %>% 
  bind_rows() %T>% # %T>% passes the value of left hand side
  export("mRNA_predictions.xlsx") # Save results to a spreadsheet

# Filter up- and down-regulated mRNAs, store as seperated vectors
# Rank value as the vector, mRNA names as vector name
mRNA_list <- predictions %>%
  split(f = predictions$regulation) %>%
  map(~ {
    vector <- pull(., rank_product)
    names(vector) <- pull(., mRNA)
    vector})

# 2. Reactome analysis -----
pathways <- mRNA_list %>%
  map( ~ enrichPathway(
    gene = names(.x), organism = "mouse", pvalueCutoff = 0.05, readable = T))

# 3. Plot -----
# Save data image
save.image("miRNA_data.RData")
# Plot to pdf
pdf("Figures.pdf", width = 16, height = 9, onefile = T)
# Plot histogram of mRNA prediction
gghistogram(
  predictions, x="miRNA", stat="count",
  title="Predicted mRNA Targets of miRNAs", fill="regulation",
  palette="lancet", color=F, xlab = "", ylab = "Count", legend="bottom")

# Plot reactome analysis results
map(pathways, dotplot, showCategory = 20)
map(pathways, cnetplot, categorySize="pvalue")
dev.off()
