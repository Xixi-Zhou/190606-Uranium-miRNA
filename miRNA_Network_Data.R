# This script is to predict miRNA-mRNA interactions, perform GO enrichment, and analyze pathways.
# Xixi Zhou, UNM CoP, Dept. of PS: xzhou@salud.unm.edu

# Load packages:
#   miRNAtap: miRNA prediction
#   topGO: Gene Ontology enrichment
#   ReactomePA: Reactome pathway analysis
#   org.Mm.eg.db: Genome wide annotation for Mouse
library(miRNAtap)
library(topGO)
library(ReactomePA)
library(org.Mm.eg.db)

# Read miRNA targets from two txt files, up- and down-regulated.
mir_up <- read.table(file="Source/miRNAs_up.txt", header = TRUE)
mir_down <- read.table(file="Source/miRNAs_down.txt", header = TRUE)

# Predict gene targets:
#   Write gene predictions from each miRNA to a txt file under "Targets" folder.
#   Save lists of predicted genes for GO analysis.
predictions_up <- vector()
for (i in mir_up$Up_Regulated) {
  targets <- getPredictedTargets(i,
                                 species = "mmu",
                                 min_src = 2,
                                 method = "geom")
  ranked_genes <- targets[, "rank_product"]
  write.table(targets, file = paste("Targets/", i, ".txt", sep = ""),
              append = FALSE, sep = "\t", dec = ".",
              row.names = TRUE, col.names = TRUE)
  predictions_up <- c(predictions_up, ranked_genes)}

predictions_down <- vector()
for (i in mir_down$Down_regulated) {
  targets <- getPredictedTargets(i,
                                 species = "mmu",
                                 min_src = 2,
                                 method = "geom")
  ranked_genes <- targets[, "rank_product"]
  write.table(targets, file = paste("Targets/", i, ".txt", sep = ""),
              append = FALSE, sep = "\t", dec = ".",
              row.names = TRUE, col.names = TRUE)
  predictions_down <- c(predictions_down, ranked_genes)
}
rm(i,targets,ranked_genes)

# Use predicted genes for GO enrichment
selection <- function(x) TRUE # We don't impose a cut off
all_GO2_genes <- annFUN.org(whichOnto = "BP",
                            feasibleGenes = NULL,
                            mapping = "org.Mm.eg.db",
                            ID = "entrez")
GO_data_up <- new('topGOdata',
                  ontology = 'BP',
                  allGenes = predictions_up,
                  annot = annFUN.GO2genes,
                  GO2genes = all_GO2_genes,
                  geneSel = selection,
                  nodeSize = 10)
ks_results_up = runTest(GO_data_up,
                        algorithm = "classic",
                        statistic = "ks")
GO_results_up = GenTable(GO_data_up,
                         KS = ks_results_up,
                         orderBy = "KS",
                         topNodes = 20)
write.table(GO_results_up[,c('GO.ID','Term','KS')],
            file = "GO/GO_Up.txt",
            append = FALSE, sep = "\t", dec = ".",
            row.names = TRUE, col.names = TRUE)

GO_data_down <- new('topGOdata',
                    ontology = 'BP',
                    allGenes = predictions_down,
                    annot = annFUN.GO2genes,
                    GO2genes = all_GO2_genes,
                    geneSel = selection,
                    nodeSize = 10)
ks_results_down = runTest(GO_data_down,
                          algorithm = "classic",
                          statistic = "ks")
GO_results_down = GenTable(GO_data_down,
                           KS = ks_results_down,
                           orderBy = "KS",
                           topNodes = 20)
write.table(GO_results_down[,c('GO.ID','Term','KS')],
            file = "GO/GO_Down.txt",
            append = FALSE, sep = "\t", dec = ".",
            row.names = TRUE, col.names = TRUE)

# Using predicted genes for reactome analysis
names_up <- names(predictions_up)
names_down <- names(predictions_down)
pathway_up <- enrichPathway(gene = names_up,
                            organism = "mouse",
                            pvalueCutoff=0.05, # Cut off at p<0.05
                            readable=T)
pathway_down <- enrichPathway(gene = names_down,
                              organism = "mouse",
                              pvalueCutoff=0.05,
                              readable=T)

# Save data image for future analysis
save.image("miRNA_data.RData")