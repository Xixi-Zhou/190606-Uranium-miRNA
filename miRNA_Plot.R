# Plot GO enrichment charts and reactome analysis charts
# Work with data from miRNA_Network_Data.R
# Xixi Zhou
library(topGO)
library(ReactomePA)
library(org.Mm.eg.db)

load(file = "miRNA_data.RData")
pdf("Plots.pdf", width = 18, height = 9)
showSigOfNodes(GO_data_up, score(ks_results_up), useInfo = 'all')
showSigOfNodes(GO_data_down, score(ks_results_down), useInfo = 'all')
dotplot(pathway_up, showCategory=15)
dotplot(pathway_down, showCategory=15)
cnetplot(pathway_up, categorySize="pvalue")
cnetplot(pathway_down, categorySize="pvalue")
dev.off()