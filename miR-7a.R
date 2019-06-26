# Analyze miR-7a
# Xixi Zhou
library(miRNAtap)
library(ReactomePA)

mir <- 'miR-7a'
targets <- getPredictedTargets(mir,
                               species = "mmu",
                               min_src = 2,
                               method = "geom")
predictions <- targets[, "rank_product"]

mRNA_names <- names(predictions)
mRNA_pathways <- enrichPathway(
  gene = mRNA_names,
  organism = "mouse",
  # 0.05 returns nothing
  pvalueCutoff = 0.9,
  qvalueCutoff = 0.9,
  readable = T
)

pdf(paste(mir, ".pdf", sep = ""), width = 18, height = 9)
dotplot(mRNA_pathways, showCategory = 20)
cnetplot(mRNA_pathways, categorySize = "pvalue")
dev.off()
