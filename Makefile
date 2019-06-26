VPATH = Source

.PHONY: all
all: Plots.pdf miR-7a.pdf

miR-7a.pdf: miR-7a.R
	Rscript --vanilla $<

Plots.pdf: miRNA_Plot.R miRNA_data.RData
	Rscript --vanilla $<

miRNA_data.RData: miRNA_Network_Data.R miRNAs_up.txt miRNAs_down.txt
	Rscript --vanilla $<
