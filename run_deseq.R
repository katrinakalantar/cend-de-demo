
#
## SCRIPT PART 1 - LOADING IN THE DATA
#

setwd("/Users/katrina.kalantar/code/cend-de-demo")
getwd()

# Section 1.1 - Read in the metadata 
metatable <- read.csv("./data/rr057_metadata.csv", row.names=1)

dim(metatable)      # What are the dimensions of the metadata table?
head(metatable)     # general view of top few rows of the metadata
colnames(metatable) # What data is available in the metadata table?


# Section 1.2 - Read in the gene counts file
genecounts <- read.csv("data/rr057_genecounts.csv", row.names=1)

dim(genecounts)    # What are the dimensions of the gene counts table?
head(genecounts[,1:5])  # general view of the top few rows, first few columns of the gene counts

# ensure that the `viral_status` variable is coded appropriately
metatable$viral_status <- factor(metatable$viral_status, levels = c("SC2","no_virus","other_virus"))

# remove samples with fewer than 1 million gene counts
genecounts <- genecounts[,colSums(genecounts) > 1000000]

# make sure that the samples in the metadata match the samples in the gene counts matrix
metatable <- metatable[colnames(genecounts),]



# Section 1.3 - Read in the gene name annotations
#     note: the gene name annotations will be used to make human readable gene names

gene_annot_tsv <- read.table("data/gene2name.txt", row.names = 1, col.names = c("gene_ID","gene_symbol", "chr"))

head(gene_annot_tsv)   # general view of the information contained in gene_annot_tsv



#
## SCRIPT PART 2 - RUNNING DIFFERENTIAL EXPRESSION ANALYSIS
#


# load in the DESeq2 library
library("DESeq2")

# combine the genecounts and metadata with the design
dds <- DESeqDataSetFromMatrix(countData = genecounts,
                              colData = metatable,
                              design = ~ viral_status)

# filter to keep only the genes with at least 10 reads total 
keep_genes <- rowSums(counts(dds)) >= 10
dds <- dds[keep_genes,]

# Filter is to ensure at least N samples with a count of 10 or more
keep_genes <- rowSums(counts(dds) >= 10) >= 25
dds <- dds[keep_genes,]

dim(dds)

# Run differential expression using DESeq2
dds <- DESeq(dds)


#
## SCRIPT PART 3 - REVIEW THE RESULTS of SC2 v. no virus comparison
#


# Section 3.1 - Generate the results for comparison of ** SC2 v. no virus **

# output the DESeq2 results to a dataframe
res_sc2_v_no_virus <- results(dds, contrast=c("viral_status","SC2","no_virus"))

# sort the results by adjusted pvalue
res_sc2_v_no_virus_ordered <- res_sc2_v_no_virus[order(res_sc2_v_no_virus$padj),]

# add in human-readable gene names
res_sc2_v_no_virus_ordered$genenames <- gene_annot_tsv[rownames(res_sc2_v_no_virus_ordered), "gene_symbol"]

# write the full table of results to a .csv file
write.csv(res_sc2_v_no_virus_ordered, "results/part3.res_sc2_v_no_virus_ordered.results.csv")

head(res_sc2_v_no_virus_ordered, 10) # show the top differentially expressed genes


# Section 3.2 - Plot a heatmap of the results for comparison of ** SC2 v. no virus **

# Plot Heatmap
library("pheatmap") # load in the pheatmap library, which will be used to generate heatmaps

# select a subset of (the top 25) differential expressed genes to plot
select <- rownames(head(res_sc2_v_no_virus_ordered, 25))

df <- as.data.frame(colData(dds)[,c("viral_status","gender")])
ntd <- normTransform(dds)
heatmap_data <- assay(ntd)[select,]
rownames(heatmap_data) <- gene_annot_tsv[rownames(heatmap_data), "gene_symbol"] # update to human-readable gene names

pdf("results/part3.heatmap.pdf")
pheatmap(heatmap_data, cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df)
dev.off()

# Section 3.3 - Plot a volcano plot of the results for comparison of ** SC2 v. no virus **

# Plot Volcano plot
devtools::install_github('kevinblighe/EnhancedVolcano') # install volcano plot library
library(EnhancedVolcano)

pdf("results/part3.volcano.pdf")
EnhancedVolcano(res_sc2_v_no_virus_ordered,
                lab = res_sc2_v_no_virus_ordered$genenames,
                x = 'log2FoldChange',
                y = 'pvalue',
                pointSize = 3.0,
                labSize = 3.0)
dev.off()

#
## SCRIPT PART 4 - REVIEW THE RESULTS of SC2 v. other virus comparison
#


# SC2 v. other virus
res_sc2_v_other_virus <- results(dds, contrast=c("viral_status","SC2","other_virus"))
res_sc2_v_other_virus_ordered <- res_sc2_v_other_virus[order(res_sc2_v_other_virus$padj),]
res_sc2_v_other_virus_ordered$genenames <- gene_annot_tsv[rownames(res_sc2_v_other_virus_ordered), "gene_symbol"]

# write the results to a .csv file
write.csv(res_sc2_v_other_virus_ordered, "results/part4.res_sc2_v_other_virus_ordered.results.csv")

head(res_sc2_v_other_virus_ordered, 10) # show the top differentially expressed genes


# Plot Heatmap

# select a subset of (the top 25) differential expressed genes to plot
select <- rownames(head(res_sc2_v_other_virus_ordered, 25)) 

df <- as.data.frame(colData(dds)[,c("viral_status","gender")])
ntd <- normTransform(dds)
heatmap_data <- assay(ntd)[select,]
rownames(heatmap_data) <- gene_annot_tsv[rownames(heatmap_data), "gene_symbol"]
pdf("results/part4.heatmap.pdf")
pheatmap(heatmap_data, cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df)
dev.off()

# Plot Volcano plot
pdf("results/part4.volcano.pdf")
EnhancedVolcano(res_sc2_v_other_virus_ordered,
                lab = res_sc2_v_other_virus_ordered$genenames,
                x = 'log2FoldChange',
                y = 'pvalue',
                pointSize = 3.0,
                labSize = 3.0)
dev.off()