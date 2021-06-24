#Carregar as bibliotecas instaladas. Deve ser executado no inicio de toda a análise
library("IHW")
library("DESeq2")
library("ggplot2")
library("RColorBrewer")
library("pheatmap")
library("tximport")
#library("biomart")

print(snakemake@output[[1]])

#Carrega os nomes dos arquivos de expressão de cada biblioteca
tx2gene <- read.csv(snakemake@params[["transcripts_to_genes"]], header = FALSE, sep=",")
head(tx2gene)
print("\n\n")
head(snakemake@params[["samples_file"]])
print("\n\n")
samples_table <- read.csv(snakemake@params[["samples_file"]], header = TRUE, sep=",")
samples <- samples_table$sample
files <- file.path("kallisto",samples, "abundance.tsv")
names(files) <- samples
files
print("\n\n")
#Importa os resultados de expressão de cada biblioteca
txi.kallisto.tsv <- tximport(files, type = "kallisto",  tx2gene = tx2gene, dropInfReps = TRUE)
head(samples)
head(txi.kallisto.tsv$counts)
colnames(txi.kallisto.tsv$counts) <- samples
head(txi.kallisto.tsv$counts)

colDatak <-data.frame(row.names=colnames(txi.kallisto.tsv$counts), condition=as.factor(samples_table$condition), type=as.factor(samples_table$type))
colDatak

ddsk <- DESeqDataSetFromTximport(txi = txi.kallisto.tsv, colData = colDatak, design = ~ type + condition + type:condition)

#Cria objeto do DESeq2
ddsk <- DESeq(ddsk)

print("Criar gráfico do Heatmap")
#Cria gráfico do Heatmap
vsdk <- rlog(ddsk, blind=FALSE)
sampleDistsk <- dist(t(assay(vsdk)))
sampleDistMatrixk <- as.matrix(sampleDistsk)
rownames(sampleDistMatrixk) <- paste(vsdk$type, vsdk$condition, sep="-")
#colnames(sampleDistMatrixk)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
print("Heatmap")
png("Heatmap.png")
pheatmap(sampleDistMatrixk,
         clustering_distance_rows=sampleDistsk,
         clustering_distance_cols=sampleDistsk,
         col=colors)
dev.off()

#Cria PCA
pcaDatak <- plotPCA(vsdk, intgroup=c("type", "condition"), returnData=TRUE)
percentVark <- round(100 * attr(pcaDatak, "percentVar"))
png("PCA.png")
ggplot(pcaDatak, aes(PC1, PC2, color=condition, shape=type)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVark[1],"% variance")) +
  ylab(paste0("PC2: ",percentVark[2],"% variance")) +
  coord_fixed()
dev.off()
print("Results")

ddsk$group <- factor(paste0(ddsk$type, "_", ddsk$condition))
design(ddsk) <- ~ group
ddsk <- DESeq(ddsk)

#Comparações que serão feitas
res_ribo_S0_S14<-results(ddsk, contrast=c("group","ribo_D14","ribo_NPC"),  filterFun=ihw)
res_rna_S0_S14<-results(ddsk, contrast=c("group","rna_D14","rna_NPC"),  filterFun=ihw)
res_ribo_S0_rna_S0<-results(ddsk, contrast=c("group","ribo_NPC","rna_NPC"),  filterFun=ihw)
res_ribo_S14_rna_S14<-results(ddsk, contrast=c("group","ribo_D14","rna_D14"),  filterFun=ihw)

#Adiciona anotação para os genes
split_gene_id <- function(word){
    strsplit(word, ".", fixed = TRUE)[[1]][1]
    #unlist(strsplit(my.string, ".", fixed = TRUE))[1]
}

gene_id = apply(as.matrix(rownames(res)), 1, split_gene_id)

rownames(res) <- gene_id

res$ensembl_id <- rownames(res)
martk <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
genes.tablek <- biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description", "gene_biotype"), mart = martk)
res <- merge(x = as.matrix(res), y = genes.tablek, by.x = "ensembl_id", by.y = "ensembl_gene_id", all.x = T, all.y = F )
counts_normalizedk<-counts(ddsk, normalized=TRUE)
counts_rawk<-counts(ddsk)
Tabelao<-cbind(res,counts_rawk,counts_normalizedk)

write.table(as.data.frame(Tabelao),file=snakemake@output[[1]], sep = "\t")
