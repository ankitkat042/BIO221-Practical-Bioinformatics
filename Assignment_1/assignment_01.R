# Load required libraries
library(ggplot2)
library(GEOquery)
library(limma)
library(tidyverse)
library(impute)
library(tibble)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)


###  1 and 2  ###
gse <- getGEO("GSE67255")
gse <- gse[[1]]
exprs  <- exprs(gse)
pData <- pData(gse)
fData <- fData(gse)


# Extract the expression data
data <- exprs(gse)

# Step 1: Check the dimensions of the dataset
dim(data)

# Step 2: Check the structure of the dataset
str(data)

# Step 3: Check summary statistics of the dataset
summary(data)

# Step 4: Check for missing values in the dataset
sum(is.na(data))

# if negative values deal with them
min_value <- min(data, na.rm = TRUE)
data_shifted <- data - min_value + 1

# Step 5: Preprocessing - Log transformation

data_log <- log2(data_shifted)
boxplot(data_log)

##### 3. State the effects after completing the log transformation of microarray data. #### 
par(mfrow = c(1, 2))
hist(data[1, ], main = "Before Log Transformation", xlab = "Expression Level")
hist(data_log[1, ], main = "After Log Transformation", xlab = "Expression Level")

##### 4 #####

# Load data using GEOquery

data <- exprs(gse)

# Define classes
class1 <- data[,1:3]
class2 <- data[,4:6]
print(class1)
print(class2)

# Perform t-test
ttest <- apply(data, 1, function(x) t.test(x[1:3], x[4:6])$p.value)

# Calculate log fold change and -log10 p-value
logFC <- apply(data, 1, function(x) log2(mean(x[1:3]) / mean(x[4:6])))
logP <- -log10(ttest)

# Create data frame for plotting
plot_data <- data.frame(logFC = logFC, logP = logP)

# Create volcano plot
ggplot(plot_data, aes(x = logFC, y = logP)) +
  geom_point() +
  xlab("Log2 Fold Change") +
  ylab("-Log10 p-value")

##### 5 #####

colnames(pData(gse))
#number of classes in ch1.1
new_levels <- make.names(levels(pData(gse)$characteristics_ch1.1))

pData(gse)$characteristics_ch1.1 <- factor(pData(gse)$characteristics_ch1.1)
Group1 <- levels(pData(gse)$characteristics_ch1.1)[1]
Group2 <- levels(pData(gse)$characteristics_ch1.1)[2]

design <- model.matrix(~ 0 + pData(gse)$characteristics_ch1.1)
colnames(design) <- new_levels
cont.matrix <- makeContrasts(dose..250.mg.hydrocortisone - dose..50.mg.hydrocortisone, levels = design)
gg <- exprs(gse)
fit <- lmFit(gg, design)

#  Calculate the moderated t-statistics:
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, trend = TRUE)



#########################
fit <- lmFit(data, design)
fit <- contrasts.fit(fit, cont.matrix)
fit <- eBayes(fit)

# Extract the log-fold change and p-value
results <- topTable(fit, adjust.method = "holm", sort.by = "none", number = Inf)
logFC <- results$logFC
PValue <- results$P.Value

## Create a volcano plot
plot_data <- data.frame(logFC, -log10(PValue))
ggplot(plot_data, aes(x = logFC, y = -log10(PValue))) +
  geom_point(aes(color = ifelse(PValue < 0.05, "red", "black"))) +
  scale_color_identity() +
  ggtitle("Volcano Plot") +
  xlab("Log Fold Change") +
 ylab("-Log10 P-Value")


#### 6 ####

results <- rownames_to_column(results, "Gene_ID")

# select significant genes based on fold change and p-value cutoffs
de_genes <- results %>%
  filter(abs(logFC) >= 1, P.Value <= 0.05) %>%
  select(Gene_ID)

#### 7 ####
results <- rownames_to_column(results, "Gene_ID")
de_genes <- results %>%
  filter(abs(logFC) >= 1, PValue <= 0.05) %>%
  dplyr::select(Gene_ID)

de_genes_entrez <- unique(results$Gene_ID)


de_genes_entrez <- as.list(factor(de_genes_entrez))
class(de_genes_entrez)

# Perform enrichment analysis using GO Biological Process
go_enrichment <- enrichGO(de_genes_entrez$ENTREZID, OrgDb="org.Hs.eg.db", keyType="ENTREZID", ont="BP", pvalueCutoff=0.05, qvalueCutoff=0.1)

de_genes_entrez <- as.list(factor(de_genes_entrez))
go_enrichment <- enrichGO(de_genes_entrez, OrgDb="org.Hs.eg.db", keyType="ENTREZID", ont="BP", pvalueCutoff=0.05, qvalueCutoff=0.1)
dotplot(go_enrichment)
dotplot(go_enrichment, showCategory=20)

barplot(go_enrichment$Count, names.arg=go_enrichment$Description,
        ylab="Count", xlab="GO Term Description", las=2,
        main="Enriched GO Terms")



# View top enriched terms
head(go_enrichment)



