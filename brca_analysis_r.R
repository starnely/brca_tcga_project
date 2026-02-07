rm(list = ls())

# Load libraries
library(survival)
library(survminer)
library(GEOquery)


# Read expression data
expression_data <- read.csv("expression_data_for_R.csv", row.names = 1)

# Read clinical data
clinical_data <- read.csv("clinical_data_for_R.csv")

#Am renaming my own data to avoid duplicate since when you edit and run the code m R will create a new object at every runtime
#Rename expression_data from expr
#expression_data <- expr
#rm(expr)
#Rename clinical_data from clinical
#clinical_data <- clinical
#rm(clinical)


# Checking the shape of my two tables in R
dim(expression_data)
dim(clinical_data)

#This check whether the columns on Expression data and row on the clinical data are same and in the same order.
#This is because Survival analysis assumes:gene expression values and clinical info belong to the same patients
all(colnames(expression_data) == clinical_data$SampleID)

#Tells R which samples are Obese and non-Obese
#Converts BMI_group into a categorical variable required for statistical modeling
clinical_data$BMI_group <- factor(clinical_data$BMI_group) 

#returns a table with the number of obese patients and non-obese patients
table(clinical_data$BMI_group)


#Installing Limma for gene expression
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
  BiocManager::install("limma")
  library(limma)#Limma is a statistical framework designed specifically for gene expression comparisons.

#Create the design matrix (this is the heart of the analysis)
#Compare gene expression between obese and non-obese samples.”
design <- model.matrix(~ BMI_group, data = clinical_data)

(design) #checking how it looks

head(design) ##checking how it looks(sample)

#Fit the linear model (gene-by-gene)
  #Fits a model for each gene
  #Estimates expression difference between groups
  #Stabilizes variance (critical for human data)
fit <- lmFit(expression_data, design)
fit <- eBayes(fit)

#| #Column    | Meaning (simple)           |
 # | --------- | -------------------------- |
  #| #logFC     | Direction & size of change |
  #| #P.Value   | Raw p-value                |
  #| #adj.P.Val | FDR-corrected p-value      |
  #| #t         | Strength of evidence       |
  
results <- topTable(
  fit,
  coef = "BMI_groupObese",
  number = Inf,
  adjust.method = "fdr"
)

(results)

#Count significant genes using an adjusted p value of 0.05:
sum(results$adj.P.Val < 0.05)

#Results
  #Only 1 gene at 0.05
  #Many more at 0.1
  #This suggests:
    #Obesity has a subtle effect on tumor gene expression
    #Only one gene shows a very strong, confident change

# From the results table, keep only the genes whose adjusted p-value is less than 0.05
significant_gene_at_0.05 <- results[results$adj.P.Val < 0.05, ]

# Check if any
if(nrow(significant_gene_at_0.05) > 0){
  # If gene names are row names
  significant_gene_name <- rownames(significant_gene_at_0.05)
  message("significant genes at pvalue 0.05 is ",significant_gene_name)
} else {
  message("No significant genes at adj.P.Val < 0.05")
}


#I Decided change the adjusted p value to 0.1

#Count significant genes:
sum(results$adj.P.Val < 0.1)

#Results
  #87 gene at 0.1

# From the results table, keep only the genes whose adjusted p-value is less than 0.1
significant_genes_at_0.1 <- results[results$adj.P.Val < 0.1, ]

#look at the name  of the genes them
if(nrow(significant_genes_at_0.1) > 0){
  # If gene names are row names
  significant_gene_names <- rownames(significant_genes_at_0.1)
  message("significant genes at pvalue 0.1 is ",significant_gene_names)
} else {
  message("No significant genes at adj.P.Val < 0.1")
}

#check the top 20genes genes with strongest evidence even if they did not pass the threshold
head(results[order(results$P.Value), ], 20)

#save results
write.csv(results, "Defferential_results_obese_vs_nonobese.csv")

#installing ggplot2
library(ggplot2)

# Add a column to label significant genes at 0.05
results$Significant_gene_at_0.05 <- ifelse(results$adj.P.Val < 0.05, "Yes", "No")

# Volcano plot for 0.05 threshold
ggplot(results, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = Significant_gene_at_0.05), alpha = 0.6) +
  scale_color_manual(values = c("No" = "grey", "Yes" = "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot (adj.P.Val < 0.05)",
       x = "log2 Fold Change", y = "-log10(adj.P.Val)")

#Results explains that
  #Obesity has a focused, not widespread impact on gene expression in BRCA
  #This means Obesity does not change many genes in breast cancer tumors.
  #Instead, it affects only a small number of specific genes.



head(clinical_data)
colnames(clinical_data)









#Install & load packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("TCGAbiolinks", "SummarizedExperiment"))
install.packages(c("dplyr", "ggplot2"))

library(TCGAbiolinks) #Package to download and prepare TCGA data
library(SummarizedExperiment)
library(dplyr) #Data handling
library(ggplot2) # Visualization

#Choose TCGA project
project <- "TCGA-BRCA"

#Download TCGA RNA-seq data
query <- GDCquery(
  project = project,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

#Download the data to your computer
GDCdownload(query) 

#Prepare data as a SummarizedExperiment object
exp.data <- GDCprepare(query)

clinical <- GDCquery_clinic(project, type = "clinical")

clinical <- clinical %>%
  select(bcr_patient_barcode, bmi)

counts <- assay(exp.data)
gene.info <- rowData(exp.data)

# Convert Ensembl IDs → gene symbols
gene.symbols <- gene.info$external_gene_name
rownames(counts) <- gene.symbols

gene <- "ADI1"

gene.exp <- counts[gene, ]
