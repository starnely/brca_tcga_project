
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
#Converts BMI_group into a categorical variable
  #Required for statistical modeling

clinical_data$BMI_group <- factor(clinical_data$BMI_group) 

table(clinical_data$BMI_group)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")

#limma is a statistical framework designed specifically for gene expression comparisons.
library(limma)

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

#Count significant genes:
sum(results$adj.P.Val < 0.05)

# Check if any
if(nrow(significant_genes) > 0){
  # If gene names are rownames
  significant_gene_names <- rownames(significant_genes)
  message("significant genes ",significant_gene_names)
} else {
  message("No significant genes at adj.P.Val < 0.05")
}

# Subset significant genes
significant_genes <- results[results$adj.P.Val < 0.05, ]

#check the top 20genes genes with strongest evidence even if they did not pass the threshold
head(results[order(results$P.Value), ], 20)

#save results
write.csv(results, "DE_results_obese_vs_nonobese.csv")



#My results are in a data frame 'results' with columns:
  # results$logFC  --> log2 fold change
  # restults$adj.P.Val --> adjusted p-value

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


# For adj.P.Val < 0.05
significant_gene_at_0.05 <- results[results$adj.P.Val < 0.05, ]
significant_gene_at_0.05$Gene <- rownames(significant_gene_at_0.05)
significant_gene_at_0.05

