import pandas as pd #tables
import numpy as np #math
import matplotlib.pyplot as plt #graphs

# 1. Data loading
df = pd.read_csv(
    r"C:\Users\hanaf\Desktop\projekt AD\variant_summary.txt.gz",
    sep="\t", #columns are separated by tabs
    compression="gzip", #pandas can automaticaly decompress the file
    low_memory=False #can read whole file at once
)

print(df.columns) #name of columns
print(df.head()) #first 5 rows as control

# 2. Filtering of AD data
alz = df[df["PhenotypeList"].str.contains("Alzheimer", na=False)] #filtering rows where "PhenotypeList" contains "Alzheimer", na=False means that if there is no data in "PhenotypeList" it will be treated as False and not included in the result

# 3. Quick check
print(alz.shape)#print how many columns and how many rows
print(alz.head())# will show first 5
print(alz["GeneSymbol"].value_counts().head(10)) #number of variants per genes, top 10 genes

# 4. Visualisation - Top Alzheimer genes
#The graph then shows the number of variants in individual genes that are associated with Alzheimer's in the ClinVar database. (Top 10 most represented genes)
top_genes = alz["GeneSymbol"].value_counts().head(10) #count how many times each gene appears in the "GeneSymbol" column, sort them in descending order and take the top 10
print(top_genes)

#this will mage graph
plt.figure(figsize=(10,5)) #10x5 inches
plt.bar(top_genes.index, top_genes.values, width=0.6)
plt.title("Top Alzheimer-associated genes (ClinVar)\n TOP 10 genes") #graph title
plt.ylabel("Number of reported variants") # y axis name
plt.xlabel("Gene") #x axis name
plt.xticks(rotation=45) #rotate the names of gene
plt.tight_layout() #Automatically adjusts the margins
plt.savefig(r"C:\Users\Desktop\projekt AD\alzheimers_top_genes.png", dpi=300)
plt.show() #show the graph

print(alz.shape)
print(alz["ClinicalSignificance"].value_counts())  # signifikance
alz.to_csv(r"C:\Users\Desktop\projekt AD\alzheimers_subset.csv", index=False) #save the Ad subset
df.to_csv(r"C:\Users\Desktop\projekt AD\clinvar_all.csv", index=False) #save the whole dataset


# Risk score (weighted)
gene_counts = alz["GeneSymbol"].value_counts() #variant count

risk = pd.DataFrame(gene_counts).reset_index() #transfer to dataframe
risk.columns = ["Gene", "VariantCount"] #rename

# Base
risk["RiskScore"] = risk["VariantCount"] # risk score == variant count

# Weighted score
weight_map = {
    "Pathogenic": 3,
    "Likely pathogenic": 2,
    "Uncertain significance": 1,
    "Benign": 0
}

alz["Weight"] = alz["ClinicalSignificance"].map(weight_map).fillna(0) #formula, N/A is change to 0

weighted = alz.groupby("GeneSymbol")["Weight"].sum().reset_index()#sum of score, grouped by gene symbol
weighted.columns = ["Gene", "WeightedScore"]

# Merge together
risk = risk.merge(weighted, on="Gene", how="left") #left join

risk["WeightedScore"] = risk["WeightedScore"].fillna(0)

print(risk.sort_values("WeightedScore", ascending=False).head(10))

top = risk.sort_values("WeightedScore", ascending=False).head(10) #arrangement, ascending, TOP 10

#Graph for weighted score
plt.figure(figsize=(10,5))
plt.bar(top["Gene"], top["WeightedScore"], width=0.6)
plt.title("Weighted Alzheimer Gene Risk Score\n TOP 10 genes")
plt.xticks(rotation=45)
plt.tight_layout()

plt.savefig(r"C:\Users\Desktop\projekt AD\weighted_risk_score.png", dpi=300)
plt.show()

# Normalized score
risk["NormRisk"] = risk["WeightedScore"] / risk["WeightedScore"].max() #normalisation to max

print(risk.sort_values("NormRisk", ascending=False).head(10))

top = risk.sort_values("NormRisk", ascending=False).head(10) #arrangement, ascending top 10

#Graph for normalized values
plt.figure(figsize=(10,5))
plt.bar(top["Gene"], top["NormRisk"], width=0.6)
plt.title("Normalized Alzheimer Gene Risk Score \n TOP 10 genes")
plt.xticks(rotation=45)
plt.tight_layout()

plt.savefig(r"C:\Users\hanaf\Desktop\projekt AD\weighted_risk_score_normalized.png", dpi=300)
plt.show()

