# Practical Bioinformatics: Assignment 1
Ankit Kumar IIIT Delhi 2021015
April 03, 2023

1. Dataset for this assignment

### Series GSE67255
| Title        | Effects of Systemically Administered Hydrocortisone on the Human Immunome |
|--------------|---------------------------------------------------------------------------|
| Organism     | Homo sapiens                                                              |
| Classes      | 1. dose..250.mg.hydrocortisone <br> 2. dose..50.mg.hydrocortisone         |
| Sample Count | 108                                                                       |
| Link         | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi                            |
| GEO2R        | https://www.ncbi.nlm.nih.gov/geo/geo2r/?acc=GSE67255                      |
|              |                                                                           |


1. Download any microarray data of interest from GEO with at least 100
samples and two classes.
    - Loading GEOquery to directly download the archive of the dataset named "GSE67255_series_matrix.txt.gz" and then extracting using getGEO() function.

    ![image](https://github.com/ankitkat042/BIO221-Practical-Bioinformatics/assets/79627254/095dc92e-1824-4783-84a0-830fea74ff30)

2. After performing EDA and preprocessing, list the data attributes, including pdata and fdata.
![image](https://github.com/ankitkat042/BIO221-Practical-Bioinformatics/assets/79627254/14985d25-1420-4184-8a6b-0e22a28a52c4)

Accessing basic details of dataset and performing EDA

![image](https://github.com/ankitkat042/BIO221-Practical-Bioinformatics/assets/79627254/9875e5e5-5977-445c-bd92-5cc2698959ec)


Boxplot of Data_log:

![image](https://github.com/ankitkat042/BIO221-Practical-Bioinformatics/assets/79627254/3e8890c7-db7f-4f23-a79b-b0e7262ea728)
dimension of data: 33297 x 108

3. State the effects after completing the log transformation of microarray data.
    - Log transformation converts the raw data values into logarithmic scale, which compresses the dynamic range of the data and brings the values closer together. This is particularly useful for data with high variability, such as gene expression data.
    - After log transformation, the data distribution becomes more symmetric and closer to a normal distribution. As a result, logtransformed data is more amenable to statistical analyses such as hypothesis testing and clustering, which assume normality and homogeneity of variance.
    - Overall, the log transformation of microarray data helps to reduce noise, improve accuracy, and increase the sensitivity of downstream analyses.
- This can be seen in the following histograms which represent
variable and distribution of data before and after Log Transform
![image](https://github.com/ankitkat042/BIO221-Practical-Bioinformatics/assets/79627254/d2fa36ab-483c-4435-bccf-7157c9ba492b)

Tough there is not much visible difference but more of the redundant values and negative/NaN  values are removed to make data more consistent.

![image](https://github.com/ankitkat042/BIO221-Practical-Bioinformatics/assets/79627254/8520fd4b-e107-473f-bf47-e5232fd75582)

4. Perform differential expression analysis using simple t-test, log fold change, and correct pvalues using Holm correction. Draw a volcano plot.

![image](https://github.com/ankitkat042/BIO221-Practical-Bioinformatics/assets/79627254/049c9889-c96c-46f3-b2e7-731eea332ef6)

performing the analysis for first and last 3 elements

![image](https://github.com/ankitkat042/BIO221-Practical-Bioinformatics/assets/79627254/0a7e811e-2f7e-47e1-8b19-67a7a882659f)

I performed a t-test on gene expression data for two defined classes, calculated the log fold  change and -log10 p-value for each gene, and created a volcano plot to visualize the results.
The volcano plot shows the log fold change on the x-axis and the negative logarithm of the p-value on the y-axis. The plot_data data frame is created with the log fold change and - log10 p-value for each gene,

5. Perform differential expression analysis using the limma package. Draw a volcano plot.
![image](https://github.com/ankitkat042/BIO221-Practical-Bioinformatics/assets/79627254/d645de1c-9584-4d41-9e72-5b53967f40f2)

Results after differential expression analysis using limma

![image](https://github.com/ankitkat042/BIO221-Practical-Bioinformatics/assets/79627254/cbe837b9-cc7a-44a1-9b04-26b0ff071879)
- I used the limma package to perform differential expression analysis. First, I defined the experimental design using model.matrix(), create the contrast matrix using makeContrasts(), fit a linear model using lmFit(), and performed empirical Bayes moderation using eBayes(). 
- I created a volcano plot using ggplot() and geom_point(). The xaxis represents the log2 fold change, while the y-axis represents the negative log10-transformed p-value. The red dots represent the genes that are significantly differentially expressed between the two groups with a corrected p-value less than 0.05. These genes are considered statistically significant and are often of greater interest in downstream analyses. The black dots, on the other hand, represent genes that are not significantly differentially expressed. The position of a dot in the plot can indicate the magnitude of the gene expression difference and the statistical significance of the difference. Typically, the farther a dot is from the origin, the greater the magnitude of the difference in expression levels, and the higher it is on the y-axis, the more significant the difference.


6. Choose a significant cutoff based on log(FC) and p-values and justify why you chose those values as the cutoff.

The cutoff for logFC is set to 1, which corresponds to a two-fold change in expression between the two classes. For p-value, I used cutoff as 0.05, which corresponds to a false discovery rate of 5%.
**Why?**
This cutoff is biologically relevant because it ensures that only genes with significant changes in expression are considered. A fold change of 2 is also a commonly used threshold in the our PB books.
7) Perform Enrichment analysis using the set of genes that you have obtained using the Gene set enrichment analysis method.
![image](https://github.com/ankitkat042/BIO221-Practical-Bioinformatics/assets/79627254/3b8a24ce-f5e4-4c2f-8d8a-fa306f6281a1)


I performed gene ontology (GO) enrichment analysis using the enrichGO function from the **clusterProfiler** package.   
I extracted the expression data from the GEO dataset and stored it in the 'data' variable. Then divided the data into two groups, 'class1' and 'class2'and performed a t-test on each gene to compare the two groups and calculates the log fold change and -log10 p-value for each gene.   
After performing steps like in step 5 and 6, I Converted the gene symbols to Entrez IDs and performed Gene Ontology (GO) enrichment analysis for biological processes using the enrichGO function.

![image](https://github.com/ankitkat042/BIO221-Practical-Bioinformatics/assets/79627254/85470df3-83c0-4797-88c1-dee058db3b38)
fig : dotplot() of the variable go_enrichment with top 10 pathways  

![image](https://github.com/ankitkat042/BIO221-Practical-Bioinformatics/assets/79627254/2d8a8036-e550-40c2-b669-34485855bff9)  
fig : dotplot() of the variable go_enrichment with top 20 pathways  

![image](https://github.com/ankitkat042/BIO221-Practical-Bioinformatics/assets/79627254/53d1ef76-b2ba-448b-966e-380ebbd1de6c)
fig : barplot() of the variable go_enrichment with all pathways

![image](https://github.com/ankitkat042/BIO221-Practical-Bioinformatics/assets/79627254/c98d1b6b-1702-4552-b8ab-e2611ac703a8)
fig : heatmap() of the variable go_enrichment with 20 pathways

9. Observe and analyze the pathways which you obtained.
![image](https://github.com/ankitkat042/BIO221-Practical-Bioinformatics/assets/79627254/12e84ca4-2201-46cd-a9ce-ebcb84da410e)

- To analyze the pathways obtained from Gene Ontology (GO) enrichment analysis, we first need to understand what these pathways represent. GO is a widely used ontology for describing the biological processes, molecular functions, and cellular components of genes and gene products.

- After performing gene ontology enrichment analysis, we obtained the top 20 enriched pathways. The pathways are ranked based on the p-value obtained from the statistical analysis. The pathways and their descriptions are as follows:

    1. Gland development - refers to the process by which glands, specialized structures that produce and secrete substances, are formed and mature.   
    2. Response to xenobiotic stimulus - refers to the organism's response to exposure to a foreign substance, such as a drug, pollutant, or toxin.   
    3. Positive regulation of cell adhesion - refers to the processes by which cells are encouraged to stick together or to other surfaces, and is important for many cellular processes, including tissue formation and wound healing.   
    4. Response to steroid hormone - refers to the cellular response to steroid hormones, which play a role in a variety of physiological processes such as development, reproduction, and metabolism.   
    5. Regulation of body fluid levels - refers to the maintenance of appropriate levels of fluids in the body, which is important for maintaining homeostasis.   
    6. Response to metal ion - refers to the cellular response to metal ions, which can have toxic effects on the body at high concentrations.   
    7. Muscle system process - refers to the processes involved in the development, function, and regulation of muscle tissue, including muscle contraction and relaxation.   8. Response to peptide hormone - refers to the cellular response to peptide hormones, which play a role in a variety of physiological processes such as growth and development, metabolism, and stress response.

and so on...