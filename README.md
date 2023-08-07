# Differential Expression RNA-seq Analysis

**DEseq2**

**Principal component analysis (PCA) analysis**

**KEGG pathway enrichment analysis**

## Project Overview:

Investigating Key Pathway Involving Mcl-1 Proteins in Mammary Glands

To explore the significant pathway associated with Mcl-1 proteins in mammary glands, we utilized publicly available RNA-Seq data from luminal (Wild-type) and basal cells (Mcl-1 deletion) obtained from virgin, 18.5-day pregnant, and 2-day lactating mice. The original study was documented in (PMID: 25730472), and the corresponding RNA-Seq data can be accessed at GSE60450. Our analysis involved several steps, including Differential **Expression Analysis** using **DEseq2**, **Principal Component Analysis (PCA)**, and **KEGG pathway enrichment analysis**, followed by visualization of the results.

-   **THE NUMBER OF DIFFERENTIAL EXPRESS GENES.**

    ![](https://github.com/chingyaousf/Differential-Expression-RNA-seq-Analysis/blob/main/plots/Bar_plot_differential_expressed_genes.png?raw=true){width="374"}

-   **GENES ARE CONSISTENTLY DIFFERENTIALLY EXPRESSED IN ALL THREE STAGES BETWEEN LUMINAL (WILD-TYPE) AND BASAL CELLS (MCL-1 DELETION).**

![](https://github.com/chingyaousf/Differential-Expression-RNA-seq-Analysis/blob/main/plots/Venn_diagram_upregulated_downregulated_genes.png?raw=true){width="639"}

![](https://github.com/chingyaousf/R-Differential-Expression-RNA-seq-Analysis/blob/main/plots/Venn_diagram_upregulated_downregulated_genes_02.png?raw=true)

-   **LUMINAL CELL DURING THE LACTATION STAGE IS FAR FROM THE REMAINING SAMPLES IN THE PCA PLOT.**![](https://github.com/chingyaousf/Differential-Expression-RNA-seq-Analysis/blob/main/plots/PCA.png?raw=true){width="432"}

-   **THE VOLCANO PLOT FOR THREE STAGES.** ![](https://github.com/chingyaousf/Differential-Expression-RNA-seq-Analysis/blob/main/plots/Volcano_virgin_pregnant_lactation.png?raw=true)

-   **THE BAR PLOT REGARDING KEGG ANALYSIS FOR THREE STAGES.** ![](https://github.com/chingyaousf/Differential-Expression-RNA-seq-Analysis/blob/main/plots/KEGG_annotation_virgin_stage.png?raw=true){width="450"}![](https://github.com/chingyaousf/Differential-Expression-RNA-seq-Analysis/blob/main/plots/KEGG_annotation_pregnant_stage.png?raw=true){width="440"}

![](https://github.com/chingyaousf/Differential-Expression-RNA-seq-Analysis/blob/main/plots/KEGG_annotation_latate_stage.png?raw=true){width="399"}

## Discussion:

•Trend in the Volcano plot and Venn diagrams displaying differentially expressed genes are consistent with Western blot analysis, immunohistochemistry and qRT--PCR analysis = more upregulation during lactation stage.  

•EGF seems to play a major role in orchestrating the increased translation of Mcl-1 through the mTOR pathway.

•The results of the differential gene expression analysis and subsequent pathway analysis complement the authors findings by highlighting the role of the mTOR pathways and the biological processes its genes are involved with. 

**Limitations and future direction:**

•Only three mice per stage (small sample size). 

Future studies could include a larger sample size.

## References:

Fu, N. Y., Rios, A. C., Pal, B., Soetanto, R., Lun, A. T. L., Liu, K., Beck, T., Best, S. A., Vaillant, F., Bouillet, P., Strasser, A., Preiss, T., Lindeman, G. J., & Visvader, J. E. (2015). EGF-mediated induction of Mcl-1 at the switch to lactation is essential for alveolar cell survival. *Nature Cell Biology*, *17*(4), 365--375. <https://doi.org/10.1038/ncb3117>

Hinck, L. & Silberstein, G. B. Key stages in mammary gland development: the mammary end bud as a motile organ. Breast Cancer Res. 7, 245--251 (2005). 

Watson, C. J. Involution: apoptosis and tissue remodelling that convert the mammary gland from milk factory to a quiescent organ. Breast Cancer Res. 8, 203 (2006).

## Blog:

<https://ssidmarine.wordpress.com/2023/08/07/differential-expression-rna-seq-analysis/>

## Access data:

<https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60450>

PMID: 25730472

GSE60450_Lactation-GenewiseCounts

kegg.pathway.gene

SampleInfo_Corrected

**Data available in the data folder**
