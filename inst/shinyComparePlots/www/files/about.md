<!-- TOC depthFrom:1 depthTo:6 withLinks:1 updateOnSave:1 orderedList:0 -->

- [Contact/Feedback](#contactfeedback)
- [Data Sources](#data-sources)
- [Usage](#usage)
- [About the Data](#about-the-data)
- [About CellMinerCDB](#about-cellminercdb)
	- [NCI-DTB Genomics and Bioinformatics Group](#nci-dtb-genomics-and-bioinformatics-group)
	- [Biostatistics and Computational Biology, Dana-Farber Cancer Institute, Harvard Medical School](#biostatistics-and-computational-biology-dana-farber-cancer-institute-harvard-medical-school)
	- [MSKCC Computational Biology](#mskcc-computational-biology)
- [References](#references)

<!-- /TOC -->

# Contact/Feedback<a name="contactfeedback"></a>
Please send comments and feedback to 
* fathi.elloumi AT nih.gov 
* aluna AT jimmy.harvard.edu 
* vinodh.rajapakse AT nih.gov

# Data Sources<a name="data-sources"></a>
CellMinerCDB integrates data from the following sources, which provide additional data and specialized analyses.
* [CellMiner NCI-60](https://discover.nci.nih.gov/cellminer/)
* [Sanger/Massachusetts General Hospital Genomics of Drug Sensitivity in Cancer (GDSC)](http://www.cancerrxgene.org/)
* [Broad/Novartis Cancer Cell Line Encyclopedia (CCLE)](https://portals.broadinstitute.org/ccle)
* [Broad Cancer Therapeutics Response Portal (CTRP)](https://portals.broadinstitute.org/ctrp/)
* [NCI/DTP Small Cell Lung Cancer Project (SCLC)](https://sclccelllines.cancer.gov/sclc/)

# Usage<a name="usage"></a>
**For detailed instructions, please refer to our [video tutorial](https://youtu.be/WJ_A_kzxO-Y).**

CellMinerCDB simplifies exploration of cancer cell line pharmacogenomic data from the above sources. Users can plot and more
broadly interrelate molecular and drug activity profiling data, both within and across data sources.

The example shown below is a scatterplot of PAXX (C9orf142) gene expression and Bleomycin (NSC294979) activity within the NCI60.
For the NCI-60 cell lines, tissues of origin are indicated using an established color scheme. Two different data sources can be 
specified for the x and y axis variables. In this case, where data are available, plots and other analyses will be restricted to
overlapping cell lines between the two data sources.

From the **'Univariate Analyses'** navbar tab at the top left of the application (shown below), users can
* plot observations for any pair of cell line data variables ('Plot Ids' tab)
* view, filter, and download the plot data in tabular, Excel-readable form ('Download Data' tab)
* search for identifiers (genes, drugs, etc.) available within the x-axis variable-associated data source ('Search Ids' tab)
* identify and rank additional variables correlated with either the x-axis or y-axis variables ('Compare Patterns' tab)

From the **'Regression Models'** navbar tab at the top left, users can incrementally develop multivariate response 
prediction models and
* view data for the most (least) responsive cell lines in an interactive heatmap ('Heatmap' tab)
* view, filter, and downloaded the modeling data in tabular, Excel-readable form ('Download Data' tab)
* view plots of the observed response values versus the predicted response values ('Plot' tab)
* view plots of the observed response values versus the 10-fold cross-validation predicted response values ('Cross-Validation' tab)
* view statistical and other technical details relating to the constructed response model ('Technical Details' tab)
* identify additional variables that could improve the existing response prediction model ('Partial Correlation' tab)



![Screenshot of CellMinerCDB Application](files/rcellminer_screenshot_anno.png)

# About the Data<a name="about-the-data"></a>
For specific information about the data made available for particular sources, please refer to the 'Metadata' navbar tab.

Drug mechanism of action details:
* [NCI60](https://raw.githubusercontent.com/cannin/rcellminer/devel/inst/extdata/Drug_MOA_Key.txt)
* [GDSC](http://www.cancerrxgene.org/translation/Drug)
* [CTRP](https://portals.broadinstitute.org/ctrp/?page=#ctd2Compounds)

Gene sets used for annotation of analysis results or algorithm input filtering were curated by the
NCI/DTB CellMiner team, based on surveys of the applicable research literature.

# About CellMinerCDB<a name="about-cellminercdb"></a>
The CellMinerCDB application is developed and maintained using R and Shiny by:

* Augustin Luna; Research Fellow, Biostatistics and Computational Biology, Dana-Farber Cancer Institute, Harvard Medical School
* Vinodh N. Rajapakse; Postdoctoral Fellow, Developmental Therapeutics Branch, National Cancer Institute
* Fathi Elloumi; Bioinformatics Software Engineer, Developmental Therapeutics Branch, National Cancer Institute

### NCI-DTB Genomics and Bioinformatics Group<a name="nci-dtb-genomics-and-bioinformatics-group"></a>
* William C. Reinhold
* Sudhir Varma
* Margot Sunshine
* Fathi Elloumi
* Lisa Loman (Special Volunteer)
* Fabricio G. Sousa
* Kurt W. Kohn
* Yves Pommier

### Biostatistics and Computational Biology, Dana-Farber Cancer Institute, Harvard Medical School<a name="biostatistics-and-computational-biology-dana-farber-cancer-institute-harvard-medical-school"></a>
* Chris Sander

### MSKCC Computational Biology<a name="mskcc-computational-biology"></a>
* Jianjiong Gao
* Nikolaus Schultz

# References<a name="references"></a>
Shankavaram UT, Varma S, Kane D, Sunshine M, Chary KK, Reinhold WC, Pommier Y, Weinstein JN. [CellMiner: a relational database and query tool for the NCI-60 cancer cell lines. ](https://www.ncbi.nlm.nih.gov/pubmed/19549304) BMC Genomics. 2009 Jun 23;10:277. doi: 10.1186/1471-2164-10-277.

Reinhold WC, Sunshine M, Liu H, Varma S, Kohn KW, Morris J, Doroshow J, Pommier Y.
[CellMiner: a web-based suite of genomic and pharmacologic tools to explore transcript and drug patterns in the NCI-60 cell line set. ](https://www.ncbi.nlm.nih.gov/pubmed/22802077) Cancer Res. 2012 Jul 15;72(14):3499-511. doi: 10.1158/0008-5472.CAN-12-1370.

Reinhold WC, Sunshine M, Varma S, Doroshow JH, Pommier Y. [Using CellMiner 1.6 for Systems Pharmacology and Genomic Analysis of the NCI-60. ](https://www.ncbi.nlm.nih.gov/pubmed/26048278) Clin Cancer Res. 2015 Sep 1;21(17):3841-52. doi: 10.1158/1078-0432.CCR-15-0335. Epub 2015 Jun 5.

Luna A, Rajapakse VN, Sousa FG, Gao J, Schultz N, Varma S, Reinhold W, Sander C, Pommier Y. [rcellminer: exploring molecular profiles and drug response of the NCI-60 cell lines in R.](https://www.ncbi.nlm.nih.gov/pubmed/26635141) Bioinformatics. 2015 Dec 3. pii: btv701.
