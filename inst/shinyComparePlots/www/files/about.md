<!-- TOC depthFrom:1 depthTo:6 withLinks:1 updateOnSave:1 orderedList:0 -->

- [Contact/Feedback](#contactfeedback)
- [Usage](#usage)
- [About the Data](#about-the-data)
	- [Mechanisms of Action](#mechanisms-of-action)
	- [Gene Sets in Regression Tab](#gene-sets-in-regression-tab)
	- [Gene Annotations](#gene-annotations)
- [About CellMinerCDB](#about-cellminercdb)
	- [MSKCC Computational Biology](#mskcc-computational-biology)
	- [Biostatistics and Computational Biology, Dana-Farber Cancer Institute, Harvard Medical School](#biostatistics-and-computational-biology-dana-farber-cancer-institute-harvard-medical-school)
	- [NCI-DTB Genomics and Bioinformatics Group](#nci-dtb-genomics-and-bioinformatics-group)
- [References](#references)

<!-- /TOC -->

# Contact/Feedback

Please send comments and feedback to aluna AT jimmy.harvard.edu and vinodh.rajapakse AT nih.gov

# Usage
This application simplifies the exploration of datasets connected to the NCI-60. Users are allows to view both molecular data and drug activity for cell lines. The example shows a scatterplot of PAXX (C9orf142) gene expression and Bleomycin (NSC294979) activity within the NCI60. In the application, point colors represent tissues of origin for the cell lines.

![Screenshot of CellMinerCDB Application](files/rcellminer_screenshot_anno.png)

# About the Data

Please see the Metadata tab for information data.

## Mechanisms of Action

* NCI60: [Mechanism of Action Description](https://raw.githubusercontent.com/cannin/rcellminer/devel/inst/extdata/Drug_MOA_Key.txt)

## Gene Sets in Regression Tab

Gene sets expertly curated by the CellMiner team.

## Gene Annotations

Gene annotations are based on membership in CellMiner gene sets or other literature-based annotations.

# About CellMinerCDB
The CellMinerCDB application was developed using R and Shiny by:

* Augustin Luna; Research Fellow, Biostatistics and Computational Biology, Dana-Farber Cancer Institute, Harvard Medical School
* Vinodh N. Rajapakse; Postdoctoral Fellow, Developmental Therapeutics Branch, National Cancer Institute
* Lisa Loman; Special Volunteer, Developmental Therapeutics Branch, National Cancer Institute
* Harini Salgado; Special Volunteer, Developmental Therapeutics Branch, National Cancer Institute

## MSKCC Computational Biology
* Jianjiong Gao
* Nikolaus Schultz

## Biostatistics and Computational Biology, Dana-Farber Cancer Institute, Harvard Medical School
* Chris Sander

## NCI-DTB Genomics and Bioinformatics Group
* William C. Reinhold
* Sudhir Varma
* Margot Sunshine
* Fabricio G. Sousa
* Kurt W. Kohn
* Yves Pommier

# References
Luna A, Rajapakse VN, Sousa FG, Gao J, Schultz N, Varma S, Reinhold W, Sander C, Pommier Y. [rcellminer: exploring molecular profiles and drug response of the NCI-60 cell lines in R.](https://www.ncbi.nlm.nih.gov/pubmed/26635141) Bioinformatics. 2015 Dec 3. pii: btv701.
