# TIPS
Trajectory Inference of Pathway Significance through Pseudotime Comparison for Functional Assessment of
single-cell RNAseq Data
## Launch TIPS
### Step 1: Install R and RStudio
Before running TIPS, you will need to have R and RStudio installed
Please check CRAN (https://cran.r-project.org/) for the installation of R.
Please check https://www.rstudio.com/ for the installation of RStudio.

### Step 2: Install required packages

Start an R session using RStudio and run these lines:
```
install.packages(c("shiny","shinydashboard","markdown","ggplot2","Seurat","kohonen","viridis"))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("monocle","switchde"))
```

**Notes**

*TIPS was written and tested with the following specific versions of packages:*

    ggplot2 3.3.0
    kohonen 3.0.10
    markdown 1.1
    monocle 2.14.0
    Seurat 2.3.4
    shiny 1.4.0.2
    shinydashboard 0.7.1
    switchde 1.12.0
    viridis 0.5.1

We recommend installing these versions to ensure compatibility with TIPS.



### Step 3: Start the app

Start an R session using RStudio and run these lines:

```
library(shiny)
shiny::runGitHub("TIPS", "qingshanni")    
```

Or you can clone or download this repository, and run:

```
shiny::runApp("TIPS")
```

## Citation
Please use the following citation:

Zheng Z, Qiu X, Wu H, Chang L, Tang X, Zou L, Li J, Wu Y, Zhou J, Jiang S, Wan Y, Ni Q. TIPS: trajectory inference of pathway significance through pseudotime comparison for functional assessment of single-cell RNAseq data. Brief Bioinform. 2021 Sep 2;22(5):bbab124. doi: [10.1093/bib/bbab124](https://doi.org/10.1093/bib/bbab124). PMID: 34370020; PMCID: PMC8425418.


## License


This program is free software and can be redistributed and/or modified under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the license, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but without any real or implied warranty of merchantability or fitness for a particular purpose. See the GNU General Public License for more details.
