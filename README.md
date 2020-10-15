# TIPS
A simplified approach for trajectory inference of pathway significance.
# Launch TIPS directly from R and GitHub
### Step 1: Install R and RStudio
Before running TIPS, you will need to have R and RStudio installed

Please check CRAN (https://cran.r-project.org/) for the installation of R.

Please check https://www.rstudio.com/ for the installation of RStudio.

### Step 2: Install required packages

Start an R session using RStudio and run these lines:
```
install.packages(c("shiny","shinydashboard","ggplot2","Seurat","kohonen","viridis"))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("monocle","switchde"))
```

### Step 3: Start the app

Start an R session using RStudio and run these lines:
```
shiny::runGitHub("TIPS", "qingshanni")    
```
# Documentation
Detailed usage instructions can be found in the user manual 

# License
Copyright(c) <2020><ZH Zheng, Y Wan, QS Ni, AMU China, All Rights Reserved

This program is free software and can be redistributed and/or modified under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the license, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but without any real or implied warranty of merchantability or fitness for a particular purpose. See the GNU General Public License for more details.