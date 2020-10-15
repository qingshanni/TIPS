# TIPS
A simplified approach for trajectory inference of pathway significance.
# Launch TIPS directly from R and GitHub
## 

### Step 1: Install R and RStudio
Before running TIPS, you will need to have R and RStudio installed

Please check CRAN (https://cran.r-project.org/) for the installation of R.

Please check https://www.rstudio.com/ for the installation of RStudio.

### Step 2: Install required packages

Start an R session using RStudio and run these lines:

- install.packages("shiny")  
- install.packages("shinydashboard")  
- install.packages("ggplot2") 
- install.packages("Seurat")
- install.packages("monocle")  
- install.packages("kohonen")
- install.packages("viridis")
- install.packages("switchde")


### Step 3: Start the app

Start an R session using RStudio and run these lines:

shiny::runGitHub("TIPS", "qingshanni")  
