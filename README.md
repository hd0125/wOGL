# wOGL
An R package: "Weighted Overlapping Group Lasso" for the selection of gene sets as groups

# Description
Gene set analysis aims to identify differentially expressed gene sets, often ignoring genetic network structure, which is less effective for sparse signals. Weighted Overlapping Group Lasso leverages network knowledge to identify interconnected genes, combining network-based regularization with overlapping group lasso, using the l2-norm of regression coefficients for individual genes, which can play a role in the weight of gene sets for group selection.

# Installation
To install the 'devtools' package in R, you can use the following steps:
1. Install 'devtools' package:
```
install.packages('devtools')
```
2. For Windows users, it's also necessary to install Rtools from the following URL:  
    - https://cran.r-project.org/bin/windows/Rtools
4. Once 'devtools' and Rtools are installed, you can load the 'devtools' library and use it to install packages from GitHub:
```
library(devtools)
install_github("statddd25/wOGL")
```

# Reference
- **Sun, H.** and Wang, S. (2012) Penalized Logistic Regression for High-dimensional DNA Methylation Data with Case-Control Studies, Bioinformatics 28(10), p.1368-1375.
- **Sun, H.** and Wang, S. (2013) Network-based Regularization for Matched Case-Control Analysis of High-dimensional DNA Methylation Data, Statistics in Medicine 32(12), p.2127-2139.
- **Yang Yi** and Hui Zou (2015) A fast unified algorithm for solving group-lasso penalize learning problems, Statistics and Computing 25, 1129-1141.

# Report bugs
- Send an email to Dan Huang at dhuang1221@gmail.com
