## Biostatistics 615 Final Project - Sparse Matrix Toolkit for package "Matrix"
### Exploratory Data Analysis of a Matrix object (dense or sparse)
- sparsity: computes sparsity of dense/sparse matrix
- sparse.colMeans: compute colMeans, either with or without zeroes included
- sparse.rowMeans: compute rowMeans, either with or without zeroes included
- sparse.quantile: compute quantile with or without zeroes included(use sparseMatrix slicing w/ "drop=FALSE" to keep the sparse form and prevent memory overflow)
- sparse.summary: similar to the summary() function in base package, compute summary statistics for a sparseMatrix (with or without zeroes) including including sparsity, mean and quartiles. This function does not work for dense matrix.

### Covariance/Correlation of a Matrix object (dense or sparse)
- sparse.covariance
- sparse.correlation (from scratch)
- sparse.correlation2 (sparse.covariance + Matrix::cov2cor(cov), two approachces with close time performance)


### Sparse Linear Regression with optional 2-way interaction
- sparse.lm
