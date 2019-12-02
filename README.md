## Biostatistics 615 Final Project - Sparse Matrix Statistics and Modeling
### Exploratory Data Analysis of a Matrix object (dense or sparse)
- sparsity: computes sparsity of dense/sparse matrix
- sparse.colMeans: compute colMeans, either with or without zeroes included
- sparse.rowMeans: compute rowMeans, either with or without zeroes included
- sparse.quantile: compute quantile for a sparseVector object with or without zeroes included(use sparseMatrix slicing w/ "drop=FALSE" to keep the sparse form and prevent memory overflow)

### Covariance/Correlation of a Matrix object (dense or sparse)
- sparse.covariance
- sparse.correlation (from scratch)
- sparse.correlation2 (sparse.covariance + Matrix::cov2cor(cov), two approachces with close time performance)
