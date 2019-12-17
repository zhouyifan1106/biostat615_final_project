## Biostatistics 615 Final Project - Sparse Matrix Toolkit for package "Matrix"

### Exploratory Data Analysis functions (basics.R)
#### All functions accept sparse and dense Matrix objects
- **sparsity**
- **sparse.colMeans** (with or without zeroes)
- **sparse.rowMeans** (with or without zeroes)
- **sparse.quantile** (with or without zeroes)
- **sparse.summary** outputs summary statistics, including sparsity,min,Q1,median,mean,Q3,and max(with or without zeroes)
- **sparse.covariance**
- **sparse.correlation** (2 implementations)

### Sparse Linear Regression with optional 2-way interaction and user-defined sparsity threshold for variable selection
- sparse.lm

### Truncated PCA for extremely sparse matrices with large number of features
- truncated_pca
