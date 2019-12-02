# import dependencies
library(Matrix)
library(microbenchmark) # for speed comparison when writing different approaches

## A sparse matrix example to start with
# n <- 10000
# p <- 5000
# x <- matrix(rnorm(n * p), n, p)
# x <- matrix(rbinom(n*p, 80, 0.4),n,p)
# iz <- sample(1:(n * p),size = n * p * 0.85,replace = FALSE)
# x[iz] <- 0
# sx <- Matrix(x, sparse = TRUE)

############# Some useful functions in "Matrix" that can be used ###############
## colSums, rowSums, colMeans & rowMeans already implemented (with 0 included for calculation)
## mean, sum already implemented

# computes sparsity of dense/sparse matrix
sparsity <- function(x) {
  n = dim(x)[1]
  m = dim(x)[2]
  return (Matrix::nnzero(x)/n/m)
}

# compute colMeans, either with or without zeroes included.
sparse.colMeans <- function(x,zero.omit = FALSE,na.rm=FALSE,dims=1,sparseResult=FALSE) {
  if (!zero.omit) {
    return (Matrix::colMeans(x,na.rm,dims,sparseResult))
  } else {
    return (Matrix::colSums(x,na.rm,dims,sparseResult)/Matrix::colSums(abs(sign(x)),na.rm,dims,sparseResult))
  }
}

# compute rowMeans, either with or without zeroes included.
sparse.rowMeans <- function(x,zero.omit = FALSE,na.rm=FALSE,dims=1,sparseResult=FALSE) {
  if (!zero.omit) {
    return (Matrix::rowMeans(x,na.rm,dims,sparseResult))
  } else {
    return (Matrix::rowSums(x,na.rm,dims,sparseResult)/Matrix::rowSums(abs(sign(x)),na.rm,dims,sparseResult))
  }
}
