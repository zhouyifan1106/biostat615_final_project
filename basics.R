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
  p = dim(x)[2]
  return (1-Matrix::nnzero(x)/n/p)
}

# compute colMeans, either with or without zeroes included.
sparse.colMeans <- function(x,zero.omit = FALSE,na.rm=FALSE,dims=1) {
  if (!zero.omit) {
    return (Matrix::colMeans(x,na.rm,dims))
  } else {
    return (Matrix::colSums(x,na.rm,dims)/Matrix::colSums(abs(sign(x)),na.rm,dims))
  }
}

# compute rowMeans, either with or without zeroes included.
sparse.rowMeans <- function(x,zero.omit = FALSE,na.rm=FALSE,dims=1) {
  if (!zero.omit) {
    return (Matrix::rowMeans(x,na.rm,dims))
  } else {
    return (Matrix::rowSums(x,na.rm,dims,sparseResult)/Matrix::rowSums(abs(sign(x)),na.rm,dims))
  }
}

# Do "sparse_matrix[,i,drop=FALSE]" to keep sparse format
sparse.quantile <- function(sparse_vector,probs = seq(0, 1, 0.25),zero.omit = FALSE) {
  if (!zero.omit) {
    # consider zeroes
    return (quantile(as.vector(sparse_vector),probs))
  } else {
    # remove zeroes before quantile
    return (quantile(sparse_vector@x,probs))
  }
}

# print out sparsity,mean and quantile for each column of the sparseMatrix (does not work for dense matrix)
sparse.summary <- function(x,zero.omit=FALSE) {
  # variable summary
  sparse.variable_summary <- function(variable,zero.omit=FALSE) {
    quantiles = sparse.quantile(variable,zero.omit = zero.omit)
    sparsity = sparsity(variable)
    return (c(quantiles,sparsity))
  }
  output = sapply(1:ncol(x),function (i) {sparse.variable_summary(x[,i,drop=FALSE],zero.omit = zero.omit)})
  mean = sparse.colMeans(x,zero.omit = zero.omit)
  output = rbind(mean,output)
  rownames(output) = c("Mean","Min.","1st Qu.","Median","3rd Qu.","Max.","sparsity")
  output = output[c(7,2,3,4,1,5,6),]
  return (output)
}

# Covariance for Matrix object (accepts both dense and sparse input, although dense matrix should probably use cov() from base)
sparse.covariance <- function(x){
  n = nrow(x)
  cSums = colSums(x)
  cMeans = colMeans(x)
  covariance = tcrossprod(cMeans,(-2*cSums+n*cMeans))+as.matrix(crossprod(x))
  return (covariance/(n-1))
}

## Correlation for Matrix object -- accepts both dense and sparse input, although dense matrix should probably use cov() from base
## two approaches have almost identital time performance

# Correlation approach 1
sparse.correlation <- function(x) {
  n = nrow(x)
  cSums = colSums(x)
  cMeans = colMeans(x)
  covariance = tcrossprod(cMeans,(-2*cSums+n*cMeans))+as.matrix(crossprod(x))
  correlation = covariance/crossprod(t(sqrt(diag(covariance))))
  return (correlation)
}

# Correlation approach 2
sparse.correlation2 <- function(x) {
  covariance = sparse.covariance(x)
  return (Matrix::cov2cor(covariance))
}
