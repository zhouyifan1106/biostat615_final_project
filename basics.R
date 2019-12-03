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

# sparse linear regression with optional 2-way interaction 
# methods available: qr, svd, chol
sparse.lm <- function(X,Y,interaction=FALSE,method = "qr"){
  Y = as.matrix(Y)
  if (interaction && ncol(X) > 1) {
    # add 2-way interaction terms
    sparse.cbind <- function (...) {
      input = lapply(list(...),as,"sparseVector" )
      l = unique(sapply(input,length))
      return(sparseMatrix(x=unlist(lapply(input,slot,"x")), i=unlist(lapply(input,slot,"i")), 
                          p=c(0,cumsum(sapply(input,function(x){length(x@x)}))),dims=c(l,length(input))))
    }
    X_interaction = unlist(lapply(1:(ncol(X)-1), function(i) 
      lapply((i+1):ncol(X), function(j) as(X[,i]*X[,j],"sparseVector"))))
    X_interaction = do.call(sparse.cbind,X_interaction)
    X = cbind(X,X_interaction)
  }
  # add intercept term
  X = cbind(rep(1,nrow(X)),X)
  # start beta computation
  p = ncol(X) 
  n = nrow(X)
  if (method == "qr") {
    QR = qr.default(X)
    b = backsolve(QR$qr, qr.qty(QR, Y))
  } else if (method == "svd") {
    tXX = crossprod(X) 
    SVD = svd(tXX)
    inv = SVD$v%*%diag(1/SVD$d)%*%t(SVD$u)
    b = inv %*% (t(X) %*% Y)
  } else { #chol
    tXX = crossprod(X) 
    tXY = crossprod(X,Y)
    R = chol(tXX)
    z = forwardsolve(R,tXY,upper.tri=TRUE,transpose=TRUE)
    b = backsolve(R,z)
  }
  return (b)
}
