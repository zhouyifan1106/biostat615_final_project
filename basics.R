# import dependencies
library(Matrix)

######## A sparse matrix example to start with ########
# x: dense form
# sx: sparse form
n <- 1000000
p <- 5
x <- matrix(rnorm(n * p), n, p)
x <- matrix(rbinom(n*p, 80, 0.4),n,p)
iz <- sample(1:(n * p),size = n * p * 0.85,replace = FALSE)
x[iz] <- 0
sx <- Matrix(x, sparse = TRUE)


########################################## FUNCTIONS ################################################

###################  sparsity  ###################
## x: dense or sparse Matrix object
sparsity <- function(x) {
  n = dim(x)[1]
  p = dim(x)[2]
  return (1-Matrix::nnzero(x)/n/p)
}

###################  colMeans  ###################
## x: dense or sparse Matrix object
## zero.omit: if TRUE, only nonzero values are used for computation
## na.rm: if TRUE, omit NA from calculations
sparse.colMeans <- function(x,zero.omit = FALSE,na.rm=FALSE) {
  if (!zero.omit) {
    return (Matrix::colMeans(x,na.rm))
  } else if (is(x,'sparseMatrix')) {
    return (Matrix::colSums(x,na.rm,sparseResult)/Matrix::colSums(abs(sign(x)),na.rm,sparseResult))
  } else {
    return (Matrix::colSums(x,na.rm)/Matrix::colSums(abs(sign(x)),na.rm))
  }
}

################### rowMeans  ###################
## x: dense or sparse Matrix object
## zero.omit: if TRUE, only nonzero values are used for computation
## na.rm: if TRUE, omit NA from calculations
sparse.rowMeans <- function(x,zero.omit = FALSE,na.rm=FALSE) {
  if (!zero.omit) {
    return (Matrix::rowMeans(x,na.rm))
  } else if (is(x,'sparseMatrix')) {
    # keep result as sparse for memory saving
    return (Matrix::rowSums(x,na.rm,sparseResult)/Matrix::rowSums(abs(sign(x)),na.rm,sparseResult))
  } else {
    # dense matrix
    return (Matrix::rowSums(x,na.rm)/Matrix::rowSums(abs(sign(x)),na.rm))
  }
}

################### quantile ###################
## x: dense or sparse Matrix object
## probs: numeric vector of probabilities with values in [0,1]
## zero.omit: if TRUE, only nonzero values are used for computation
sparse.quantile <- function(x,probs = seq(0,1,0.25),zero.omit = FALSE) {
  if (!zero.omit) {
    return (quantile(as.vector(x),probs))
  } else if (is(x, 'sparseMatrix')) {
    # remove zero for sparse matrix
    return (quantile(x@x,probs))
  } else {
    # remove zero for dense matrix
    return (quantile(x[-which(x==0)],probs))
  }
}

################### summary (sparsity,min,Q1,median,mean,Q3,max) ###################
## x: dense or sparse Matrix object
## zero.omit: if TRUE, include zero values for summary statistics
sparse.summary <- function(x,zero.omit=FALSE) {
  # computate QUANTILE & SPARSITY for each variable
  sparse.variable_summary <- function(variable,zero.omit=FALSE) {
    quantiles = sparse.quantile(variable,zero.omit = zero.omit)
    sparsity = sparsity(variable)
    return (c(quantiles,sparsity))
  }
  # use "drop=FALSE" to keep sparse format during slicing
  output = sapply(1:ncol(x),function (i) {sparse.variable_summary(x[,i,drop=FALSE],zero.omit = zero.omit)})
  
  # compute MEAN for each variable
  mean = sparse.colMeans(x,zero.omit = zero.omit)
  # combine outputs and sort columns
  output = rbind(mean,output)
  rownames(output) = c("Mean","Min.","1st Qu.","Median","3rd Qu.","Max.","sparsity")
  output = output[c(7,2,3,4,1,5,6),]
  return (output)
}

################### covariance  ###################
## x: dense or sparse Matrix object
sparse.covariance <- function(x){
  if (is(x,'sparseMatrix')) {
    n = nrow(x)
    cSums = colSums(x)
    cMeans = colMeans(x)
    covariance = tcrossprod(cMeans,(-2*cSums+n*cMeans))+as.matrix(crossprod(x))
    return (covariance/(n-1))
  } else {
    return (cov(x))
  }
}

################### correlation  ###################
## Two implementations, almost identical runtime
## x: dense or sparse Matrix object
# Approach 1 -- implement from scratch
sparse.correlation <- function(x) {
  if (is(x,'sparseMatrix')) {
    n = nrow(x) 
    cSums = colSums(x)
    cMeans = colMeans(x)
    covariance = tcrossprod(cMeans,(-2*cSums+n*cMeans))+as.matrix(crossprod(x))
    return (covariance/crossprod(t(sqrt(diag(covariance)))))
  } else {
    return (cor(x))
  }
}
# Approach 2 -- use previously implemented "sparse.covariance" and Matrix::cov2cor
sparse.correlation2 <- function(x) {
  if (is(x,'sparseMatrix')) {
    covariance = sparse.covariance(x)
    return (Matrix::cov2cor(covariance))
  } else {
    return (cor(x))
  }
}
