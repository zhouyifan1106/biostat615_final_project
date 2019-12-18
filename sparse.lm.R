set.seed(1000L) 
n = 10000L
p = 10L 
q = 6L
R2 = 0.9
X = rsparsematrix(nrow = n,ncol = p,nnz = n*0.1*p)
beta = sample(c(runif(q/2,min=2,max=5),runif(q/2,min=-5,max=-2),rep(0,p-q)))
mu = X%*%beta
sigma = sqrt(var(mu[,1])/R2*(1-R2)) 
Y = mu + rnorm(n,sd=sigma)

########################  Sparse Linear Regression  ########################
########## methods available: qr, svd, chol
sparse.lm <- function(X,Y,interaction=FALSE,method = "qr",interaction_threshold = 0.8){
  Y = as.matrix(Y)
  if (interaction && ncol(X) > 1) {
    ### add 2-way interaction terms
    # "apply" for sparse matrix
    sparse.apply <- function(X, FUN) {
      res = as.numeric(ncol(X))
      X2 = as(X, "dgTMatrix")
      tmp = tapply(X2@x, X2@j, FUN)
      res[as.integer(names(tmp))+1] = tmp
      res
    }
    # find columns of X with sparsity>interaction_threshold -- only these are used for interaction
    X_less_sparse = X[,which(sparse.apply(X,function(x) {1-Matrix::nnzero(x)/nrow(X)})<=interaction_threshold),drop=FALSE]
    
    X_interaction = lapply(1:(ncol(X_less_sparse)-1), function(i) {
      tmp = drop0(X_less_sparse[,(i+1):(ncol(X_less_sparse))]*X_less_sparse[,i],tol=1e-15)
      tmp = tmp[,which(sparse.apply(tmp,function(x) {1-Matrix::nnzero(x)/nrow(tmp)})<=interaction_threshold),drop=FALSE]
      return (tmp)})
    # if p is too large, do.call(cbind,X_interaction) might cause recursion max depth error; we use iterative approach to prevent this from happening
    X_interaction = Reduce(cbind, X_interaction[-1], X_interaction[[1]])
    X = cbind(X,X_interaction)
  }
  # add intercept term
  X = cbind(rep(1,nrow(X)),X)
  
  # start beta computation
  p = ncol(X) 
  n = nrow(X)
  #### We want to use Cholesky if we can, but crossprod(X) isn't always PD.
  #### One way to check PD is to do eigen(), and see if any eigenvalue is nonpositive. 
  # However, this process might be computational expensive. (O(p^3))
  # Also, eigenvalues might be numerically unstable, so we need to incorporate tolerance.
  check_PD <- function(x,tolerance = 1e-10) {
    e_values <- eigen(x,symmetric=TRUE,only.values=TRUE)$values
    # numerical stability
    e_values[which(abs(e_values) < tolerance)] <- 0
    return (!min(e_values <=0))
  }
  #### Another way to check is to directly compute chol(), use "try except" and see if there's a non-PD error message
  # the matrix is PD, the cost of checking is 0; if it's not, the cost is O(p^3) when doing cholesky.
  # However, while cholesky is accurate for exact precision arithmetic, it does not always work with floating point arithmetic
  # Here is an example with mathematical proof, where chol() fails but the matrix is indeed PD: https://math.stackexchange.com/questions/2147339/method-to-check-for-positive-definite-matrices/2149453#2149453
  
  ## We decided to use the second approach for both numerical stability and computation efficiency.
  tXX = crossprod(X)
  b = 0
  tryCatch({
      # cholesky
      R = chol(tXX)
      tXY = crossprod(X,Y)
      z = forwardsolve(R,tXY,upper.tri=TRUE,transpose=TRUE)
      b = backsolve(R,z)},
    error = function(e) {
      if (method == "svd") {
        SVD = sparsesvd(tXX)
        inv = SVD$v%*%diag(1/SVD$d)%*%t(SVD$u)
        b <<- inv %*% (t(X) %*% Y)
      } else {
        # here we still use the qr function from "base" (written in Fortran)
        # in "Matrix", there is qr for sparse matrix (in C), but it's significantly slower than the "base qr"
        # if memory saving is the priority, use qr(); otherwise, use qr.default().
        b <<- qr.solve(X, as.matrix(Y), tol = 1e-10)
      }}
    )
  #se =  (sqrt(diag(sum((Y-X%*%as.matrix(b))^2)/(nrow(X)-ncol(X)-1)*qr.solve(tXX))))
  return (b)
}

     
