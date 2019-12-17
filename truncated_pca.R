################### Trunccated PCA ###################
## This is a modification of the prcomp.default() function inside stats package.

## for a n*p matrix, the number of PCs in the output will be min(n,p). When p is large, the computation is very expensive.
## With a too-sparse matrix, not much information is inside each variable. We want to find a way to cap the number of eigenvalues/eigenvectors in the output and ignore unimportant PCs. 
## However,currently, prcomp() does not support the desired cap method.


## There are two ways of computing PCA once covariance matrix is ready -- getting eigenvectors and eigenvalues (package princomp) and doing SVD (package prcomp).
## We use the SVD approach, since svd() allows you to specify the maxinum number of eigenvectors/eigenvalues.

truncated.pca <- function(x, number_of_pc = 5,retx = TRUE, center = TRUE, scale. = FALSE, tol = NULL, ...){
    ## Sparse matrix after centering will be dense, so we do not worry too much and directly convert the matrix to dense form
    x <- as.matrix(x)
    # scale and center the data as needed
    x <- scale(x, center = center, scale = scale.)
    cen <- attr(x, "scaled:center")
    sc <- attr(x, "scaled:scale")
    if(any(sc == 0))
      stop("cannot rescale a constant/zero column to unit variance")
    # only compute the top "number_of_pc" PCs for the output
    s <- svd(x,number_of_pc,number_of_pc)
    s$d <- s$d/sqrt(max(1, nrow(x) - 1))
    if (!is.null(tol)) {
      rank <- sum(s$d > (s$d[1L]*tol))
      if (rank < ncol(x)) {
        s$v <- s$v[, 1L:rank, drop = FALSE]
        s$d <- s$d[1L:rank]
      }
    }
    dimnames(s$v) <- list(colnames(x), paste0("PC", seq_len(ncol(s$v))))
    r <- list(sdev = s$d, rotation = s$v,center = if(is.null(cen)) FALSE else cen,scale = if(is.null(sc)) FALSE else sc)
    if (retx) {
      r$x <- x %*% s$v
    }
    r
  }

plot.prcomp <- function(x, main = deparse(substitute(x)), ...)
  screeplot.default(x, main = main, ...)

print.prcomp <- function(x, print.x = FALSE, ...) {
  cat("Standard deviations:\n")
  print(x$sdev, ...)
  cat("\nRotation:\n")
  print(x$rotation, ...)
  if (print.x && length(x$x)) {
    cat("\nRotated variables:\n")
    print(x$x, ...)
  }
  invisible(x)
}

summary.prcomp <- function(object, ...){
  vars <- object$sdev^2
  vars <- vars/sum(vars)
  importance <- rbind("Standard deviation" = object$sdev,
                      "Proportion of Variance" = round(vars, 5),
                      "Cumulative Proportion" = round(cumsum(vars), 5))
  colnames(importance) <- colnames(object$rotation)
  object$importance <- importance
  class(object) <- "summary.prcomp"
  object
}

print.summary.prcomp <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
    cat("Importance of components:\n")
    print(x$importance, digits = digits, ...)
    invisible(x)
  }

predict.prcomp <- function(object, newdata, ...) {
  if (missing(newdata)) {
    if(!is.null(object$x)) return(object$x)
    else stop("no scores are available: refit with 'retx=TRUE'")
  }
  if(length(dim(newdata)) != 2L)
    stop("'newdata' must be a matrix or data frame")
  nm <- rownames(object$rotation)
  if(!is.null(nm)) {
    if(!all(nm %in% colnames(newdata)))
      stop("'newdata' does not have named columns matching one or more of the original columns")
    newdata <- newdata[, nm, drop = FALSE]
  } else {
    if(NCOL(newdata) != NROW(object$rotation) )
      stop("'newdata' does not have the correct number of columns")
  }
  scale(newdata, object$center, object$scale) %*% object$rotation
}

.check_vars_numeric <- function(mf) {
  mt <- attr(mf, "terms")
  mterms <- attr(mt, "factors")
  mterms <- rownames(mterms)[apply(mterms, 1L, function(x) any(x > 0L))]
  any(sapply(mterms, function(x) is.factor(mf[,x]) || !is.numeric(mf[,x])))
}
