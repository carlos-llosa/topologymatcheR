

#' @export
Styp  <- function(x,inv = F){
  ei <- eigen(x,symmetric = T)
  if(inv) ei$values <- 1/ei$values
  matgen <- function(r) tcrossprod(t(ei$values^(r/2) * t(ei$vectors)))
  inv <- matgen(-1)
  list(invsqr = matgen(-0.5)  ,
       inv =   inv,
       mhalf = matgen(0.5),
       orig  =matgen(1),
       det = prod(ei$values),
       ldet = sum(log(ei$values)),
       sinv = sum(inv),
       eigen = ei,
       relnorm = sqrt(sum(ei$values^2))/ncol(x)
       )
}




#' Maximum Likelihood Estimation of the Matrix-t Distribution
#'
#' This function fits the matrix-variate-t distribution with pre-specified degrees of freedom, and under identical
#' column matrix and AR(1) correlation in the row-scale matrix
#'
#' This function is not as versatile as the one in MixMatrix but it works well with AR(1).
#' Will opt for MixMatrix in the future when it gets fixed
#'
#' @param Yall p x q x n array that contain the n indepedent p x q matrices
#' @param maxit integer - maximum number of iterations
#' @param df numeric - degrees of freedom for the Matrix-variate t
#' @param err  numeric - tolerance for convergence
#' @param init list of initial values with init$mean a p x q mean matrix, init$U p x p row-scale matrix and
#' init$V qxq column-scale matrix.
#' @return a list with elements\cr\cr
#'  mean - estimated p x q mean matrix \cr\cr
#'  U - estimated p x p row-scale matrix \cr\cr
#'  V - estimated q x q column-scale matrix \cr\cr
#'  df - degrees of freedom that were used \cr\cr
#'  iter - iterations it took \cr\cr
#'  tol - the error that was used to assess convergence  \cr\cr
#'  logLik - loglikelihood of the final fitted model \cr\cr
#'  covergence - T or F indicating convergence \cr\cr
#'  allconv - vector containing all the sum of relative norms, which were used to assess convergence \cr\cr
#' @export
#' @examples
#' set.seed(1234)
#' p <- 2
#' q <- 6
#' U <- drop(rWishart(1,p,diag(p)))
#' V <- toeplitz(0.2^(0:(q-1)))
#' m <- rgamma(p,2) %o% rep(1,q)
#' df <- 5
#' Y <- MixMatrix::rmatrixt(1000,df,mean = m,U=U,V=V)
#' fit <- MLmatrixt_ar1(Y,df = df)
#
#' list("true" = U , "est" = fit$U)
#' sapply(list("true" = V , "est" = fit$V),function(x) x[2,1])
#' sapply(list("true" = m , "est" = fit$mean),function(x) x[,1])
#' plot(ts(fit$allconv))
#'
#' @author Carlos Llosa-Vite
#' @references \url{https://doi.org/10.1080/10618600.2019.1696208}
#' @import MixMatrix
MLmatrixt_ar1 <- function(Yall,maxit = 1e4,df = 5,err = 1e-5,init = NULL){
  #assigning sizes
  di <- dim(Yall)
  p <- di[1]; q <- di[2]; n <- di[3]
  covgen <- function(rho) toeplitz(rho^(0:(q-1)))
  rang <- c(-1+1e-3,1-1e-3)
  #this is used for the optimization of AR(1) or CS
  Vgen <- function(rho,S){
    Vd <- chol(covgen(rho))
    p*n*sum(2*log(diag(Vd))) + sum(chol2inv(Vd)*S)
  }
  if(is.null(init)){
    mv <- apply(Yall,1,mean)
    Vtyp <- Styp(diag(q))
    Utyp <- Styp(diag(p))
  } else {
    mv <- apply(init$mean,1,mean)
    Vtyp <- Styp(init$V)
    Utyp <- Styp(init$U)
  }
  #not iterative
  allconv <- 0
  convergence <- F
  #estimate
  for(k in 1:maxit){
    #center the data (needed iteratively whenever the mean has a special structure)
    cent <- Yall - mv
    #E-step:  This finds the E(S_i|X_i) for all i = 1,..,n
    E1 <- lapply(1:n,function(i)chol2inv(chol(cent[,,i] %*% Vtyp$inv %*% t(cent[,,i]) + Utyp$orig)))
    E1 <- (df+p+q-1)*aperm(simplify2array(E1),c(3,1,2))
    #M-step for U : This is the pxp matrix that is left unstructured
    SU <- apply(E1,2:3,sum)
    Utyp <- Styp(SU/(n*(p+df-1)),inv = T)
    #M-step for V : This is the qxq matrix with AR(1) structure
    S <-  Reduce( "+" , lapply(1:n,function(i)t(cent[,,i]) %*% E1[i,,] %*% cent[,,i] ))
    rho_optim <- optimize(Vgen,rang,S=S)
    Vtyp <- Styp(covgen(rho_optim$minimum ))
    #M-step for m : This is a px1 vector. Recall the pxq matrix has q identical columns m
    mv <- Reduce("+",lapply(1:n,function(i)E1[i,,] %*% Yall[,,i]))
    mv <- apply(Utyp$orig %*% mv %*% Vtyp$inv,1,sum)/(Vtyp$sinv*n*(p+df-1))
    #calculate convergence criteria, which is the sum of relative (to size) norms
    conv <- Utyp$relnorm + Vtyp$relnorm + sqrt(mean(mv^2))
    allconv <- c(allconv,conv)
    #asses convergence
    if(abs((allconv[k+1]-allconv[k])/allconv[k])<err){
      #print(paste("converged at",k,"iterations with a relative tolerance of ", err))
      convergence <- T
      break
    }
  }
  allconv <- allconv[-1] #the first term was 0, so remove it
  llik <- sum(dmatrixt(x=Yall,df = df,mean = replicate(q,mv),U = Utyp$orig,V = Vtyp$orig,log = T))
  return(list(mean = replicate(q,mv),
              U = Utyp$orig,
              V = Vtyp$orig,
              df = df,
              iter = k,
              tol = err,
              logLik = llik,
              convergence = convergence,
              allconv = allconv))
}


# p <- 2
# q <- 6
# U <- drop(rWishart(1,p,diag(p)))
# V <- toeplitz(0.2^(0:(q-1)))
# m <- rgamma(p,2) %o% rep(1,q)
# df <- 5
# Y <- MixMatrix::rmatrixt(100,df,mean = m,U=U,V=V)
# fit <- MLmatrixt_ar1(Y,df = df)
#
# list("true" = U , "est" = fit$U)
# sapply(list("true" = V , "est" = fit$V),function(x) x[2,1])
# sapply(list("true" = m , "est" = fit$mean),function(x) x[,1])
# plot(ts(fit$allconv))

