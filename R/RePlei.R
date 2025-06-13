#' @useDynLib RePlei, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @import RcppArmadillo
NULL
#' Fast EM for RePlei
#'
#' Run the Expectation-Maximization algorithm to estimate replicability
#' @param betajk Array of standardized effect sizes (dimensions M x K x Q).
#' @param sjk Array of standard errors corresponding to betajk (M x K x Q).
#' @param C_g Genetic correlation matrix between traits.
#' @param C_e Environmental correlation matrix between traits.
#' @param h2 Vector of heritabilities for each trait.
#' @param p Initial guess for the proportion of replicable SNPs.
#' @param lambda Initial guess for the proportion of outlier null SNPs.
#' @param alpha Initial guess for the variance inflation of outliers.
#' @param tau2 Initial guess for the variance of true effects.
#' @param min_iter Minimum number of EM iterations before checking convergence.
#' @param max_iter Maximum number of EM iterations.
#' @param rel.eps Relative tolerance for EM convergence.
#' @param verbose Logical; print progress during EM iterations.
#' @param nthreads Number of CPU threads for OpenMP parallelization.
#' @return A list containing:
#' \describe{
#'   \item{p}{Estimated proportion of replicable SNPs}
#'   \item{lambda}{Estimated proportion of outlier null SNPs}
#'   \item{alpha}{Estimated inflation factor for outliers}
#'   \item{tau2}{Estimated variance of true effects}
#'   \item{R_hat}{Posterior probability of replicability (PPR) for each SNP}
#'   \item{O_hat}{Posterior outlier probabilities for each study/SNP}
#'   \item{mu_hat}{Posterior mean effect size vectors for each SNP}
#'   \item{post_means}{PPR-weighted posterior effect estimates}
#'   \item{FDR}{Estimated false discovery rate for each SNP}
#'   \item{logLik}{Log-likelihood values across iterations}
#' }
#' @export
RePlei <- function(betajk, sjk, C_g, C_e, h2,
                   p = 0.1, lambda = 0.10, alpha = 2, tau2 = 5e-4,
                   min_iter = 100, max_iter = 1000,
                   rel.eps = 1e-6, verbose = TRUE,
                   nthreads = parallel::detectCores()) {
  
  stopifnot(length(dim(betajk)) == 3, all(dim(betajk) == dim(sjk)))
  Sys.setenv(OMP_NUM_THREADS = nthreads)
  
  Omega  <- diag(sqrt(h2)) %*% C_g %*% diag(sqrt(h2))
  OmInv  <- solve(Omega)
  
  M <- dim(betajk)[1]; K <- dim(betajk)[2]; Q <- dim(betajk)[3]
  ll_hist <- numeric(max_iter)
  
  for (it in seq_len(max_iter)) {
    #E step
    E_res <- E_step(
      aperm(betajk, c(3, 2, 1)),
      aperm(sjk,    c(3, 2, 1)),
      C_e, Omega, OmInv, lambda, alpha, tau2, p)
    
    R_hat    <- E_res$Rhat
    ll_total <- sum(E_res$ll)
    alph_num <- sum(E_res$alph_num)
    alph_den <- sum(E_res$alph_den)
    tau_num  <- sum(E_res$tau_num)
    
    #M step
    p_new      <- mean(R_hat)
    lambda_new <- alph_den / (K * sum(1 - R_hat))
    alpha_new  <- alph_num / (Q * alph_den)
    tau2_new   <- max(tau_num / (Q * sum(R_hat)), 1e-6)
    
    # eps
    eps <- 1e-9
    p_new      <- max(min(p_new,      1 - eps), eps)
    lambda_new <- max(min(lambda_new, 1 - eps), eps)
    alpha_new  <- max(alpha_new, eps)
    tau2_new   <- max(tau2_new,  eps)
    
    #Convergence
    ll_hist[it] <- ll_total
    if (verbose) {
      msg <- sprintf(
        "iter %03d  ll=%.3f  p=%.4g  lambda=%.4g  alpha=%.3f  tau2=%.2e",
        it, ll_total, p_new, lambda_new, alpha_new, tau2_new)
      message(msg)
    }
    
    # Check log likelihood criteria after the first iteration
    if (it > 1) {
      ll_rel  <- abs((ll_hist[it] - ll_hist[it - 1]) / ll_hist[it])
      ll_diff <- abs(ll_hist[it] - ll_hist[it - 1])
      
      converged <- (
        (ll_rel  < rel.eps && ll_diff < 0.01) ||
          (ll_diff < 0.001    && it > min_iter)   ||
          (ll_rel  < rel.eps  && it > min_iter))
      
      if (converged) {
        ll_hist <- ll_hist[seq_len(it)]
        break
      }
    }
    
    # Update parameters for next iteration ---------------------------------
    p      <- p_new
    lambda <- lambda_new
    alpha  <- alpha_new
    tau2   <- tau2_new
  }
  
  list(p = p, lambda = lambda, alpha = alpha, tau2 = tau2,
       R_hat = as.vector(E_res$Rhat), O_hat = E_res$Ohat, mu_hat = E_res$mu_hat,
       post_means = E_res$mu_hat * as.vector(E_res$Rhat),
       FDR = Get_FDR(E_res$Rhat),
       logLik = ll_hist)
  
}

#' Calculate FDR
#' @param ppr Posterior probabilities
#' @return FDR estimates
#' @export
Get_FDR <- function(ppr){
  o <- order(1 - ppr); f <- cumsum(1 - ppr[o]) / seq_along(ppr)
  f[order(o)]
}
