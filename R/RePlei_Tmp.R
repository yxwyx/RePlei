#' @useDynLib RePlei, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @import RcppArmadillo
NULL

#' RePlei EM 
#'
#' @param betajk 3-D array of effect sizes (M × K × Q). NA where missing.
#' @param sjk    3-D array of standard errors (same dims). NA where missing.
#' @param C_g    Trait genetic correlation matrix (Q × Q).
#' @param C_e    Trait environmental correlation matrix (Q × Q).
#' @param h2     Length-Q vector of trait heritabilities.
#' @param p,lambda,alpha,tau2 Starting values (see paper).
#' @param min_iter,max_iter,rel.eps Convergence controls.
#' @param verbose Print progress?
#' @param nthreads OpenMP threads to use.
#'
#' @return List with estimates, posteriors, logLik trajectory (unchanged).
#' @export
RePlei <- function(betajk, sjk, C_g, C_e, h2,
                   p = 0.10, lambda = 0.10, alpha = 2,
                   tau2 = 5e-4,
                   min_iter = 100, max_iter = 1000,
                   rel.eps = 1e-6, verbose = TRUE,
                   nthreads = parallel::detectCores()) {
  
  ## --- basic sanity -------------------------------------------------------
  stopifnot(length(dim(betajk)) == 3,
            all(dim(betajk) == dim(sjk)),
            ncol(C_g) == length(h2))
  
  Sys.setenv(OMP_NUM_THREADS = nthreads)
  
  Omega  <- diag(sqrt(h2)) %*% C_g %*% diag(sqrt(h2))
  OmInv  <- solve(Omega)
  
  M <- dim(betajk)[1]; K <- dim(betajk)[2]; Q <- dim(betajk)[3]
  ll_hist <- numeric(max_iter)
  
  for (it in seq_len(max_iter)) {
    ## --- E-step -----------------------------------------------------------
    E <- E_step(
      aperm(betajk, c(3, 2, 1)),
      aperm(sjk,    c(3, 2, 1)),
      C_e, Omega, OmInv, lambda, alpha, tau2, p)
    
    R_hat    <- E$Rhat
    ll_total <- sum(E$ll)
    alph_num <- sum(E$alph_num)
    alph_den <- sum(E$alph_den)
    tau_num  <- sum(E$tau_num)
    
    ## --- M-step -----------------------------------------------------------
    p_new      <- mean(R_hat)
    lambda_new <- alph_den / (K * sum(1 - R_hat))
    alpha_new  <- alph_num / (Q * alph_den)
    tau2_new   <- max(tau_num / (Q * sum(R_hat)), 1e-6)
    
    eps <- 1e-9
    p_new      <- max(min(p_new,      1 - eps), eps)
    lambda_new <- max(min(lambda_new, 1 - eps), eps)
    alpha_new  <- max(alpha_new, eps)
    tau2_new   <- max(tau2_new,  eps)
    
    ll_hist[it] <- ll_total
    if (verbose) {
      message(sprintf(
        "iter %03d  ll=%.3f  p=%.4g  lambda=%.4g  alpha=%.3f  tau2=%.2e",
        it, ll_total, p_new, lambda_new, alpha_new, tau2_new))
    }
    
    if (it > 1) {
      ll_rel  <- abs((ll_hist[it] - ll_hist[it - 1]) / ll_hist[it])
      ll_diff <- abs(ll_hist[it] - ll_hist[it - 1])
      converged <- (ll_rel < rel.eps && ll_diff < 0.01) ||
        (ll_diff < 0.001 && it > min_iter)   ||
        (ll_rel < rel.eps && it > min_iter)
      if (converged) {
        ll_hist <- ll_hist[seq_len(it)]
        break
      }
    }
    
    ## update parameters
    p      <- p_new
    lambda <- lambda_new
    alpha  <- alpha_new
    tau2   <- tau2_new
  }
  
  list(p = p, lambda = lambda, alpha = alpha, tau2 = tau2,
       R_hat = as.vector(E$Rhat),
       O_hat = E$Ohat,
       mu_hat = E$mu_hat,
       post_means = E$mu_hat * as.vector(E$Rhat),
       FDR = Get_FDR(E$Rhat),
       logLik = ll_hist)
}

#' Benjamini–Hochberg style local FDR estimate
#' @param ppr Posterior probability of replicability (length M)
#' @export
Get_FDR <- function(ppr) {
  o <- order(1 - ppr)
  f <- cumsum(1 - ppr[o]) / seq_along(ppr)
  f[order(o)]
}

setwd("~/thesis/RePlei/R")
Rcpp::sourceCpp("replei.cpp")
load("~/thesis/RePlei_Analysis/simulations/RePlei/new_sims_M5000/Cg_S1/Ce_S1/Sim_M5000_K2_Q3_p0.1_lambda0.1_alpha2_tau25e-04_CgS1_CeS1_rep1.RData")
betajk <- sim_data$betajk
sjk    <- sim_data$sjk
C_g    <- sim_data$C_g
C_e    <- sim_data$C_e
h2     <- sim_data$h2
replei_fit <- RePlei(
  betajk  = betajk,
  sjk     = sjk,
  C_g     = C_g,
  C_e     = C_e,
  h2      = h2
)

set.seed(1)
M <- 200; K <- 3; Q <- 4
B <- array(rnorm(M*K*Q), c(M,K,Q))
S <- array(runif(M*K*Q, 0.05, 0.12), dim(B))

drop <- runif(M*K*Q) < .25           # 25 % missing completely at random
B[drop] <- NA
S[drop] <- NA

fit <- RePlei(B, S, C_g = diag(Q), C_e = diag(Q), h2 = rep(.5, Q),
              nthreads = 1)          # turn OMP off for first test
print(tail(fit$logLik))
