#' Simulate Data from the Full RePlei Generative Model
#'
#' @description
#' This function simulates summary statistics from the full RePlei generative model,
#' which includes pleiotropic effects, study-specific heterogeneity, and
#' non-replicable outlier effects among null SNPs.
#'
#' @details
#' This function implements the complete data generation process for RePlei. It simulates
#' a proportion of true pleiotropic SNPs (`p`) whose effects are drawn from a
#' distribution governed by genetic correlation (`C_g`). A proportion of the remaining
#' null SNPs (`lambda`) are simulated as outliers with inflated variance (`alpha`).
#'
#' By setting `lambda = 0`, this function can also simulate data under a simpler
#' fixed-effect pleiotropy model without outliers.
#'
#' @note
#' For this function to work seamlessly within your package, the `res2ref_cpd` data object
#' should be included in the package's `data/` directory. You can do this by running
#' `usethis::use_data(res2ref_cpd, overwrite = TRUE)` from your package's root directory.
#'
#' @param M An integer. The number of SNPs to simulate. Default: `5000`.
#' @param p A numeric value. The proportion of SNPs with true non-zero pleiotropic effects. Default: `0.01`.
#' @param tau2 A numeric value. The variance scaling factor for the true effect sizes (`mu_j`). Default: `5e-04`.
#' @param C_g A `Q x Q` numeric matrix. The genetic correlation matrix between the `Q` traits.
#' @param C_e A `Q x Q` numeric matrix. The environmental (or residual) correlation matrix between the `Q` traits.
#' @param h2 A numeric vector of length `Q`. The heritability for each of the `Q` traits.
#' @param Q An integer. The number of traits. Default: `3`.
#' @param n A numeric vector of length `K`. The sample sizes for each of the `K` studies. Default: `rep(50000, 10)`.
#' @param lambda A numeric value. The proportion of null SNPs that are outliers in a given study. Default: `0.025`.
#' @param alpha A numeric value. The inflation factor for the variance of outlier SNPs. Default: `5`.
#' @param resref A numeric vector. A reference pool of residual standard deviations to sample from. If `NULL` (default), the function attempts to use the `res2ref_cpd` data object.
#'
#' @return A list containing the full simulated data, including `betajk`, `sjk`, `Rj`, and `Ojk`.
#'
#' @export
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats rbinom rnorm
#'
#' @examples
#' # --- Setup parameters for a 3-trait simulation ---
#' Q_val <- 3
#' Cg_mat <- matrix(c(1, 0.4, 0.2, 0.4, 1, 0.3, 0.2, 0.3, 1), nrow = Q_val)
#' Ce_mat <- matrix(c(1, 0.1, 0.05, 0.1, 1, 0.15, 0.05, 0.15, 1), nrow = Q_val)
#' h2_vec <- c(0.3, 0.4, 0.25)
#' sample_sizes <- rep(60000, 5) # 5 studies
#'
#' # --- Simulate data with outlier effects ---
#' set.seed(123)
#' sim_data_replei <- simulate_replei_data(
#'   M = 1000, p = 0.02, C_g = Cg_mat, C_e = Ce_mat, h2 = h2_vec,
#'   Q = Q_val, n = sample_sizes, lambda = 0.03, alpha = 5
#' )
#'
simulate_replei_data <- function(M = 5000, p = 0.01, tau2 = 5e-04, C_g, C_e, h2,
                                 Q = 3, n = rep(50000, 10), lambda = 0.025,
                                 alpha = 5, resref = NULL) {
  # --- 1. Input Validation ---
  if (missing(C_g) || missing(C_e)) {
    stop("Error: Both C_g and C_e must be provided.")
  }
  if (length(h2) != Q) {
    stop("Error: Length of h2 must match Q.")
  }
  
  if (is.null(resref)) {
    if (exists("res2ref_cpd")) {
      resref <- sqrt(res2ref_cpd$res2)
    } else {
      warning("Using a default 'resref'. Provide 'resref' or add 'res2ref_cpd' to package data.")
      resref <- stats::rnorm(1000, mean = 1, sd = 0.1)
    }
  }
  
  K <- length(n)
  
  # --- 2. Simulate True Genetic Effects (mu_j) ---
  h_sqrt <- sqrt(h2)
  Omega <- diag(h_sqrt) %*% C_g %*% diag(h_sqrt)
  Rj <- stats::rbinom(M, size = 1, prob = p)
  Rj_1_indices <- which(Rj == 1)
  Rj_0_indices <- which(Rj == 0)
  num_causal <- length(Rj_1_indices)
  
  mu_j <- matrix(0, nrow = M, ncol = Q)
  if (num_causal > 0) {
    mu_j[Rj_1_indices, ] <- mvtnorm::rmvnorm(n = num_causal, mean = rep(0, Q), sigma = tau2 * Omega)
  }
  
  # --- 3. Simulate Variances and Outliers ---
  res_sd_array <- array(sample(resref, M * K * Q, replace = TRUE), dim = c(M, K, Q))
  sjk2 <- (res_sd_array^2) / array(rep(n, each = M), dim = c(M, K, Q))
  sjk <- sqrt(sjk2)
  
  Ojk <- matrix(0, nrow = M, ncol = K)
  num_null <- length(Rj_0_indices)
  if (num_null > 0 && lambda > 0) {
    Ojk[Rj_0_indices, ] <- matrix(stats::rbinom(num_null * K, size = 1, prob = lambda),
                                  nrow = num_null, ncol = K)
  }
  
  # --- 4. Simulate Observed Effects (betajk) ---
  betajk <- array(0, dim = c(M, K, Q))
  for (k in 1:K) {
    if (num_causal > 0) {
      for (j in Rj_1_indices) {
        Sigma_e_jk <- diag(sjk[j, k, ]) %*% C_e %*% diag(sjk[j, k, ])
        betajk[j, k, ] <- mvtnorm::rmvnorm(1, mean = mu_j[j, ], sigma = Sigma_e_jk)
      }
    }
    if (num_null > 0) {
      for (j in Rj_0_indices) {
        inflation_sqrt <- if (Ojk[j, k] == 1) sqrt(alpha) else 1.0
        sjk_jk_inflated <- sjk[j, k, ] * inflation_sqrt
        Sigma_e_jk <- diag(sjk_jk_inflated) %*% C_e %*% diag(sjk_jk_inflated)
        betajk[j, k, ] <- mvtnorm::rmvnorm(1, mean = rep(0, Q), sigma = Sigma_e_jk)
      }
    }
  }
  
  # --- 5. Return Results ---
  return(list(betajk = betajk, sjk2 = sjk2, sjk = sjk, Rj = Rj, Ojk = Ojk,
              mu_j = mu_j, C_g = C_g, C_e = C_e, h2 = h2))
}
