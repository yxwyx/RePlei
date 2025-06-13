// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#include <RcppArmadillo.h>
#ifdef _OPENMP
# include <omp.h>
#endif
#include <vector>

using namespace Rcpp;
using arma::mat;   using arma::vec;   using arma::rowvec;
using arma::uvec;  using arma::uword; using arma::cube;

/* ------------------------------------------------------------------ */
/*                              helpers                                */
/* ------------------------------------------------------------------ */

inline double lse(double a, double b) noexcept {
  return (a > b) ? a + std::log1p(std::exp(b - a))
    : b + std::log1p(std::exp(a - b));
}

inline bool safe_chol(mat &L, const mat &A) {
  if (arma::chol(L, A, "lower")) return true;
  mat Ap = A + 1e-10 * arma::eye(A.n_rows, A.n_cols);
  return arma::chol(L, Ap, "lower");
}

/* extract full row k of a matrix into a *plain* rowvec */
inline rowvec copy_row(const mat &M, uword k)
{
  rowvec out(M.n_cols);
  for (uword j = 0; j < M.n_cols; ++j) out[j] = M(k, j);
  return out;
}

/* pull selected elements into a new rowvec */
inline rowvec subset(const rowvec &r, const uvec &idx)
{
  rowvec out(idx.n_elem);
  for (uword i = 0; i < idx.n_elem; ++i) out[i] = r[idx[i]];
  return out;
}

/* indices of finite elements in a rowvec */
inline uvec finite_idx(const rowvec &r)
{
  std::vector<uword> tmp;
  for (uword i = 0; i < r.n_elem; ++i)
    if (arma::is_finite(r[i])) tmp.push_back(i);
    return uvec(tmp);
}

/* Σ_s = D Ce_sub D, return Σ and its chol L */
inline void make_sigma(const rowvec &s_row, const uvec &idx,
                       const mat &Ce, mat &Sig, mat &L)
{
  mat D(idx.n_elem, idx.n_elem, arma::fill::zeros);
  for (uword i = 0; i < idx.n_elem; ++i) D(i, i) = s_row[idx[i]];
  Sig = D * Ce.submat(idx, idx) * D;
  safe_chol(L, Sig);
}

/* log-dens of one MVN obs given chol(L), logdet */
inline double dmvnrm1(const rowvec &x,
                      const mat    &L,
                      double        logdet)
{
  const double k2Pi = std::log(2.0 * M_PI);
  rowvec z = (x) * arma::inv( arma::trimatu(L.t()) );
  return -0.5 * (x.n_elem * k2Pi + logdet + arma::dot(z, z));
}

/* ------------------------------------------------------------------ */
/*                log-likelihood pieces (missing aware)                */
/* ------------------------------------------------------------------ */

double ll_R0_j(const mat &b, const mat &s, const mat &Ce,
               double lambda, double alpha)
{
  double loglik = 0.0;
  
  for (uword k = 0; k < b.n_rows; ++k) {
    rowvec srow = copy_row(s, k);
    uvec   idx  = finite_idx(srow);
    if (idx.is_empty()) continue;
    
    mat Sig, L;  make_sigma(srow, idx, Ce, Sig, L);
    double ld = 2.0 * arma::sum(log(L.diag()));
    
    rowvec x  = subset(copy_row(b, k), idx);
    double ln0 = dmvnrm1(x, L, ld);
    
    mat   La  = std::sqrt(alpha) * L;
    double ldA = ld + idx.n_elem * std::log(alpha);
    double lna = dmvnrm1(x, La, ldA);
    
    loglik += lse(std::log(1.0 - lambda) + ln0,
                  std::log(lambda)       + lna);
  }
  return loglik;
}

double ll_R1_j(const mat &b, const mat &s, const mat &Ce,
               const mat &Omega, double tau2)
{
  std::vector<uvec>   idcs;
  std::vector<mat>    Sigks;
  std::vector<uword>  starts, obs_row;
  uword big = 0;
  
  for (uword k = 0; k < b.n_rows; ++k) {
    rowvec srow = copy_row(s, k);
    uvec   idx  = finite_idx(srow);
    if (idx.is_empty()) continue;
    
    mat Sig, dummy;  make_sigma(srow, idx, Ce, Sig, dummy);
    
    idcs.push_back(idx);
    Sigks.push_back(std::move(Sig));
    starts.push_back(big);
    obs_row.push_back(k);
    big += idx.n_elem;
  }
  if (big == 0) return 0.0;
  
  mat SigTot(big, big, arma::fill::zeros);
  /* diag sampling blocks */
  for (std::size_t bidx = 0; bidx < Sigks.size(); ++bidx) {
    uword st = starts[bidx];
    SigTot.submat(st, st,
                  st + idcs[bidx].n_elem - 1,
                  st + idcs[bidx].n_elem - 1) = Sigks[bidx];
  }
  /* τ² Ω prior blocks */
  for (std::size_t i = 0; i < idcs.size(); ++i)
    for (std::size_t j = 0; j < idcs.size(); ++j) {
      uword r1 = starts[i], r2 = starts[j];
      SigTot.submat(r1, r2,
                    r1 + idcs[i].n_elem - 1,
                    r2 + idcs[j].n_elem - 1)
        += tau2 * Omega.submat(idcs[i], idcs[j]);
    }
    
    mat L;  safe_chol(L, SigTot);
    double ld = 2.0 * arma::sum(log(L.diag()));
    
    rowvec bvec(big);
    for (std::size_t k = 0; k < idcs.size(); ++k) {
      rowvec brow = copy_row(b, obs_row[k]);
      uword st = starts[k];
      for (uword c = 0; c < idcs[k].n_elem; ++c)
        bvec[st + c] = brow[idcs[k][c]];
    }
    
    return dmvnrm1(bvec, L, ld);
}

vec deltis(const mat &b, const mat &s, const mat &Ce,
           double lambda, double alpha)
{
  vec out(b.n_rows, arma::fill::zeros);
  
  for (uword k = 0; k < b.n_rows; ++k) {
    rowvec srow = copy_row(s, k);
    uvec   idx  = finite_idx(srow);
    if (idx.is_empty()) continue;
    
    mat Sig, L;  make_sigma(srow, idx, Ce, Sig, L);
    double ld = 2.0 * arma::sum(log(L.diag()));
    
    rowvec x = subset(copy_row(b, k), idx);
    double ln0 = dmvnrm1(x, L, ld);
    
    mat   La  = std::sqrt(alpha) * L;
    double ldA = ld + idx.n_elem * std::log(alpha);
    double lna = dmvnrm1(x, La, ldA);
    
    double denom = lse(std::log(1.0 - lambda) + ln0,
                       std::log(lambda)       + lna);
    out[k] = std::exp(std::log(lambda) + lna - denom);
  }
  return out;
}

/* ------------------------------------------------------------------ */
/*                        one-SNP  E-step                              */
/* ------------------------------------------------------------------ */

struct SNPResult {
  double Rhat, ll, alph_num, alph_den, tau_num;
  vec    Ohat;
  rowvec mu_hat;
};

SNPResult E_step_one(const mat &b, const mat &s,
                     const mat &Ce, const mat &Om, const mat &OmInv,
                     double lambda, double alpha, double tau2, double p)
{
  const uword Q = b.n_cols, K = b.n_rows;
  
  double ll1 = ll_R1_j(b, s, Ce, Om, tau2);
  double ll0 = ll_R0_j(b, s, Ce, lambda, alpha);
  double l1  = std::log(p)       + ll1;
  double l0  = std::log(1.0 - p) + ll0;
  double Rhat = std::exp(l1 - lse(l1, l0));
  
  vec Ohat = deltis(b, s, Ce, lambda, alpha);
  
  mat Lambda = OmInv / tau2;
  rowvec tmp(Q, arma::fill::zeros);
  double alph_num = 0.0, alph_den = 0.0;
  
  for (uword k = 0; k < K; ++k) {
    rowvec srow = copy_row(s, k);
    uvec   idx  = finite_idx(srow);
    if (idx.is_empty()) { Ohat[k] = 0.0; continue; }
    
    mat Sig, L;  make_sigma(srow, idx, Ce, Sig, L);
    mat invS = arma::inv_sympd(Sig);
    
    rowvec b_obs = subset(copy_row(b, k), idx);
    double quad  = arma::as_scalar(b_obs * invS * b_obs.t());
    
    alph_num += (1.0 - Rhat) * Ohat[k] * quad;
    alph_den += (1.0 - Rhat) * Ohat[k];
    
    mat invS_full(Q, Q, arma::fill::zeros);
    invS_full.submat(idx, idx) = invS;
    Lambda += invS_full;
    
    rowvec contrib = b_obs * invS;
    for (uword c = 0; c < idx.n_elem; ++c)
      tmp[idx[c]] += contrib[c];
  }
  
  mat Vj     = arma::inv_sympd(Lambda);
  rowvec muj = tmp * Vj.t();
  double tau_comp = Rhat * arma::trace(OmInv * (Vj + muj.t() * muj));
  
  SNPResult res;
  res.Rhat = Rhat;   res.ll = lse(l1, l0);
  res.alph_num = alph_num;   res.alph_den = alph_den;
  res.tau_num  = tau_comp;   res.Ohat = Ohat;   res.mu_hat = muj;
  return res;
}

/* ------------------------------------------------------------------ */
/*                    full E-step  (serial + OMP)                      */
/* ------------------------------------------------------------------ */

// [[Rcpp::export]]
List E_step(const cube &B, const cube &S,
            const mat  &Ce, const mat  &Om, const mat  &OmInv,
            double lambda, double alpha, double tau2, double p)
{
  const uword M = B.n_slices, K = B.n_cols, Q = B.n_rows;
  
  vec  Rhat(M),  ll(M), alph_num(M), alph_den(M), tau_num(M);
  mat  Ohat(M, K);   mat mu_hat(M, Q);
  
  auto loop = [&](uword j) {
    SNPResult r = E_step_one(B.slice(j).t(), S.slice(j).t(),
                             Ce, Om, OmInv, lambda, alpha, tau2, p);
    Rhat[j]      = r.Rhat;    ll[j]       = r.ll;
    alph_num[j]  = r.alph_num; alph_den[j] = r.alph_den;
    tau_num[j]   = r.tau_num;
    for (uword k = 0; k < K; ++k) Ohat(j, k) = r.Ohat[k];
    for (uword q = 0; q < Q; ++q) mu_hat(j, q) = r.mu_hat[q];
  };
  
#ifdef _OPENMP
  if (omp_get_max_threads() > 1) {
#pragma omp parallel for schedule(static)
    for (uword j = 0; j < M; ++j) loop(j);
  } else
#endif
    for (uword j = 0; j < M; ++j) loop(j);
  
  return List::create(_["Rhat"]=Rhat,
                      _["Ohat"]=Ohat,
                      _["mu_hat"]=mu_hat,
                      _["alph_num"]=alph_num,
                      _["alph_den"]=alph_den,
                      _["tau_num"]=tau_num,
                      _["ll"]=ll);
}
