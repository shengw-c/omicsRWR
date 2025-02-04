#include <RcppArmadillo.h>
#include <RcppEigen.h>    // For efficient matrix operations
#include <progress.hpp>   // For progress bar

// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::SelfAdjointEigenSolver;

// [[Rcpp::export]]
NumericVector Cpp_rowSums(const NumericMatrix& x) {
  int nr = x.nrow(), nc = x.ncol();
  NumericVector ans(nr);
  for (int j = 0; j < nc; j++) {
    for (int i = 0; i < nr; i++) {
      ans[i] += x(i, j);
    }
  }
  return ans;
}

// [[Rcpp::export]]
NumericVector Cpp_colSums(const NumericMatrix& x) {
  int nr = x.nrow(), nc = x.ncol();
  NumericVector ans(nc);
  for (int i = 0; i < nr; i++) {
    for (int j = 0; j < nc; j++) {
      ans[i] += x(i, j);
    }
  }
  return ans;
}

// [[Rcpp::export]]
arma::mat Arma_Inv(const arma::mat & x) {
  return arma::inv(x);
}

// [[Rcpp::export]]
List rwr_cpp(NumericMatrix adj_mat,
             NumericVector init_scores,
             Nullable<NumericVector> int_scores = R_NilValue,
             double restart_prop = 0.15,
             int num_iter = 1000,
             double delta = 1e-6) {
  const Map<MatrixXd> AdjMat(as<Map<MatrixXd>>(adj_mat));
  const Map<VectorXd> Init_Scores(as<Map<VectorXd>>(init_scores));

  VectorXd Init_Scores_norm;
  if (std::abs(Init_Scores.sum() - 1.0) < 1e-10) {
    Init_Scores_norm = Init_Scores;
  } else {
    warning("Normalizing scores to sum as 1...");
    Init_Scores_norm = Init_Scores.array() / Init_Scores.sum();
  }

  Progress pb(num_iter, true);

  MatrixXd W_p = AdjMat * (1.0 - restart_prop);
  W_p.diagonal().setZero();

  VectorXd restart_term = restart_prop * Init_Scores_norm;

  VectorXd F;
  if (int_scores.isNotNull()) {
    warning("Using provided int scores to continue propagation...");
    F = as<VectorXd>(int_scores);
  } else {
    F = Init_Scores_norm;
  }

  double prevNorm = F.norm();
  bool converged = false;

  MatrixXd W_p_t = W_p.transpose();
  for (int i = 0; i < num_iter; ++i) {
    VectorXd newF = (W_p_t * F) + restart_term;
    double currentNorm = newF.norm();

    if (std::abs(currentNorm - prevNorm) < delta) {
      converged = true;
      Rcout << "Convergence at iteration " << ++i << std::endl;
      break;
    }

    F = newF;
    prevNorm = currentNorm;
    pb.increment();
  }

  if (!converged) {
    warning("RWR did not converge within the specified number of iterations.");
  }

  List result;
  result["vector"] = wrap(F);
  result["isconverged"] = converged;

  return result;
}
