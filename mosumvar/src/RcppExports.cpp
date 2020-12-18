// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// H_ik
arma::vec H_ik(arma::mat& x, int i, int k, int p, arma::mat Phi, arma::mat eps);
RcppExport SEXP _mosumvar_H_ik(SEXP xSEXP, SEXP iSEXP, SEXP kSEXP, SEXP pSEXP, SEXP PhiSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Phi(PhiSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(H_ik(x, i, k, p, Phi, eps));
    return rcpp_result_gen;
END_RCPP
}
// H_ik_univ
arma::vec H_ik_univ(arma::mat& x, int i, int k, int p, arma::mat Phi, arma::mat eps);
RcppExport SEXP _mosumvar_H_ik_univ(SEXP xSEXP, SEXP iSEXP, SEXP kSEXP, SEXP pSEXP, SEXP PhiSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Phi(PhiSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(H_ik_univ(x, i, k, p, Phi, eps));
    return rcpp_result_gen;
END_RCPP
}
// H_k
arma::vec H_k(arma::mat& x, int k, int p, arma::mat Phi, arma::mat eps);
RcppExport SEXP _mosumvar_H_k(SEXP xSEXP, SEXP kSEXP, SEXP pSEXP, SEXP PhiSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Phi(PhiSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(H_k(x, k, p, Phi, eps));
    return rcpp_result_gen;
END_RCPP
}
// H_k_univ
arma::vec H_k_univ(arma::mat& x, int k, int p, arma::field<arma::mat> PhiList, arma::mat eps);
RcppExport SEXP _mosumvar_H_k_univ(SEXP xSEXP, SEXP kSEXP, SEXP pSEXP, SEXP PhiListSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< arma::field<arma::mat> >::type PhiList(PhiListSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(H_k_univ(x, k, p, PhiList, eps));
    return rcpp_result_gen;
END_RCPP
}
// H_all
arma::mat H_all(arma::mat& x, int p, int G, arma::mat Phi, arma::mat eps);
RcppExport SEXP _mosumvar_H_all(SEXP xSEXP, SEXP pSEXP, SEXP GSEXP, SEXP PhiSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type G(GSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Phi(PhiSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(H_all(x, p, G, Phi, eps));
    return rcpp_result_gen;
END_RCPP
}
// H_all_univ
arma::mat H_all_univ(arma::mat& x, int p, int G, arma::field<arma::mat> PhiList, arma::mat eps);
RcppExport SEXP _mosumvar_H_all_univ(SEXP xSEXP, SEXP pSEXP, SEXP GSEXP, SEXP PhiListSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type G(GSEXP);
    Rcpp::traits::input_parameter< arma::field<arma::mat> >::type PhiList(PhiListSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(H_all_univ(x, p, G, PhiList, eps));
    return rcpp_result_gen;
END_RCPP
}
// DiagH
arma::sp_mat DiagH(arma::mat x, int k, int G, int p, arma::mat h_all);
RcppExport SEXP _mosumvar_DiagH(SEXP xSEXP, SEXP kSEXP, SEXP GSEXP, SEXP pSEXP, SEXP h_allSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type G(GSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type h_all(h_allSEXP);
    rcpp_result_gen = Rcpp::wrap(DiagH(x, k, G, p, h_all));
    return rcpp_result_gen;
END_RCPP
}
// FullH
arma::mat FullH(arma::mat x, int k, int G, arma::mat h_all);
RcppExport SEXP _mosumvar_FullH(SEXP xSEXP, SEXP kSEXP, SEXP GSEXP, SEXP h_allSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type G(GSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type h_all(h_allSEXP);
    rcpp_result_gen = Rcpp::wrap(FullH(x, k, G, h_all));
    return rcpp_result_gen;
END_RCPP
}
// getA
arma::vec getA(arma::mat x, int k, int G, int p, arma::mat eps, arma::mat h_all);
RcppExport SEXP _mosumvar_getA(SEXP xSEXP, SEXP kSEXP, SEXP GSEXP, SEXP pSEXP, SEXP epsSEXP, SEXP h_allSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type G(GSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type h_all(h_allSEXP);
    rcpp_result_gen = Rcpp::wrap(getA(x, k, G, p, eps, h_all));
    return rcpp_result_gen;
END_RCPP
}
// getsigma_iGlobal
double getsigma_iGlobal(arma::mat eps, int p, int i);
RcppExport SEXP _mosumvar_getsigma_iGlobal(SEXP epsSEXP, SEXP pSEXP, SEXP iSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    rcpp_result_gen = Rcpp::wrap(getsigma_iGlobal(eps, p, i));
    return rcpp_result_gen;
END_RCPP
}
// getsigma_dGlobal
arma::mat getsigma_dGlobal(arma::mat eps, int p);
RcppExport SEXP _mosumvar_getsigma_dGlobal(SEXP epsSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(getsigma_dGlobal(eps, p));
    return rcpp_result_gen;
END_RCPP
}
// getsigma_dLocal
arma::mat getsigma_dLocal(arma::mat eps, int k, int p, int G);
RcppExport SEXP _mosumvar_getsigma_dLocal(SEXP epsSEXP, SEXP kSEXP, SEXP pSEXP, SEXP GSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type G(GSEXP);
    rcpp_result_gen = Rcpp::wrap(getsigma_dLocal(eps, k, p, G));
    return rcpp_result_gen;
END_RCPP
}
// DiagC
arma::mat DiagC(arma::mat x, int p, arma::mat sigma_d, int k, int G);
RcppExport SEXP _mosumvar_DiagC(SEXP xSEXP, SEXP pSEXP, SEXP sigma_dSEXP, SEXP kSEXP, SEXP GSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma_d(sigma_dSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type G(GSEXP);
    rcpp_result_gen = Rcpp::wrap(DiagC(x, p, sigma_d, k, G));
    return rcpp_result_gen;
END_RCPP
}
// DiagC_univ
arma::mat DiagC_univ(arma::mat x, int p, arma::mat sigma_d, int k, int G);
RcppExport SEXP _mosumvar_DiagC_univ(SEXP xSEXP, SEXP pSEXP, SEXP sigma_dSEXP, SEXP kSEXP, SEXP GSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma_d(sigma_dSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type G(GSEXP);
    rcpp_result_gen = Rcpp::wrap(DiagC_univ(x, p, sigma_d, k, G));
    return rcpp_result_gen;
END_RCPP
}
// Tkn
double Tkn(arma::mat x, int k, int p, int G, arma::mat Phi, arma::mat eps, arma::mat h_all, String estim, String var_estim, arma::mat sgd, bool univariate);
RcppExport SEXP _mosumvar_Tkn(SEXP xSEXP, SEXP kSEXP, SEXP pSEXP, SEXP GSEXP, SEXP PhiSEXP, SEXP epsSEXP, SEXP h_allSEXP, SEXP estimSEXP, SEXP var_estimSEXP, SEXP sgdSEXP, SEXP univariateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type G(GSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Phi(PhiSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type h_all(h_allSEXP);
    Rcpp::traits::input_parameter< String >::type estim(estimSEXP);
    Rcpp::traits::input_parameter< String >::type var_estim(var_estimSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sgd(sgdSEXP);
    Rcpp::traits::input_parameter< bool >::type univariate(univariateSEXP);
    rcpp_result_gen = Rcpp::wrap(Tkn(x, k, p, G, Phi, eps, h_all, estim, var_estim, sgd, univariate));
    return rcpp_result_gen;
END_RCPP
}
// T
arma::vec T(arma::mat x, int p, int G, arma::mat Phi, arma::mat eps, arma::field<arma::mat> PhiList, String estim, String var_estim, bool univariate);
RcppExport SEXP _mosumvar_T(SEXP xSEXP, SEXP pSEXP, SEXP GSEXP, SEXP PhiSEXP, SEXP epsSEXP, SEXP PhiListSEXP, SEXP estimSEXP, SEXP var_estimSEXP, SEXP univariateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type G(GSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Phi(PhiSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< arma::field<arma::mat> >::type PhiList(PhiListSEXP);
    Rcpp::traits::input_parameter< String >::type estim(estimSEXP);
    Rcpp::traits::input_parameter< String >::type var_estim(var_estimSEXP);
    Rcpp::traits::input_parameter< bool >::type univariate(univariateSEXP);
    rcpp_result_gen = Rcpp::wrap(T(x, p, G, Phi, eps, PhiList, estim, var_estim, univariate));
    return rcpp_result_gen;
END_RCPP
}
// test_Score
List test_Score(arma::mat x, int p, int G, arma::mat Phi, arma::mat eps, double alpha, String estim);
RcppExport SEXP _mosumvar_test_Score(SEXP xSEXP, SEXP pSEXP, SEXP GSEXP, SEXP PhiSEXP, SEXP epsSEXP, SEXP alphaSEXP, SEXP estimSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type G(GSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Phi(PhiSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< String >::type estim(estimSEXP);
    rcpp_result_gen = Rcpp::wrap(test_Score(x, p, G, Phi, eps, alpha, estim));
    return rcpp_result_gen;
END_RCPP
}
// MFA_Score
List MFA_Score(arma::mat x, int p, arma::vec Gset, arma::mat Phi, arma::mat eps, String estim, double alpha);
RcppExport SEXP _mosumvar_MFA_Score(SEXP xSEXP, SEXP pSEXP, SEXP GsetSEXP, SEXP PhiSEXP, SEXP epsSEXP, SEXP estimSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Gset(GsetSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Phi(PhiSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< String >::type estim(estimSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(MFA_Score(x, p, Gset, Phi, eps, estim, alpha));
    return rcpp_result_gen;
END_RCPP
}
// VAR_sim
arma::mat VAR_sim(int n, arma::vec mu, arma::mat Sigma, arma::field<arma::mat> coeffs, std::string error_dist, arma::mat P1, arma::mat Q1, int df);
RcppExport SEXP _mosumvar_VAR_sim(SEXP nSEXP, SEXP muSEXP, SEXP SigmaSEXP, SEXP coeffsSEXP, SEXP error_distSEXP, SEXP P1SEXP, SEXP Q1SEXP, SEXP dfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Sigma(SigmaSEXP);
    Rcpp::traits::input_parameter< arma::field<arma::mat> >::type coeffs(coeffsSEXP);
    Rcpp::traits::input_parameter< std::string >::type error_dist(error_distSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type P1(P1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Q1(Q1SEXP);
    Rcpp::traits::input_parameter< int >::type df(dfSEXP);
    rcpp_result_gen = Rcpp::wrap(VAR_sim(n, mu, Sigma, coeffs, error_dist, P1, Q1, df));
    return rcpp_result_gen;
END_RCPP
}
// H_ik_Wald
arma::vec H_ik_Wald(arma::mat& x, int i, int k, int p, arma::vec a);
RcppExport SEXP _mosumvar_H_ik_Wald(SEXP xSEXP, SEXP iSEXP, SEXP kSEXP, SEXP pSEXP, SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type a(aSEXP);
    rcpp_result_gen = Rcpp::wrap(H_ik_Wald(x, i, k, p, a));
    return rcpp_result_gen;
END_RCPP
}
// H_k_Wald
arma::vec H_k_Wald(arma::mat& x, int k, int p, arma::vec a);
RcppExport SEXP _mosumvar_H_k_Wald(SEXP xSEXP, SEXP kSEXP, SEXP pSEXP, SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type a(aSEXP);
    rcpp_result_gen = Rcpp::wrap(H_k_Wald(x, k, p, a));
    return rcpp_result_gen;
END_RCPP
}
// H_l_u
arma::mat H_l_u(arma::mat& x, int p, int l, int u, arma::vec a);
RcppExport SEXP _mosumvar_H_l_u(SEXP xSEXP, SEXP pSEXP, SEXP lSEXP, SEXP uSEXP, SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type l(lSEXP);
    Rcpp::traits::input_parameter< int >::type u(uSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type a(aSEXP);
    rcpp_result_gen = Rcpp::wrap(H_l_u(x, p, l, u, a));
    return rcpp_result_gen;
END_RCPP
}
// write_rows2
arma::mat write_rows2(Rcpp::List data, int nrows, int ncols);
RcppExport SEXP _mosumvar_write_rows2(SEXP dataSEXP, SEXP nrowsSEXP, SEXP ncolsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type nrows(nrowsSEXP);
    Rcpp::traits::input_parameter< int >::type ncols(ncolsSEXP);
    rcpp_result_gen = Rcpp::wrap(write_rows2(data, nrows, ncols));
    return rcpp_result_gen;
END_RCPP
}
// a_lu_i
arma::vec a_lu_i(arma::mat& x, int i, int p, int l, int u);
RcppExport SEXP _mosumvar_a_lu_i(SEXP xSEXP, SEXP iSEXP, SEXP pSEXP, SEXP lSEXP, SEXP uSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type l(lSEXP);
    Rcpp::traits::input_parameter< int >::type u(uSEXP);
    rcpp_result_gen = Rcpp::wrap(a_lu_i(x, i, p, l, u));
    return rcpp_result_gen;
END_RCPP
}
// make_a_lu
arma::vec make_a_lu(arma::mat& x, int p, int l, int u);
RcppExport SEXP _mosumvar_make_a_lu(SEXP xSEXP, SEXP pSEXP, SEXP lSEXP, SEXP uSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type l(lSEXP);
    Rcpp::traits::input_parameter< int >::type u(uSEXP);
    rcpp_result_gen = Rcpp::wrap(make_a_lu(x, p, l, u));
    return rcpp_result_gen;
END_RCPP
}
// V_nk
arma::sp_mat V_nk(arma::mat x, int p, int l, int u);
RcppExport SEXP _mosumvar_V_nk(SEXP xSEXP, SEXP pSEXP, SEXP lSEXP, SEXP uSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type l(lSEXP);
    Rcpp::traits::input_parameter< int >::type u(uSEXP);
    rcpp_result_gen = Rcpp::wrap(V_nk(x, p, l, u));
    return rcpp_result_gen;
END_RCPP
}
// sigma_i_k
arma::vec sigma_i_k(arma::mat x, int i, int k, int G, int p, arma::vec a_upper);
RcppExport SEXP _mosumvar_sigma_i_k(SEXP xSEXP, SEXP iSEXP, SEXP kSEXP, SEXP GSEXP, SEXP pSEXP, SEXP a_upperSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type G(GSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type a_upper(a_upperSEXP);
    rcpp_result_gen = Rcpp::wrap(sigma_i_k(x, i, k, G, p, a_upper));
    return rcpp_result_gen;
END_RCPP
}
// sigma_d_k
arma::mat sigma_d_k(arma::mat x, int k, int G, int p, arma::vec a_upper, arma::vec a_lower);
RcppExport SEXP _mosumvar_sigma_d_k(SEXP xSEXP, SEXP kSEXP, SEXP GSEXP, SEXP pSEXP, SEXP a_upperSEXP, SEXP a_lowerSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type G(GSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type a_upper(a_upperSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type a_lower(a_lowerSEXP);
    rcpp_result_gen = Rcpp::wrap(sigma_d_k(x, k, G, p, a_upper, a_lower));
    return rcpp_result_gen;
END_RCPP
}
// diagH_Wald
arma::sp_mat diagH_Wald(arma::mat x, int G, int p, arma::mat H_l, arma::mat H_u);
RcppExport SEXP _mosumvar_diagH_Wald(SEXP xSEXP, SEXP GSEXP, SEXP pSEXP, SEXP H_lSEXP, SEXP H_uSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type G(GSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type H_l(H_lSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type H_u(H_uSEXP);
    rcpp_result_gen = Rcpp::wrap(diagH_Wald(x, G, p, H_l, H_u));
    return rcpp_result_gen;
END_RCPP
}
// FullH_Wald
arma::mat FullH_Wald(arma::mat x, int G, arma::mat H_l, arma::mat H_u);
RcppExport SEXP _mosumvar_FullH_Wald(SEXP xSEXP, SEXP GSEXP, SEXP H_lSEXP, SEXP H_uSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type G(GSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type H_l(H_lSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type H_u(H_uSEXP);
    rcpp_result_gen = Rcpp::wrap(FullH_Wald(x, G, H_l, H_u));
    return rcpp_result_gen;
END_RCPP
}
// DiagC_Wald
arma::mat DiagC_Wald(arma::mat x, int p, arma::mat sigma_d, int k, int G);
RcppExport SEXP _mosumvar_DiagC_Wald(SEXP xSEXP, SEXP pSEXP, SEXP sigma_dSEXP, SEXP kSEXP, SEXP GSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma_d(sigma_dSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type G(GSEXP);
    rcpp_result_gen = Rcpp::wrap(DiagC_Wald(x, p, sigma_d, k, G));
    return rcpp_result_gen;
END_RCPP
}
// Wkn
double Wkn(arma::mat x, int p, int k, int G, String estim);
RcppExport SEXP _mosumvar_Wkn(SEXP xSEXP, SEXP pSEXP, SEXP kSEXP, SEXP GSEXP, SEXP estimSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type G(GSEXP);
    Rcpp::traits::input_parameter< String >::type estim(estimSEXP);
    rcpp_result_gen = Rcpp::wrap(Wkn(x, p, k, G, estim));
    return rcpp_result_gen;
END_RCPP
}
// W
arma::vec W(arma::mat x, int p, int G, String estim);
RcppExport SEXP _mosumvar_W(SEXP xSEXP, SEXP pSEXP, SEXP GSEXP, SEXP estimSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type G(GSEXP);
    Rcpp::traits::input_parameter< String >::type estim(estimSEXP);
    rcpp_result_gen = Rcpp::wrap(W(x, p, G, estim));
    return rcpp_result_gen;
END_RCPP
}
// cps_Wald
arma::vec cps_Wald(arma::vec Wn, double D_n, int G, double nu);
RcppExport SEXP _mosumvar_cps_Wald(SEXP WnSEXP, SEXP D_nSEXP, SEXP GSEXP, SEXP nuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type Wn(WnSEXP);
    Rcpp::traits::input_parameter< double >::type D_n(D_nSEXP);
    Rcpp::traits::input_parameter< int >::type G(GSEXP);
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    rcpp_result_gen = Rcpp::wrap(cps_Wald(Wn, D_n, G, nu));
    return rcpp_result_gen;
END_RCPP
}
// test_Wald
List test_Wald(arma::mat x, int p, int G, double alpha, String estim);
RcppExport SEXP _mosumvar_test_Wald(SEXP xSEXP, SEXP pSEXP, SEXP GSEXP, SEXP alphaSEXP, SEXP estimSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type G(GSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< String >::type estim(estimSEXP);
    rcpp_result_gen = Rcpp::wrap(test_Wald(x, p, G, alpha, estim));
    return rcpp_result_gen;
END_RCPP
}
// sim_data
arma::mat sim_data(List pars, int n, int d, double sd);
RcppExport SEXP _mosumvar_sim_data(SEXP parsSEXP, SEXP nSEXP, SEXP dSEXP, SEXP sdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< double >::type sd(sdSEXP);
    rcpp_result_gen = Rcpp::wrap(sim_data(pars, n, d, sd));
    return rcpp_result_gen;
END_RCPP
}
// MFA_Wald
List MFA_Wald(arma::mat x, int p, arma::vec Gset, String estim, double alpha);
RcppExport SEXP _mosumvar_MFA_Wald(SEXP xSEXP, SEXP pSEXP, SEXP GsetSEXP, SEXP estimSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Gset(GsetSEXP);
    Rcpp::traits::input_parameter< String >::type estim(estimSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(MFA_Wald(x, p, Gset, estim, alpha));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_mosumvar_H_ik", (DL_FUNC) &_mosumvar_H_ik, 6},
    {"_mosumvar_H_ik_univ", (DL_FUNC) &_mosumvar_H_ik_univ, 6},
    {"_mosumvar_H_k", (DL_FUNC) &_mosumvar_H_k, 5},
    {"_mosumvar_H_k_univ", (DL_FUNC) &_mosumvar_H_k_univ, 5},
    {"_mosumvar_H_all", (DL_FUNC) &_mosumvar_H_all, 5},
    {"_mosumvar_H_all_univ", (DL_FUNC) &_mosumvar_H_all_univ, 5},
    {"_mosumvar_DiagH", (DL_FUNC) &_mosumvar_DiagH, 5},
    {"_mosumvar_FullH", (DL_FUNC) &_mosumvar_FullH, 4},
    {"_mosumvar_getA", (DL_FUNC) &_mosumvar_getA, 6},
    {"_mosumvar_getsigma_iGlobal", (DL_FUNC) &_mosumvar_getsigma_iGlobal, 3},
    {"_mosumvar_getsigma_dGlobal", (DL_FUNC) &_mosumvar_getsigma_dGlobal, 2},
    {"_mosumvar_getsigma_dLocal", (DL_FUNC) &_mosumvar_getsigma_dLocal, 4},
    {"_mosumvar_DiagC", (DL_FUNC) &_mosumvar_DiagC, 5},
    {"_mosumvar_DiagC_univ", (DL_FUNC) &_mosumvar_DiagC_univ, 5},
    {"_mosumvar_Tkn", (DL_FUNC) &_mosumvar_Tkn, 11},
    {"_mosumvar_T", (DL_FUNC) &_mosumvar_T, 9},
    {"_mosumvar_test_Score", (DL_FUNC) &_mosumvar_test_Score, 7},
    {"_mosumvar_MFA_Score", (DL_FUNC) &_mosumvar_MFA_Score, 7},
    {"_mosumvar_VAR_sim", (DL_FUNC) &_mosumvar_VAR_sim, 8},
    {"_mosumvar_H_ik_Wald", (DL_FUNC) &_mosumvar_H_ik_Wald, 5},
    {"_mosumvar_H_k_Wald", (DL_FUNC) &_mosumvar_H_k_Wald, 4},
    {"_mosumvar_H_l_u", (DL_FUNC) &_mosumvar_H_l_u, 5},
    {"_mosumvar_write_rows2", (DL_FUNC) &_mosumvar_write_rows2, 3},
    {"_mosumvar_a_lu_i", (DL_FUNC) &_mosumvar_a_lu_i, 5},
    {"_mosumvar_make_a_lu", (DL_FUNC) &_mosumvar_make_a_lu, 4},
    {"_mosumvar_V_nk", (DL_FUNC) &_mosumvar_V_nk, 4},
    {"_mosumvar_sigma_i_k", (DL_FUNC) &_mosumvar_sigma_i_k, 6},
    {"_mosumvar_sigma_d_k", (DL_FUNC) &_mosumvar_sigma_d_k, 6},
    {"_mosumvar_diagH_Wald", (DL_FUNC) &_mosumvar_diagH_Wald, 5},
    {"_mosumvar_FullH_Wald", (DL_FUNC) &_mosumvar_FullH_Wald, 4},
    {"_mosumvar_DiagC_Wald", (DL_FUNC) &_mosumvar_DiagC_Wald, 5},
    {"_mosumvar_Wkn", (DL_FUNC) &_mosumvar_Wkn, 5},
    {"_mosumvar_W", (DL_FUNC) &_mosumvar_W, 4},
    {"_mosumvar_cps_Wald", (DL_FUNC) &_mosumvar_cps_Wald, 4},
    {"_mosumvar_test_Wald", (DL_FUNC) &_mosumvar_test_Wald, 5},
    {"_mosumvar_sim_data", (DL_FUNC) &_mosumvar_sim_data, 4},
    {"_mosumvar_MFA_Wald", (DL_FUNC) &_mosumvar_MFA_Wald, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_mosumvar(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
