// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

// rowsum_Mat
vec rowsum_Mat(mat M);
RcppExport SEXP scEMundir_rowsum_Mat(SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(rowsum_Mat(M));
    return rcpp_result_gen;
END_RCPP
}
// colsum_Mat
vec colsum_Mat(mat M);
RcppExport SEXP scEMundir_colsum_Mat(SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(colsum_Mat(M));
    return rcpp_result_gen;
END_RCPP
}
// epan
float epan(float input);
RcppExport SEXP scEMundir_epan(SEXP inputSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< float >::type input(inputSEXP);
    rcpp_result_gen = Rcpp::wrap(epan(input));
    return rcpp_result_gen;
END_RCPP
}
// gamma_update_sc_undir
cube gamma_update_sc_undir(mat gamma, vec pi, vec theta, mat network, int N, int K);
RcppExport SEXP scEMundir_gamma_update_sc_undir(SEXP gammaSEXP, SEXP piSEXP, SEXP thetaSEXP, SEXP networkSEXP, SEXP NSEXP, SEXP KSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< vec >::type pi(piSEXP);
    Rcpp::traits::input_parameter< vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< mat >::type network(networkSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    rcpp_result_gen = Rcpp::wrap(gamma_update_sc_undir(gamma, pi, theta, network, N, K));
    return rcpp_result_gen;
END_RCPP
}
// grad_sc_EM_undir
mat grad_sc_EM_undir(vec theta, mat gamma, mat network, int N, int K);
RcppExport SEXP scEMundir_grad_sc_EM_undir(SEXP thetaSEXP, SEXP gammaSEXP, SEXP networkSEXP, SEXP NSEXP, SEXP KSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< mat >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< mat >::type network(networkSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    rcpp_result_gen = Rcpp::wrap(grad_sc_EM_undir(theta, gamma, network, N, K));
    return rcpp_result_gen;
END_RCPP
}
// hess_sc_EM_undir
mat hess_sc_EM_undir(vec theta, mat gamma, int N, int K);
RcppExport SEXP scEMundir_hess_sc_EM_undir(SEXP thetaSEXP, SEXP gammaSEXP, SEXP NSEXP, SEXP KSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< mat >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    rcpp_result_gen = Rcpp::wrap(hess_sc_EM_undir(theta, gamma, N, K));
    return rcpp_result_gen;
END_RCPP
}
// ELBO_conv_sc_EM_undir
float ELBO_conv_sc_EM_undir(mat gamma, vec pi, vec theta, mat network, int N, int K);
RcppExport SEXP scEMundir_ELBO_conv_sc_EM_undir(SEXP gammaSEXP, SEXP piSEXP, SEXP thetaSEXP, SEXP networkSEXP, SEXP NSEXP, SEXP KSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< vec >::type pi(piSEXP);
    Rcpp::traits::input_parameter< vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< mat >::type network(networkSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    rcpp_result_gen = Rcpp::wrap(ELBO_conv_sc_EM_undir(gamma, pi, theta, network, N, K));
    return rcpp_result_gen;
END_RCPP
}
