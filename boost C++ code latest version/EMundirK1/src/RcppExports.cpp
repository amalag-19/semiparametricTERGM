// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

// rowsum_Mat
vec rowsum_Mat(mat M);
RcppExport SEXP EMundirK1_rowsum_Mat(SEXP MSEXP) {
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
RcppExport SEXP EMundirK1_colsum_Mat(SEXP MSEXP) {
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
RcppExport SEXP EMundirK1_epan(SEXP inputSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< float >::type input(inputSEXP);
    rcpp_result_gen = Rcpp::wrap(epan(input));
    return rcpp_result_gen;
END_RCPP
}
// grad_EM_undir_K1
vec grad_EM_undir_K1(float theta_u, cube network, int N, int K, int T_data);
RcppExport SEXP EMundirK1_grad_EM_undir_K1(SEXP theta_uSEXP, SEXP networkSEXP, SEXP NSEXP, SEXP KSEXP, SEXP T_dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< float >::type theta_u(theta_uSEXP);
    Rcpp::traits::input_parameter< cube >::type network(networkSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type T_data(T_dataSEXP);
    rcpp_result_gen = Rcpp::wrap(grad_EM_undir_K1(theta_u, network, N, K, T_data));
    return rcpp_result_gen;
END_RCPP
}
// hess_EM_undir_K1
float hess_EM_undir_K1(float theta_u, int N, int K, int T_data, float bandwidth, vec data_ids, vec grid_ids, int grid_id_index);
RcppExport SEXP EMundirK1_hess_EM_undir_K1(SEXP theta_uSEXP, SEXP NSEXP, SEXP KSEXP, SEXP T_dataSEXP, SEXP bandwidthSEXP, SEXP data_idsSEXP, SEXP grid_idsSEXP, SEXP grid_id_indexSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< float >::type theta_u(theta_uSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type T_data(T_dataSEXP);
    Rcpp::traits::input_parameter< float >::type bandwidth(bandwidthSEXP);
    Rcpp::traits::input_parameter< vec >::type data_ids(data_idsSEXP);
    Rcpp::traits::input_parameter< vec >::type grid_ids(grid_idsSEXP);
    Rcpp::traits::input_parameter< int >::type grid_id_index(grid_id_indexSEXP);
    rcpp_result_gen = Rcpp::wrap(hess_EM_undir_K1(theta_u, N, K, T_data, bandwidth, data_ids, grid_ids, grid_id_index));
    return rcpp_result_gen;
END_RCPP
}
// ELBO_conv_EM_undir_K1
float ELBO_conv_EM_undir_K1(float theta_u, cube network, int N, int K, int T_data, float bandwidth, vec data_ids, vec grid_ids, int grid_id_index);
RcppExport SEXP EMundirK1_ELBO_conv_EM_undir_K1(SEXP theta_uSEXP, SEXP networkSEXP, SEXP NSEXP, SEXP KSEXP, SEXP T_dataSEXP, SEXP bandwidthSEXP, SEXP data_idsSEXP, SEXP grid_idsSEXP, SEXP grid_id_indexSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< float >::type theta_u(theta_uSEXP);
    Rcpp::traits::input_parameter< cube >::type network(networkSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type T_data(T_dataSEXP);
    Rcpp::traits::input_parameter< float >::type bandwidth(bandwidthSEXP);
    Rcpp::traits::input_parameter< vec >::type data_ids(data_idsSEXP);
    Rcpp::traits::input_parameter< vec >::type grid_ids(grid_idsSEXP);
    Rcpp::traits::input_parameter< int >::type grid_id_index(grid_id_indexSEXP);
    rcpp_result_gen = Rcpp::wrap(ELBO_conv_EM_undir_K1(theta_u, network, N, K, T_data, bandwidth, data_ids, grid_ids, grid_id_index));
    return rcpp_result_gen;
END_RCPP
}
// grad_ELBO_K1_CV
vec grad_ELBO_K1_CV(float theta_u, vec gamma, cube network, int N, int K, int T_data, int test_node_set_len);
RcppExport SEXP EMundirK1_grad_ELBO_K1_CV(SEXP theta_uSEXP, SEXP gammaSEXP, SEXP networkSEXP, SEXP NSEXP, SEXP KSEXP, SEXP T_dataSEXP, SEXP test_node_set_lenSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< float >::type theta_u(theta_uSEXP);
    Rcpp::traits::input_parameter< vec >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< cube >::type network(networkSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type T_data(T_dataSEXP);
    Rcpp::traits::input_parameter< int >::type test_node_set_len(test_node_set_lenSEXP);
    rcpp_result_gen = Rcpp::wrap(grad_ELBO_K1_CV(theta_u, gamma, network, N, K, T_data, test_node_set_len));
    return rcpp_result_gen;
END_RCPP
}
// hess_ELBO_K1_CV
float hess_ELBO_K1_CV(float theta_u, vec gamma, cube network, int N, int K, int T_data, float bandwidth, vec data_ids, vec grid_ids, int grid_id_index, int test_node_set_len);
RcppExport SEXP EMundirK1_hess_ELBO_K1_CV(SEXP theta_uSEXP, SEXP gammaSEXP, SEXP networkSEXP, SEXP NSEXP, SEXP KSEXP, SEXP T_dataSEXP, SEXP bandwidthSEXP, SEXP data_idsSEXP, SEXP grid_idsSEXP, SEXP grid_id_indexSEXP, SEXP test_node_set_lenSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< float >::type theta_u(theta_uSEXP);
    Rcpp::traits::input_parameter< vec >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< cube >::type network(networkSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type T_data(T_dataSEXP);
    Rcpp::traits::input_parameter< float >::type bandwidth(bandwidthSEXP);
    Rcpp::traits::input_parameter< vec >::type data_ids(data_idsSEXP);
    Rcpp::traits::input_parameter< vec >::type grid_ids(grid_idsSEXP);
    Rcpp::traits::input_parameter< int >::type grid_id_index(grid_id_indexSEXP);
    Rcpp::traits::input_parameter< int >::type test_node_set_len(test_node_set_lenSEXP);
    rcpp_result_gen = Rcpp::wrap(hess_ELBO_K1_CV(theta_u, gamma, network, N, K, T_data, bandwidth, data_ids, grid_ids, grid_id_index, test_node_set_len));
    return rcpp_result_gen;
END_RCPP
}
// ELBO_conv_K1_CV
float ELBO_conv_K1_CV(vec gamma, float pi, float theta_u, cube network, int N, int K, int T_data, float bandwidth, vec data_ids, vec grid_ids, int grid_id_index, int test_node_set_len);
RcppExport SEXP EMundirK1_ELBO_conv_K1_CV(SEXP gammaSEXP, SEXP piSEXP, SEXP theta_uSEXP, SEXP networkSEXP, SEXP NSEXP, SEXP KSEXP, SEXP T_dataSEXP, SEXP bandwidthSEXP, SEXP data_idsSEXP, SEXP grid_idsSEXP, SEXP grid_id_indexSEXP, SEXP test_node_set_lenSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< vec >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< float >::type pi(piSEXP);
    Rcpp::traits::input_parameter< float >::type theta_u(theta_uSEXP);
    Rcpp::traits::input_parameter< cube >::type network(networkSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type T_data(T_dataSEXP);
    Rcpp::traits::input_parameter< float >::type bandwidth(bandwidthSEXP);
    Rcpp::traits::input_parameter< vec >::type data_ids(data_idsSEXP);
    Rcpp::traits::input_parameter< vec >::type grid_ids(grid_idsSEXP);
    Rcpp::traits::input_parameter< int >::type grid_id_index(grid_id_indexSEXP);
    Rcpp::traits::input_parameter< int >::type test_node_set_len(test_node_set_lenSEXP);
    rcpp_result_gen = Rcpp::wrap(ELBO_conv_K1_CV(gamma, pi, theta_u, network, N, K, T_data, bandwidth, data_ids, grid_ids, grid_id_index, test_node_set_len));
    return rcpp_result_gen;
END_RCPP
}
