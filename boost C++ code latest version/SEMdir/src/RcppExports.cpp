// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

// rowsum_Mat
vec rowsum_Mat(mat M);
RcppExport SEXP SEMdir_rowsum_Mat(SEXP MSEXP) {
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
RcppExport SEXP SEMdir_colsum_Mat(SEXP MSEXP) {
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
RcppExport SEXP SEMdir_epan(SEXP inputSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< float >::type input(inputSEXP);
    rcpp_result_gen = Rcpp::wrap(epan(input));
    return rcpp_result_gen;
END_RCPP
}
// gamma_update_dir
cube gamma_update_dir(mat gamma, vec pi, cube theta, cube network, int N, int K, int T_grid, vec grid_ids);
RcppExport SEXP SEMdir_gamma_update_dir(SEXP gammaSEXP, SEXP piSEXP, SEXP thetaSEXP, SEXP networkSEXP, SEXP NSEXP, SEXP KSEXP, SEXP T_gridSEXP, SEXP grid_idsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< vec >::type pi(piSEXP);
    Rcpp::traits::input_parameter< cube >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< cube >::type network(networkSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type T_grid(T_gridSEXP);
    Rcpp::traits::input_parameter< vec >::type grid_ids(grid_idsSEXP);
    rcpp_result_gen = Rcpp::wrap(gamma_update_dir(gamma, pi, theta, network, N, K, T_grid, grid_ids));
    return rcpp_result_gen;
END_RCPP
}
// grad_SEM_dir_oe
mat grad_SEM_dir_oe(mat theta_u, cube network, int N, int K, vec nonzero_ids, int nonzero_ids_len, vec cluster_ids);
RcppExport SEXP SEMdir_grad_SEM_dir_oe(SEXP theta_uSEXP, SEXP networkSEXP, SEXP NSEXP, SEXP KSEXP, SEXP nonzero_idsSEXP, SEXP nonzero_ids_lenSEXP, SEXP cluster_idsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat >::type theta_u(theta_uSEXP);
    Rcpp::traits::input_parameter< cube >::type network(networkSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< vec >::type nonzero_ids(nonzero_idsSEXP);
    Rcpp::traits::input_parameter< int >::type nonzero_ids_len(nonzero_ids_lenSEXP);
    Rcpp::traits::input_parameter< vec >::type cluster_ids(cluster_idsSEXP);
    rcpp_result_gen = Rcpp::wrap(grad_SEM_dir_oe(theta_u, network, N, K, nonzero_ids, nonzero_ids_len, cluster_ids));
    return rcpp_result_gen;
END_RCPP
}
// grad_SEM_dir_re
mat grad_SEM_dir_re(mat theta_u, cube network, int N, int K, vec nonzero_ids, int nonzero_ids_len, vec cluster_ids);
RcppExport SEXP SEMdir_grad_SEM_dir_re(SEXP theta_uSEXP, SEXP networkSEXP, SEXP NSEXP, SEXP KSEXP, SEXP nonzero_idsSEXP, SEXP nonzero_ids_lenSEXP, SEXP cluster_idsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat >::type theta_u(theta_uSEXP);
    Rcpp::traits::input_parameter< cube >::type network(networkSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< vec >::type nonzero_ids(nonzero_idsSEXP);
    Rcpp::traits::input_parameter< int >::type nonzero_ids_len(nonzero_ids_lenSEXP);
    Rcpp::traits::input_parameter< vec >::type cluster_ids(cluster_idsSEXP);
    rcpp_result_gen = Rcpp::wrap(grad_SEM_dir_re(theta_u, network, N, K, nonzero_ids, nonzero_ids_len, cluster_ids));
    return rcpp_result_gen;
END_RCPP
}
// hess_SEM_dir_oe
mat hess_SEM_dir_oe(mat theta_u, int N, int K, int T_data, float bandwidth, vec data_ids, vec grid_ids, int grid_id_index, vec cluster_ids);
RcppExport SEXP SEMdir_hess_SEM_dir_oe(SEXP theta_uSEXP, SEXP NSEXP, SEXP KSEXP, SEXP T_dataSEXP, SEXP bandwidthSEXP, SEXP data_idsSEXP, SEXP grid_idsSEXP, SEXP grid_id_indexSEXP, SEXP cluster_idsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat >::type theta_u(theta_uSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type T_data(T_dataSEXP);
    Rcpp::traits::input_parameter< float >::type bandwidth(bandwidthSEXP);
    Rcpp::traits::input_parameter< vec >::type data_ids(data_idsSEXP);
    Rcpp::traits::input_parameter< vec >::type grid_ids(grid_idsSEXP);
    Rcpp::traits::input_parameter< int >::type grid_id_index(grid_id_indexSEXP);
    Rcpp::traits::input_parameter< vec >::type cluster_ids(cluster_idsSEXP);
    rcpp_result_gen = Rcpp::wrap(hess_SEM_dir_oe(theta_u, N, K, T_data, bandwidth, data_ids, grid_ids, grid_id_index, cluster_ids));
    return rcpp_result_gen;
END_RCPP
}
// hess_SEM_dir_re
mat hess_SEM_dir_re(mat theta_u, int N, int K, int T_data, float bandwidth, vec data_ids, vec grid_ids, int grid_id_index, vec cluster_ids);
RcppExport SEXP SEMdir_hess_SEM_dir_re(SEXP theta_uSEXP, SEXP NSEXP, SEXP KSEXP, SEXP T_dataSEXP, SEXP bandwidthSEXP, SEXP data_idsSEXP, SEXP grid_idsSEXP, SEXP grid_id_indexSEXP, SEXP cluster_idsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat >::type theta_u(theta_uSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type T_data(T_dataSEXP);
    Rcpp::traits::input_parameter< float >::type bandwidth(bandwidthSEXP);
    Rcpp::traits::input_parameter< vec >::type data_ids(data_idsSEXP);
    Rcpp::traits::input_parameter< vec >::type grid_ids(grid_idsSEXP);
    Rcpp::traits::input_parameter< int >::type grid_id_index(grid_id_indexSEXP);
    Rcpp::traits::input_parameter< vec >::type cluster_ids(cluster_idsSEXP);
    rcpp_result_gen = Rcpp::wrap(hess_SEM_dir_re(theta_u, N, K, T_data, bandwidth, data_ids, grid_ids, grid_id_index, cluster_ids));
    return rcpp_result_gen;
END_RCPP
}
// hess_SEM_dir_oe_re
mat hess_SEM_dir_oe_re(mat theta_u, int N, int K, int T_data, float bandwidth, vec data_ids, vec grid_ids, int grid_id_index, vec cluster_ids);
RcppExport SEXP SEMdir_hess_SEM_dir_oe_re(SEXP theta_uSEXP, SEXP NSEXP, SEXP KSEXP, SEXP T_dataSEXP, SEXP bandwidthSEXP, SEXP data_idsSEXP, SEXP grid_idsSEXP, SEXP grid_id_indexSEXP, SEXP cluster_idsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat >::type theta_u(theta_uSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type T_data(T_dataSEXP);
    Rcpp::traits::input_parameter< float >::type bandwidth(bandwidthSEXP);
    Rcpp::traits::input_parameter< vec >::type data_ids(data_idsSEXP);
    Rcpp::traits::input_parameter< vec >::type grid_ids(grid_idsSEXP);
    Rcpp::traits::input_parameter< int >::type grid_id_index(grid_id_indexSEXP);
    Rcpp::traits::input_parameter< vec >::type cluster_ids(cluster_idsSEXP);
    rcpp_result_gen = Rcpp::wrap(hess_SEM_dir_oe_re(theta_u, N, K, T_data, bandwidth, data_ids, grid_ids, grid_id_index, cluster_ids));
    return rcpp_result_gen;
END_RCPP
}
// ELBO_conv_SEM_dir
float ELBO_conv_SEM_dir(mat gamma, vec pi, mat theta_u, cube network, int N, int K, int T_data, float bandwidth, vec data_ids, vec grid_ids, int grid_id_index, vec cluster_ids);
RcppExport SEXP SEMdir_ELBO_conv_SEM_dir(SEXP gammaSEXP, SEXP piSEXP, SEXP theta_uSEXP, SEXP networkSEXP, SEXP NSEXP, SEXP KSEXP, SEXP T_dataSEXP, SEXP bandwidthSEXP, SEXP data_idsSEXP, SEXP grid_idsSEXP, SEXP grid_id_indexSEXP, SEXP cluster_idsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< vec >::type pi(piSEXP);
    Rcpp::traits::input_parameter< mat >::type theta_u(theta_uSEXP);
    Rcpp::traits::input_parameter< cube >::type network(networkSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type T_data(T_dataSEXP);
    Rcpp::traits::input_parameter< float >::type bandwidth(bandwidthSEXP);
    Rcpp::traits::input_parameter< vec >::type data_ids(data_idsSEXP);
    Rcpp::traits::input_parameter< vec >::type grid_ids(grid_idsSEXP);
    Rcpp::traits::input_parameter< int >::type grid_id_index(grid_id_indexSEXP);
    Rcpp::traits::input_parameter< vec >::type cluster_ids(cluster_idsSEXP);
    rcpp_result_gen = Rcpp::wrap(ELBO_conv_SEM_dir(gamma, pi, theta_u, network, N, K, T_data, bandwidth, data_ids, grid_ids, grid_id_index, cluster_ids));
    return rcpp_result_gen;
END_RCPP
}
// VK_dir_oe
mat VK_dir_oe(cube theta, cube network, int N, int K, int T_data, vec cluster_ids);
RcppExport SEXP SEMdir_VK_dir_oe(SEXP thetaSEXP, SEXP networkSEXP, SEXP NSEXP, SEXP KSEXP, SEXP T_dataSEXP, SEXP cluster_idsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< cube >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< cube >::type network(networkSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type T_data(T_dataSEXP);
    Rcpp::traits::input_parameter< vec >::type cluster_ids(cluster_idsSEXP);
    rcpp_result_gen = Rcpp::wrap(VK_dir_oe(theta, network, N, K, T_data, cluster_ids));
    return rcpp_result_gen;
END_RCPP
}
// VK_dir_re
mat VK_dir_re(cube theta, cube network, int N, int K, int T_data, vec cluster_ids);
RcppExport SEXP SEMdir_VK_dir_re(SEXP thetaSEXP, SEXP networkSEXP, SEXP NSEXP, SEXP KSEXP, SEXP T_dataSEXP, SEXP cluster_idsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< cube >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< cube >::type network(networkSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type T_data(T_dataSEXP);
    Rcpp::traits::input_parameter< vec >::type cluster_ids(cluster_idsSEXP);
    rcpp_result_gen = Rcpp::wrap(VK_dir_re(theta, network, N, K, T_data, cluster_ids));
    return rcpp_result_gen;
END_RCPP
}
// HK_dir_oe
mat HK_dir_oe(cube theta, int N, int K, int T_data, vec cluster_ids);
RcppExport SEXP SEMdir_HK_dir_oe(SEXP thetaSEXP, SEXP NSEXP, SEXP KSEXP, SEXP T_dataSEXP, SEXP cluster_idsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< cube >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type T_data(T_dataSEXP);
    Rcpp::traits::input_parameter< vec >::type cluster_ids(cluster_idsSEXP);
    rcpp_result_gen = Rcpp::wrap(HK_dir_oe(theta, N, K, T_data, cluster_ids));
    return rcpp_result_gen;
END_RCPP
}
// HK_dir_re
mat HK_dir_re(cube theta, int N, int K, int T_data, vec cluster_ids);
RcppExport SEXP SEMdir_HK_dir_re(SEXP thetaSEXP, SEXP NSEXP, SEXP KSEXP, SEXP T_dataSEXP, SEXP cluster_idsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< cube >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type T_data(T_dataSEXP);
    Rcpp::traits::input_parameter< vec >::type cluster_ids(cluster_idsSEXP);
    rcpp_result_gen = Rcpp::wrap(HK_dir_re(theta, N, K, T_data, cluster_ids));
    return rcpp_result_gen;
END_RCPP
}
// HK_dir_oe_re
mat HK_dir_oe_re(cube theta, int N, int K, int T_data, vec cluster_ids);
RcppExport SEXP SEMdir_HK_dir_oe_re(SEXP thetaSEXP, SEXP NSEXP, SEXP KSEXP, SEXP T_dataSEXP, SEXP cluster_idsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< cube >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type T_data(T_dataSEXP);
    Rcpp::traits::input_parameter< vec >::type cluster_ids(cluster_idsSEXP);
    rcpp_result_gen = Rcpp::wrap(HK_dir_oe_re(theta, N, K, T_data, cluster_ids));
    return rcpp_result_gen;
END_RCPP
}
