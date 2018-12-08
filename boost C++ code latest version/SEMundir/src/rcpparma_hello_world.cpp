// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
using namespace arma;
using namespace std;

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R

// [[Rcpp::export]]
vec rowsum_Mat(mat M) {
    int nr=M.n_rows;
    vec out(nr);
    for(int i=0;i<nr;i++){
        out(i)=sum(M.row(i));
    }
    return out;
}

// [[Rcpp::export]]
vec colsum_Mat(mat M) {
    int nc=M.n_cols;
    vec out(nc);
    for(int i=0;i<nc;i++){
        out(i)=sum(M.col(i));
    }
    return out;
}

// [[Rcpp::export]]
float epan(float input){
    float output;
    if(abs(input)<=1){
        output=0.75*(1-pow(input,2));
    }
    else{
        output=0;
    }
    return output;
}

// [[Rcpp::export]]
cube gamma_update_undir(mat gamma, vec pi, mat theta, cube network, int N, int K, int T_grid, vec grid_ids){
    cube quad_lin_coeff(N,K,2);
    for(int i = 0; i < N; i++){
        if(i!=(N-1)){
            for(int k = 0; k < K; k++){
                float t1=0;
                for (int t=0; t<T_grid; t++){
                    for(int j = i+1; j < N; j++){
                        for(int l = 0; l < K; l++){
                            float exp_val=exp(theta(t,k)+theta(t,l));
                            t1+=(gamma(j,l)/(2*gamma(i,k)))*((network(i,j,(grid_ids(t)-1))*(theta(t,k)+theta(t,l)))-log(1+exp_val));
                        }
                    }
                }
                quad_lin_coeff(i,k,0)=t1-(T_grid/gamma(i,k));
                quad_lin_coeff(i,k,1)=T_grid*(log(pi(k))-log(gamma(i,k))+1);
            }
        } else if(i==(N-1)){
            for(int k = 0; k < K; k++){
                quad_lin_coeff(i,k,0)=-(T_grid/gamma((N-1),k));
                quad_lin_coeff(i,k,1)=T_grid*(log(pi(k))-log(gamma((N-1),k))+1);
            }
        }
    }
    return quad_lin_coeff;
}


// [[Rcpp::export]]
mat grad_SEM_undir(vec theta_u, cube network, int N, int K, vec nonzero_ids, int nonzero_ids_len, vec cluster_ids){
    mat grad_vector(nonzero_ids_len,K);
    for (int t=0; t<nonzero_ids_len; t++){
        mat grad_mat(K,K,fill::zeros);
        for(int i = 0; i < (N-1); i++){
            for(int j = i+1; j < N; j++){
                float exp_val=exp(theta_u(cluster_ids(i)-1)+theta_u(cluster_ids(j)-1));
                grad_mat(cluster_ids(i)-1,cluster_ids(j)-1)+=(network(i,j,(nonzero_ids(t)-1))-(exp_val/(1+exp_val)));
            }
        }
        vec rsum=rowsum_Mat(grad_mat);
        vec csum=colsum_Mat(grad_mat);
        for(int k = 0; k < K; k++){
            grad_vector(t,k)=rsum(k)+csum(k);
        }
    }
    return grad_vector;
}

// [[Rcpp::export]]
mat hess_SEM_undir(vec theta_u, int N, int K, int T_data, float bandwidth, vec data_ids, vec grid_ids, int grid_id_index, vec cluster_ids){
    cube t1(K,K,T_data);
    mat hess_mat(K,K,fill::zeros);
    for(int i = 0; i < (N-1); i++){
        for(int j = i+1; j < N; j++){
            float exp_val=exp(theta_u(cluster_ids(i)-1)+theta_u(cluster_ids(j)-1));
            hess_mat(cluster_ids(i)-1,cluster_ids(j)-1)+=-(exp_val/pow((1+exp_val),2));
        }
    }
    for (int t=0; t<T_data; t++){
        float kernel_val=(1/bandwidth)*epan((data_ids(t)-grid_ids(grid_id_index))/bandwidth);
        for(int k = 0; k < K; k++){
            for(int l = 0; l < K; l++){
                if(k!=l){
                    t1(k,l,t)=(hess_mat(k,l)+hess_mat(l,k))*kernel_val;
                }
            }
        }
        vec rsum=rowsum_Mat(hess_mat);
        vec csum=colsum_Mat(hess_mat);
        for(int k = 0; k < K; k++){
            t1(k,k,t)=(csum(k)+rsum(k)+(2*hess_mat(k,k)))*kernel_val;
        }
    }
    mat t2(K,K,fill::zeros);
    for(int k = 0; k < K; k++){
        for(int l = 0; l < K; l++){
            for (int t=0; t<T_data; t++){
                t2(k,l)+=t1(k,l,t);
            }
        }
    }
    return t2;
}

// [[Rcpp::export]]
float ELBO_conv_SEM_undir(mat gamma, vec pi, vec theta_u, cube network, int N, int K, int T_data, float bandwidth, vec data_ids, vec grid_ids, int grid_id_index, vec cluster_ids){
    vec t1(T_data,fill::zeros);
    for (int t=0; t<T_data; t++){
        for(int i = 0; i < (N-1); i++){
            for(int j = i+1; j < N; j++){
                float exp_val=exp(theta_u(cluster_ids(i)-1)+theta_u(cluster_ids(j)-1));
                t1(t)+=(network(i,j,t)*(theta_u(cluster_ids(i)-1)+theta_u(cluster_ids(j)-1))-log(1+exp_val));
            }
        }
        t1(t)=t1(t)*(1/bandwidth)*epan((data_ids(t)-grid_ids(grid_id_index))/bandwidth);
    }
    vec t2(T_data,fill::zeros);
    for (int t=0; t<T_data; t++){
        for(int i = 0; i < N; i++){
            for(int k = 0; k < K; k++){
                t2(t)+=gamma(i,k)*(log(pi(k))-log(gamma(i,k)));
            }
        }
        t2(t)=t2(t)*(1/bandwidth)*epan((data_ids(t)-grid_ids(grid_id_index))/bandwidth);
    }
    float ELBO_val=sum(t1)+sum(t2);
    return ELBO_val;
}

