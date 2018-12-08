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
cube gamma_update_dir(mat gamma, vec pi, cube theta, cube network, int N, int K, int T_grid, vec grid_ids){
    cube quad_lin_coeff(N,K,2);
    for(int i = 0; i < N; i++){
        if(i!=(N-1)){
            for(int k = 0; k < K; k++){
                float t1=0;
                for (int t=0; t<T_grid; t++){
                    for(int j = i+1; j < N; j++){
                        for(int l = 0; l < K; l++){
                            float exp_val_1=exp(theta(t,k,0));
                            float exp_val_2=exp(theta(t,l,0));
                            float exp_val_3=exp(theta(t,k,1)+theta(t,l,1));
                            vec alpha(4,fill::zeros);
                            alpha(1)=exp_val_1;
                            alpha(2)=exp_val_2;
                            alpha(3)=exp_val_3;
                            float alpha_max=alpha.max();
                            float exp_val_1_mod=exp(theta(t,k,0)-alpha_max);
                            float exp_val_2_mod=exp(theta(t,l,0)-alpha_max);
                            float exp_val_3_mod=exp(theta(t,k,1)+theta(t,l,1)-alpha_max);
                            float log_exp_val=alpha_max+log(exp(-alpha_max)+exp_val_1_mod+exp_val_2_mod+exp_val_3_mod);
                            int indicator_10=(network(i,j,(grid_ids(t)-1))==1)&(network(j,i,(grid_ids(t)-1))==0);
                            int indicator_01=(network(i,j,(grid_ids(t)-1))==0)&(network(j,i,(grid_ids(t)-1))==1);
                            int indicator_11=(network(i,j,(grid_ids(t)-1))==1)&(network(j,i,(grid_ids(t)-1))==1);
                            t1+=((gamma(j,l)/(2*gamma(i,k)))*((indicator_10*theta(t,k,0))+(indicator_01*theta(t,l,0))+((indicator_11)*(theta(t,k,1)+theta(t,l,1)))-log_exp_val));
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
mat grad_SEM_dir_oe(mat theta_u, cube network, int N, int K, vec nonzero_ids, int nonzero_ids_len, vec cluster_ids){
    mat grad_vector(nonzero_ids_len,K);
    for (int t=0; t<nonzero_ids_len; t++){
        vec grad_vec_t(K,fill::zeros);
        for(int i = 0; i < (N-1); i++){
            for(int j = i+1; j < N; j++){
                float exp_val_1=exp(theta_u(cluster_ids(i)-1,0));
                float exp_val_2=exp(theta_u(cluster_ids(j)-1,0));
                float exp_val_3=exp(theta_u(cluster_ids(i)-1,1)+theta_u(cluster_ids(j)-1,1));
                int indicator_10=(network(i,j,(nonzero_ids(t)-1))==1)&(network(j,i,(nonzero_ids(t)-1))==0);
                int indicator_01=(network(i,j,(nonzero_ids(t)-1))==0)&(network(j,i,(nonzero_ids(t)-1))==1);
                float grad_val_1=(indicator_10-(exp_val_1/(1+exp_val_1+exp_val_2+exp_val_3)));
                float grad_val_2=(indicator_01-(exp_val_2/(1+exp_val_1+exp_val_2+exp_val_3)));
                grad_vec_t(cluster_ids(i)-1)+=grad_val_1;
                grad_vec_t(cluster_ids(j)-1)+=grad_val_2;
            }
        }
        for (int k = 0; k < K; k++){
            grad_vector(t,k)=grad_vec_t(k);
        }
    }
    return grad_vector;
}

// [[Rcpp::export]]
mat grad_SEM_dir_re(mat theta_u, cube network, int N, int K, vec nonzero_ids, int nonzero_ids_len, vec cluster_ids){
    mat grad_vector(nonzero_ids_len,K);
    for (int t=0; t<nonzero_ids_len; t++){
        vec grad_vec_t(K,fill::zeros);
        for(int i = 0; i < (N-1); i++){
            for(int j = i+1; j < N; j++){
                float exp_val_1=exp(theta_u(cluster_ids(i)-1,0));
                float exp_val_2=exp(theta_u(cluster_ids(j)-1,0));
                float exp_val_3=exp(theta_u(cluster_ids(i)-1,1)+theta_u(cluster_ids(j)-1,1));
                int indicator_11=(network(i,j,(nonzero_ids(t)-1))==1)&(network(j,i,(nonzero_ids(t)-1))==1);
                float grad_val=(indicator_11-(exp_val_3/(1+exp_val_1+exp_val_2+exp_val_3)));
                grad_vec_t(cluster_ids(i)-1)+=grad_val;
                grad_vec_t(cluster_ids(j)-1)+=grad_val;
            }
        }
        for (int k = 0; k < K; k++){
            grad_vector(t,k)=grad_vec_t(k);
        }
    }
    return grad_vector;
}

// [[Rcpp::export]]
mat hess_SEM_dir_oe(mat theta_u, int N, int K, int T_data, float bandwidth, vec data_ids, vec grid_ids, int grid_id_index, vec cluster_ids){
    cube t1(K,K,T_data);
    mat hess_mat_1(K,K,fill::zeros);
    mat hess_mat_2(K,K,fill::zeros);
    mat hess_mat_3(K,K,fill::zeros);
    for(int i = 0; i < (N-1); i++){
        for(int j = i+1; j < N; j++){
            float exp_val_1=exp(theta_u(cluster_ids(i)-1,0));
            float exp_val_2=exp(theta_u(cluster_ids(j)-1,0));
            float exp_val_3=exp(theta_u(cluster_ids(i)-1,1)+theta_u(cluster_ids(j)-1,1));
            hess_mat_1((cluster_ids(i)-1),(cluster_ids(j)-1))+=((exp_val_1*exp_val_2)/(pow((1+exp_val_1+exp_val_2+exp_val_3),2)));
            hess_mat_2((cluster_ids(i)-1),(cluster_ids(j)-1))+=((exp_val_1+exp_val_1*exp_val_3)/(pow((1+exp_val_1+exp_val_2+exp_val_3),2)));
            hess_mat_3((cluster_ids(i)-1),(cluster_ids(j)-1))+=((exp_val_2+exp_val_2*exp_val_3)/(pow((1+exp_val_1+exp_val_2+exp_val_3),2)));
        }
    }
    
    for (int t=0; t<T_data; t++){
        float kernel_val=(1/bandwidth)*epan((data_ids(t)-grid_ids(grid_id_index))/bandwidth);
        for(int k = 0; k < K; k++){
            for(int l = 0; l < K; l++){
                if(k!=l){
                    t1(k,l,t)=(hess_mat_1(k,l)+hess_mat_1(l,k))*kernel_val;
                }
            }
        }
        vec rsum_1=rowsum_Mat(hess_mat_1);
        vec csum_1=colsum_Mat(hess_mat_1);
        vec rsum_2=rowsum_Mat(hess_mat_2);
        vec csum_3=colsum_Mat(hess_mat_3);
        for(int k = 0; k < K; k++){
            t1(k,k,t)=(-csum_3(k)-rsum_2(k)-(rsum_1(k)+csum_1(k)-2*hess_mat_1(k,k)))*kernel_val;
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
mat hess_SEM_dir_re(mat theta_u, int N, int K, int T_data, float bandwidth, vec data_ids, vec grid_ids, int grid_id_index, vec cluster_ids){
    cube t1(K,K,T_data);
    mat hess_mat(K,K,fill::zeros);
    for(int i = 0; i < (N-1); i++){
        for(int j = i+1; j < N; j++){
            float exp_val_1=exp(theta_u(cluster_ids(i)-1,0));
            float exp_val_2=exp(theta_u(cluster_ids(j)-1,0));
            float exp_val_3=exp(theta_u(cluster_ids(i)-1,1)+theta_u(cluster_ids(j)-1,1));
            hess_mat((cluster_ids(i)-1),(cluster_ids(j)-1))+=(((1+exp_val_1+exp_val_2)*exp_val_3)/(pow((1+exp_val_1+exp_val_2+exp_val_3),2)));
        }
    }
    
    for (int t=0; t<T_data; t++){
        float kernel_val=(1/bandwidth)*epan((data_ids(t)-grid_ids(grid_id_index))/bandwidth);
        for(int k = 0; k < K; k++){
            for(int l = 0; l < K; l++){
                if(k!=l){
                    t1(k,l,t)=-(hess_mat(k,l)+hess_mat(l,k))*kernel_val;
                }
            }
        }
        vec rsum=rowsum_Mat(hess_mat);
        vec csum=colsum_Mat(hess_mat);
        for(int k = 0; k < K; k++){
            t1(k,k,t)=(-csum(k)-rsum(k)-2*hess_mat(k,k))*kernel_val;
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
mat hess_SEM_dir_oe_re(mat theta_u, int N, int K, int T_data, float bandwidth, vec data_ids, vec grid_ids, int grid_id_index, vec cluster_ids){
    cube t1(K,K,T_data);
    mat hess_mat_1(K,K,fill::zeros);
    mat hess_mat_2(K,K,fill::zeros);
    for(int i = 0; i < (N-1); i++){
        for(int j = i+1; j < N; j++){
            float exp_val_1=exp(theta_u(cluster_ids(i)-1,0));
            float exp_val_2=exp(theta_u(cluster_ids(j)-1,0));
            float exp_val_3=exp(theta_u(cluster_ids(i)-1,1)+theta_u(cluster_ids(j)-1,1));
            hess_mat_1((cluster_ids(i)-1),(cluster_ids(j)-1))+=((exp_val_1*exp_val_3)/(pow((1+exp_val_1+exp_val_2+exp_val_3),2)));
            hess_mat_2((cluster_ids(i)-1),(cluster_ids(j)-1))+=((exp_val_2*exp_val_3)/(pow((1+exp_val_1+exp_val_2+exp_val_3),2)));
        }
    }
    
    for (int t=0; t<T_data; t++){
        float kernel_val=(1/bandwidth)*epan((data_ids(t)-grid_ids(grid_id_index))/bandwidth);
        for(int k = 0; k < K; k++){
            for(int l = 0; l < K; l++){
                if(k!=l){
                    t1(k,l,t)=(hess_mat_1(k,l)+hess_mat_2(l,k))*kernel_val;
                }
            }
        }
        vec rsum_1=rowsum_Mat(hess_mat_1);
        vec csum_2=colsum_Mat(hess_mat_2);
        for(int k = 0; k < K; k++){
            t1(k,k,t)=(csum_2(k)+rsum_1(k)+(hess_mat_1(k,k)+hess_mat_2(k,k)))*kernel_val;
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
float ELBO_conv_SEM_dir(mat gamma, vec pi, mat theta_u, cube network, int N, int K, int T_data, float bandwidth, vec data_ids, vec grid_ids, int grid_id_index, vec cluster_ids){
    vec t1(T_data,fill::zeros);
    for (int t=0; t<T_data; t++){
        for(int i = 0; i < (N-1); i++){
            for(int j = i+1; j < N; j++){
                float exp_val_1=exp(theta_u(cluster_ids(i)-1,0));
                float exp_val_2=exp(theta_u(cluster_ids(j)-1,0));
                float exp_val_3=exp(theta_u(cluster_ids(i)-1,1)+theta_u(cluster_ids(j)-1,1));
                int indicator_10=(network(i,j,t)==1)&(network(j,i,t)==0);
                int indicator_01=(network(i,j,t)==0)&(network(j,i,t)==1);
                int indicator_11=(network(i,j,t)==1)&(network(j,i,t)==1);
                t1(t)+=((indicator_10*theta_u(cluster_ids(i)-1,0))+(indicator_01*theta_u(cluster_ids(j)-1,0))+(indicator_11*(theta_u(cluster_ids(i)-1,1)+theta_u(cluster_ids(j)-1,1)))-log(1+exp_val_1+exp_val_2+exp_val_3));
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
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Model selection

// [[Rcpp::export]]
mat VK_dir_oe(cube theta, cube network, int N, int K, int T_data, vec cluster_ids){
    mat grad_vector(T_data,K);
    for (int t=0; t<T_data; t++){
        vec grad_vec_t(K,fill::zeros);
        for(int i = 0; i < (N-1); i++){
            for(int j = i+1; j < N; j++){
                float exp_val_1=exp(theta(t,cluster_ids(i)-1,0));
                float exp_val_2=exp(theta(t,cluster_ids(j)-1,0));
                float exp_val_3=exp(theta(t,cluster_ids(i)-1,1)+theta(t,cluster_ids(j)-1,1));
                int indicator_10=(network(i,j,t)==1)&(network(j,i,t)==0);
                int indicator_01=(network(i,j,t)==0)&(network(j,i,t)==1);
                float grad_val_1=(indicator_10-(exp_val_1/(1+exp_val_1+exp_val_2+exp_val_3)));
                float grad_val_2=(indicator_01-(exp_val_2/(1+exp_val_1+exp_val_2+exp_val_3)));
                grad_vec_t(cluster_ids(i)-1)+=grad_val_1;
                grad_vec_t(cluster_ids(j)-1)+=grad_val_2;
            }
        }
        for (int k = 0; k < K; k++){
            grad_vector(t,k)=grad_vec_t(k);
        }
    }
    return grad_vector;
}

// [[Rcpp::export]]
mat VK_dir_re(cube theta, cube network, int N, int K, int T_data, vec cluster_ids){
    mat grad_vector(T_data,K);
    for (int t=0; t<T_data; t++){
        vec grad_vec_t(K,fill::zeros);
        for(int i = 0; i < (N-1); i++){
            for(int j = i+1; j < N; j++){
                float exp_val_1=exp(theta(t,cluster_ids(i)-1,0));
                float exp_val_2=exp(theta(t,cluster_ids(j)-1,0));
                float exp_val_3=exp(theta(t,cluster_ids(i)-1,1)+theta(t,cluster_ids(j)-1,1));
                int indicator_11=(network(i,j,t)==1)&(network(j,i,t)==1);
                float grad_val=(indicator_11-(exp_val_3/(1+exp_val_1+exp_val_2+exp_val_3)));
                grad_vec_t(cluster_ids(i)-1)+=grad_val;
                grad_vec_t(cluster_ids(j)-1)+=grad_val;
            }
        }
        for (int k = 0; k < K; k++){
            grad_vector(t,k)=grad_vec_t(k);
        }
    }
    return grad_vector;
}

// [[Rcpp::export]]
mat HK_dir_oe(cube theta, int N, int K, int T_data, vec cluster_ids){
    cube t1(K,K,T_data);
    mat hess_mat_1(K,K,fill::zeros);
    mat hess_mat_2(K,K,fill::zeros);
    mat hess_mat_3(K,K,fill::zeros);
    for (int t=0; t<T_data; t++){
        for(int i = 0; i < (N-1); i++){
            for(int j = i+1; j < N; j++){
                float exp_val_1=exp(theta(t,cluster_ids(i)-1,0));
                float exp_val_2=exp(theta(t,cluster_ids(j)-1,0));
                float exp_val_3=exp(theta(t,cluster_ids(i)-1,1)+theta(t,cluster_ids(j)-1,1));
                hess_mat_1((cluster_ids(i)-1),(cluster_ids(j)-1))+=((exp_val_1*exp_val_2)/(pow((1+exp_val_1+exp_val_2+exp_val_3),2)));
                hess_mat_2((cluster_ids(i)-1),(cluster_ids(j)-1))+=((exp_val_1+exp_val_1*exp_val_3)/(pow((1+exp_val_1+exp_val_2+exp_val_3),2)));
                hess_mat_3((cluster_ids(i)-1),(cluster_ids(j)-1))+=((exp_val_2+exp_val_2*exp_val_3)/(pow((1+exp_val_1+exp_val_2+exp_val_3),2)));
            }
        }
    }
    
    for (int t=0; t<T_data; t++){
        for(int k = 0; k < K; k++){
            for(int l = 0; l < K; l++){
                if(k!=l){
                    t1(k,l,t)=(hess_mat_1(k,l)+hess_mat_1(l,k));
                }
            }
        }
        vec rsum_1=rowsum_Mat(hess_mat_1);
        vec csum_1=colsum_Mat(hess_mat_1);
        vec rsum_2=rowsum_Mat(hess_mat_2);
        vec csum_3=colsum_Mat(hess_mat_3);
        for(int k = 0; k < K; k++){
            t1(k,k,t)=(-csum_3(k)-rsum_2(k)-(rsum_1(k)+csum_1(k)-2*hess_mat_1(k,k)));
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
    return -t2;
}

// [[Rcpp::export]]
mat HK_dir_re(cube theta, int N, int K, int T_data, vec cluster_ids){
    cube t1(K,K,T_data);
    mat hess_mat(K,K,fill::zeros);
    for (int t=0; t<T_data; t++){
        for(int i = 0; i < (N-1); i++){
            for(int j = i+1; j < N; j++){
                float exp_val_1=exp(theta(t,cluster_ids(i)-1,0));
                float exp_val_2=exp(theta(t,cluster_ids(j)-1,0));
                float exp_val_3=exp(theta(t,cluster_ids(i)-1,1)+theta(t,cluster_ids(j)-1,1));
                hess_mat((cluster_ids(i)-1),(cluster_ids(j)-1))+=(((1+exp_val_1+exp_val_2)*exp_val_3)/(pow((1+exp_val_1+exp_val_2+exp_val_3),2)));
            }
        }
    }
    
    for (int t=0; t<T_data; t++){
        for(int k = 0; k < K; k++){
            for(int l = 0; l < K; l++){
                if(k!=l){
                    t1(k,l,t)=-(hess_mat(k,l)+hess_mat(l,k));
                }
            }
        }
        vec rsum=rowsum_Mat(hess_mat);
        vec csum=colsum_Mat(hess_mat);
        for(int k = 0; k < K; k++){
            t1(k,k,t)=(-csum(k)-rsum(k)-2*hess_mat(k,k));
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
    return -t2;
}

// [[Rcpp::export]]
mat HK_dir_oe_re(cube theta, int N, int K, int T_data, vec cluster_ids){
    cube t1(K,K,T_data);
    mat hess_mat_1(K,K,fill::zeros);
    mat hess_mat_2(K,K,fill::zeros);
    for (int t=0; t<T_data; t++){
        for(int i = 0; i < (N-1); i++){
            for(int j = i+1; j < N; j++){
                float exp_val_1=exp(theta(t,cluster_ids(i)-1,0));
                float exp_val_2=exp(theta(t,cluster_ids(j)-1,0));
                float exp_val_3=exp(theta(t,cluster_ids(i)-1,1)+theta(t,cluster_ids(j)-1,1));
                hess_mat_1((cluster_ids(i)-1),(cluster_ids(j)-1))+=((exp_val_1*exp_val_3)/(pow((1+exp_val_1+exp_val_2+exp_val_3),2)));
                hess_mat_2((cluster_ids(i)-1),(cluster_ids(j)-1))+=((exp_val_2*exp_val_3)/(pow((1+exp_val_1+exp_val_2+exp_val_3),2)));
            }
        }
    }
    
    for (int t=0; t<T_data; t++){
        for(int k = 0; k < K; k++){
            for(int l = 0; l < K; l++){
                if(k!=l){
                    t1(k,l,t)=(hess_mat_1(k,l)+hess_mat_2(l,k));
                }
            }
        }
        vec rsum_1=rowsum_Mat(hess_mat_1);
        vec csum_2=colsum_Mat(hess_mat_2);
        for(int k = 0; k < K; k++){
            t1(k,k,t)=(csum_2(k)+rsum_1(k)+(hess_mat_1(k,k)+hess_mat_2(k,k)));
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
    return -t2;
}



