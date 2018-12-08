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
//
// [[Rcpp::export]]
arma::mat rcpparma_hello_world() {
    arma::mat m1 = arma::eye<arma::mat>(3, 3);
    arma::mat m2 = arma::eye<arma::mat>(3, 3);
    
    return m1 + 3 * (m1 + m2);
}


// another simple example: outer product of a vector,
// returning a matrix
//
// [[Rcpp::export]]
arma::mat rcpparma_outerproduct(const arma::colvec & x) {
    arma::mat m = x * x.t();
    return m;
}

// and the inner product returns a scalar
//
// [[Rcpp::export]]
double rcpparma_innerproduct(const arma::colvec & x) {
    double v = arma::as_scalar(x.t() * x);
    return v;
}


// and we can use Rcpp::List to return both at the same time
//
// [[Rcpp::export]]
Rcpp::List rcpparma_bothproducts(const arma::colvec & x) {
    arma::mat op = x * x.t();
    double    ip = arma::as_scalar(x.t() * x);
    return Rcpp::List::create(Rcpp::Named("outer")=op,
                              Rcpp::Named("inner")=ip);
}

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
vec grad_EM_dir_K1_oe(vec theta_u, cube network, int N, int K, int T_data){
    vec grad_vector(T_data);
    for (int t=0; t<T_data; t++){
        float grad_mat=0;
        for(int i = 0; i < (N-1); i++){
            for(int j = i+1; j < N; j++){
                float exp_val_1=exp(theta_u(0));
                float exp_val_2=exp(2*theta_u(1));
                int indicator_10=(network(i,j,t)==1)&(network(j,i,t)==0);
                int indicator_01=(network(i,j,t)==0)&(network(j,i,t)==1);
                float grad_matsub=((indicator_10+indicator_01)-((2*exp_val_1)/(1+2*exp_val_1+exp_val_2)));
                grad_mat+=grad_matsub;
            }
        }
        grad_vector(t)=grad_mat;
    }
    return grad_vector;
}

// [[Rcpp::export]]
vec grad_EM_dir_K1_re(vec theta_u, cube network, int N, int K, int T_data){
    vec grad_vector(T_data);
    for (int t=0; t<T_data; t++){
        float grad_mat=0;
        for(int i = 0; i < (N-1); i++){
            for(int j = i+1; j < N; j++){
                float exp_val_1=exp(theta_u(0));
                float exp_val_2=exp(2*theta_u(1));
                int indicator_11=(network(i,j,t)==1)&(network(j,i,t)==1);
                float grad_matsub=((2*indicator_11)-((2*exp_val_2)/(1+2*exp_val_1+exp_val_2)));
                grad_mat+=grad_matsub;
            }
        }
        grad_vector(t)=grad_mat;
    }
    return grad_vector;
}

// [[Rcpp::export]]
float hess_EM_dir_K1_oe(vec theta_u, int N, int K, int T_data, float bandwidth, vec data_ids, vec grid_ids, int grid_id_index){
    vec t2(T_data);
    for (int t=0; t<T_data; t++){
        float hess_mat=0;
        for(int i = 0; i < (N-1); i++){
            for(int j = i+1; j < N; j++){
                float exp_val_1=exp(theta_u(0));
                float exp_val_2=exp(2*theta_u(1));
                float Num=2*(exp_val_1+exp_val_1*exp_val_2);
                float Denom=pow((1+2*exp_val_1+exp_val_2),2);
                hess_mat+=-(Num/Denom);
            }
        }
        float kernel_val=(1/bandwidth)*epan((data_ids(t)-grid_ids(grid_id_index))/bandwidth);
        t2(t)=hess_mat*kernel_val;
    }
    float t1=accu(t2);
    return t1;
}

// [[Rcpp::export]]
float hess_EM_dir_K1_re(vec theta_u, int N, int K, int T_data, float bandwidth, vec data_ids, vec grid_ids, int grid_id_index){
    vec t2(T_data);
    for (int t=0; t<T_data; t++){
        float hess_mat=0;
        for(int i = 0; i < (N-1); i++){
            for(int j = i+1; j < N; j++){
                float exp_val_1=exp(theta_u(0));
                float exp_val_2=exp(2*theta_u(1));
                float Num=4*(exp_val_2+(2*exp_val_1*exp_val_2));
                float Denom=pow((1+2*exp_val_1+exp_val_2),2);
                hess_mat+=-(Num/Denom);
            }
        }
        float kernel_val=(1/bandwidth)*epan((data_ids(t)-grid_ids(grid_id_index))/bandwidth);
        t2(t)=hess_mat*kernel_val;
    }
    float t1=accu(t2);
    return t1;
}

// [[Rcpp::export]]
float hess_EM_dir_K1_oe_re(vec theta_u, int N, int K, int T_data, float bandwidth, vec data_ids, vec grid_ids, int grid_id_index){
    vec t2(T_data);
    for (int t=0; t<T_data; t++){
        float hess_mat=0;
        for(int i = 0; i < (N-1); i++){
            for(int j = i+1; j < N; j++){
                float exp_val_1=exp(theta_u(0));
                float exp_val_2=exp(2*theta_u(1));
                float Num=4*exp_val_1*exp_val_2;
                float Denom=pow((1+2*exp_val_1+exp_val_2),2);
                hess_mat+=(Num/Denom);
            }
        }
        float kernel_val=(1/bandwidth)*epan((data_ids(t)-grid_ids(grid_id_index))/bandwidth);
        t2(t)=hess_mat*kernel_val;
    }
    float t1=accu(t2);
    return t1;
}

// [[Rcpp::export]]
float ELBO_conv_EM_dir_K1(vec theta_u, cube network, int N, int K, int T_data, float bandwidth, vec data_ids, vec grid_ids, int grid_id_index){
    vec t1(T_data,fill::zeros);
    for (int t=0; t<T_data; t++){
        for(int i = 0; i < (N-1); i++){
            for(int j = i+1; j < N; j++){
                float exp_val_1=exp(theta_u(0));
                float exp_val_2=exp(2*theta_u(1));
                int indicator_10=(network(i,j,t)==1)&(network(j,i,t)==0);
                int indicator_01=(network(i,j,t)==0)&(network(j,i,t)==1);
                int indicator_11=(network(i,j,t)==1)&(network(j,i,t)==1);
                t1(t)+=((indicator_10*theta_u(0))+(indicator_01*theta_u(0))+(indicator_11*(2*theta_u(1)))-log(1+2*exp_val_1+exp_val_2));
            }
        }
        t1(t)=t1(t)*(1/bandwidth)*epan((data_ids(t)-grid_ids(grid_id_index))/bandwidth);
    }
    float ELBO_val=sum(t1);
    return ELBO_val;
}

// [[Rcpp::export]]
vec grad_ELBO_K1_CV(float theta_u, vec gamma, cube network, int N, int K, int T_data, int test_node_set_len){
    vec grad_vector(T_data);
    for (int t=0; t<T_data; t++){
        float grad_mat=0;
        for(int i = 0; i < test_node_set_len; i++){
            for(int j = test_node_set_len; j < N; j++){
                float exp_val=exp(theta_u+theta_u);
                float grad_matsub=gamma(i)*gamma(j)*(network(i,j,t)-(exp_val/(1+exp_val)));
                grad_mat+=grad_matsub;
            }
        }
        for(int i = test_node_set_len; i < N; i++){
            if(i!=(N-1)){
                for(int j = i+1; j < N; j++){
                    float exp_val=exp(theta_u+theta_u);
                    float grad_matsub=gamma(i)*gamma(j)*(network(i,j,t)-(exp_val/(1+exp_val)));
                    grad_mat+=grad_matsub;
                }
            } else if(i==(N-1)){
                float exp_val=exp(theta_u+theta_u);
                float grad_matsub=gamma(i)*gamma(N-1)*(network(i,(N-1),t)-(exp_val/(1+exp_val)));
                grad_mat+=grad_matsub;
            }
        }
        grad_vector(t)=2*grad_mat;
    }
    return grad_vector;
}

// [[Rcpp::export]]
float hess_ELBO_K1_CV(float theta_u, vec gamma, cube network, int N, int K, int T_data, float bandwidth, vec data_ids, vec grid_ids, int grid_id_index, int test_node_set_len){
    vec t1(T_data);
    for (int t=0; t<T_data; t++){
        float hess_mat=0;
        for(int i = 0; i < test_node_set_len; i++){
            for(int j = test_node_set_len; j < N; j++){
                float exp_val=exp(theta_u+theta_u);
                float hess_matsub=-(gamma(i)*gamma(j)*(exp_val/(pow((1+exp_val),2))));
                hess_mat+=hess_matsub;
            }
        }
        for(int i = test_node_set_len; i < N; i++){
            if(i!=(N-1)){
                for(int j = i+1; j < N; j++){
                    float exp_val=exp(theta_u+theta_u);
                    float hess_matsub=-(gamma(i)*gamma(j)*(exp_val/(pow((1+exp_val),2))));
                    hess_mat+=hess_matsub;
                }
            } else if(i==(N-1)){
                float exp_val=exp(theta_u+theta_u);
                float hess_matsub=-(gamma(i)*gamma(N-1)*(exp_val/(pow((1+exp_val),2))));
                hess_mat+=hess_matsub;
            }
        }
        float kernel_val=(1/bandwidth)*epan((data_ids(t)-grid_ids(grid_id_index))/bandwidth);
        t1(t)=4*hess_mat*kernel_val;
    }
    float t2=0;
    for (int t=0; t<T_data; t++){
        t2+=t1(t);
    }
    return t2;
}

// [[Rcpp::export]]
float ELBO_conv_K1_CV(vec gamma, float pi, float theta_u, cube network, int N, int K, int T_data, float bandwidth, vec data_ids, vec grid_ids, int grid_id_index, int test_node_set_len){
    vec t1(T_data,fill::zeros);
    for (int t=0; t<T_data; t++){
        for(int i = 0; i < test_node_set_len; i++){
            for(int j = test_node_set_len; j < N; j++){
                float exp_val=exp(theta_u+theta_u);
                t1(t)+=(gamma(i)*gamma(j)*(network(i,j,t)*(theta_u+theta_u)-log(1+exp_val)));
            }
        }
        for(int i = test_node_set_len; i < N; i++){
            if(i!=(N-1)){
                for(int j = i+1; j < N; j++){
                    float exp_val=exp(theta_u+theta_u);
                    t1(t)+=(gamma(i)*gamma(j)*(network(i,j,t)*(theta_u+theta_u)-log(1+exp_val)));
                }
            } else if(i==(N-1)){
                float exp_val=exp(theta_u+theta_u);
                t1(t)+=(gamma(i)*gamma(N-1)*(network(i,(N-1),t)*(theta_u+theta_u)-log(1+exp_val)));
            }
        }
        t1(t)=t1(t)*(1/bandwidth)*epan((data_ids(t)-grid_ids(grid_id_index))/bandwidth);
    }
    vec t2(T_data,fill::zeros);
    for (int t=0; t<T_data; t++){
        for(int i = 0; i < N; i++){
            t2(t)+=gamma(i)*(log(pi)-log(gamma(i)));
        }
        t2(t)=t2(t)*(1/bandwidth)*epan((data_ids(t)-grid_ids(grid_id_index))/bandwidth);
    }
    float ELBO_val=sum(t1)+sum(t2);
    return ELBO_val;
}
