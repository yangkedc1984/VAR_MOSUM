
//#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <iostream>
using namespace Rcpp;
using namespace arma;
using namespace std;
                                    

// [[Rcpp::export(getH_ik_Wald_RCPP)]] //getH_ik_Wald
vec H_ik(mat& x, int i, int k, int p, vec a)
{
  int n = x.n_rows;
  int d = x.n_cols;
  vec ai = a.rows( span((i-1)*(d*p+1)+1 -1 , (i-1)*(d*p+1)+ d*p +1 -1) ); //intercept - index
  
  mat X = x.rows( span(k-p,k-1)) ;
  //arma::mat XX = x.rows( span(k-p,k-1));
  //XX.attr("dim") = Dimension(1, d*p);
  vec V = vectorise(X);
  vec O = ones(1);
  vec VV = join_cols(O,V);
  //V.rows(span(1  ,d*p +1) ) = XX;
  //Numericarma::vector Xx = V.push_front(1); //add intercept
  //vec xr = x.row(k);
  //for(int ii=1; ii< (p+1); i++){
  // xr = x.row(k-ii);
  // V.rows( span((ii-1)*(d+1)+1  , (ii-1)*(d+1)+ d +1 ) ) =  xr.t();
  //} 
  
  double y = x[k,i];
  double  e = y - dot(ai, VV);  //residual
  //vec Vt = transpose(V)
  vec  H_ik = - y*VV + ai * VV.t() * VV- e*VV ;
  
  return H_ik;
  
}


// [[Rcpp::export(makeH_k_Wald_RCPP)]] //makeH_k_Wald
vec H_k(mat& x, int k, int p, vec a)
{
  int d = x.n_cols;
  vec H = zeros(d+  d*d*p); //accounts for intercept 

  for (int ii=1; ii< d+1; ii++) {
  //H.subvec((1-1)*(d*p+1)+1-1, (1-1)*(d*p+1)+ d*p +1-1 )= H_ik(x,1,k,p,a) ;
      H.subvec((ii-1)*(d*p+1)+1-1, (ii-1)*(d*p+1)+ d*p +1-1 )= H_ik(x,ii,k,p,a) ;
      //H = join_cols(H1, H_ik(x,ii,k,p,a));
  };
  return H;
}


// [[Rcpp::export(makeH_l_u_RCPP)]] //makeH_l_u
mat H_l_u(mat& x,  int p, int l,int u, vec a)
{
  int n = x.n_rows;
  int d = x.n_cols;
  int nr = d+ d*d*p;
  int nc = u-l+1;
  mat H; H.zeros(nr,nc); //matrix of H values #accounts for intercept
  for (int t=0; t <(u-l+1); t++ ) {
          H.col(t) = H_k(x, l+t-1, p, a) ;
  };
  return H;
    
}

struct Matrix_types {
  arma::mat m;
  arma::sp_mat M;
}; 
// [[Rcpp::export]] // INTERLEAVE
arma::mat write_rows2(Rcpp::List data,  int nrows, int ncols) {//Rcpp::CharacterVector clss,
  
  const int len = data.length();
  std::vector<Matrix_types> matr(len);
  std::vector<bool> is_dense(len);
  arma::mat result(nrows*len, ncols);
  
  // populate the structs
  for (int j = 0; j < len; j++) {
    //is_dense[j] = (clss[j] == "matrix");
    //if (is_dense[j]) {
      matr[j].m = Rcpp::as<arma::mat>(data[j]);
    //}
    //else {
    //  matr[j].M = Rcpp::as<arma::sp_mat>(data[j]);//SPARSE M
    //}
  }
  
  // populate the result
  for (int i = 0, k = 0; i < nrows; i++) {
    for (int j = 0; j < len; j++, k++) {
      //if (is_dense[j]) {
        result.row(k) = matr[j].m.row(i);
      //}
      //else {
      //  arma::rowvec r(matr[j].M.row(i));
      //  result.row(k) = r;
      //}
    }
  }
  return result;
}

// [[Rcpp::export(get_a_lu_i_RCPP)]] //get_a_lu_i
vec a_lu_i(mat& x, int i, int p, int l, int u)
{
    int d = x.n_cols;
    vec y =   x( span(l-1,u-1),i-1 )  ;
    mat Y = diagmat(y);
    mat O = ones(y.n_elem,1) ;
    mat X;
    rowvec y_soln;  
    if(p==1)  {
       X = join_rows(O, x.rows( l-1-1,u-1-1) ); //regressors ##p = 1 case only
       y_soln = sum(Y * X);
    }  
     if(p==2){
       List L = List::create(x.rows( l-1-1,u-1-1).t(), x.rows( l-2-1,u-2-1).t() ); //construct data list
       X = join_rows(O, write_rows2(L,d,y.n_elem ).t() ) ;
       y_soln = sum(Y * X);
     }  
    mat X_soln = (X.t() * X ).i();
    vec a_out = (y_soln * X_soln).t(); 
    return a_out;
}

// [[Rcpp::export(make_a_lu_RCPP)]] //make_a_lu
vec make_a_lu(mat& x, int p, int l, int u)
{
  int d = x.n_cols;
  vec a_lu = zeros(d + d*d*p);
  for (int ii=1; ii< d+1; ii++) {
    a_lu.subvec( (ii-1)*(d*p+1)+1-1, (ii-1)*(d*p+1)+ d*p +1-1 ) = a_lu_i(x,ii,p,l,u)  ;
  };
  return a_lu;
}


// [[Rcpp::export]]
arma::sp_mat blockDiag( arma::field<mat>& xlist ) {
  
  //xlist: list of matrices 
  
  unsigned int n = xlist.n_rows ;
  int dimen = 0 ;
  arma::ivec dimvec(n) ;
  
  for(unsigned int i=0; i<n; i++) {
    dimvec(i) = xlist(i,0).n_rows ;
    dimen += dimvec(i) ;
  }

  sp_mat X(dimen,dimen);
  int idx=0;

  for(unsigned int i=0; i<n; i++) {
    X.submat( idx, idx, idx + dimvec(i) - 1, idx + dimvec(i) - 1 ) = xlist(i,0) ;
    idx = idx + dimvec(i) ;
  }
  return(X);
}


// [[Rcpp::export(get_V_nk_RCPP)]] //get_V_nk
sp_mat V_nk(mat x, int p, int l, int u)
 {
  int d = x.n_cols;
  mat xk;
  mat O = ones(u-l+1,1) ;
  if(p==1) xk = join_rows(O, x.rows( l-1-1,u-1-1) );
  if(p==2) xk = join_rows(O, x.rows( l-1-1,u-1-1) , x.rows( l-2-1,u-2-1));
  mat C =  xk.t() * xk /(u-l);
  field<mat> Vlist(d);
  for (int ii=0; ii< d; ii++) {
    Vlist(ii) = C; //coerce into block diagonal form
  };
  sp_mat out = blockDiag(Vlist);
  return out;
 }

// [[Rcpp::export(getsigma_i_kLOCAL1_RCPP)]] //getsigma_i_kLOCAL1
double sigma_i_k(mat x, int i,int k,int G,int p,vec a_upper,vec a_lower)
{
  mat x_upper; mat x_lower;
  mat O = ones(G,1) ;
  if(p==1) x_upper = join_rows(O, x.rows( k-1,k+G-1-1) ); //upper sample
  if(p==2) x_upper = join_rows(O, x.rows( k-1,k+G-1-1), x.rows(k-1-1, k+G-2-1 ) ); 
  rowvec res_upper =  x( span(k+1-1,k+G-1), i).t() - a_upper.t() * x_upper.t(); //upper residuals
    
  if(p==1)x_lower = join_rows(O, x.rows( k-G-1,k-1-1) );
  if(p==2)x_lower  = join_rows(O, x.rows( k-G-1,k-1-1), x.rows(k-G-2, k-2-1));
                    
  rowvec res_lower =  x( span(k-G+1-1,k-1), i).t() - a_lower.t() * x_lower.t(); //lower residuals
  double sigma_i =  (sum(square(res_upper)) + sum(square(res_lower)) ) /(2*G);
  return sigma_i;
  
}

// [[Rcpp::export(getsigma_d_kLOCAL1_RCPP)]] //getsigma_d_kLOCAL1
vec sigma_d_k(mat x, int k,int G,int p,vec a_upper,vec a_lower)
{
  int d = x.n_cols;
  vec sigma_d = zeros(d);
  // for (int ii=1; ii< d+1; ii++) {
  //     sigma_d(ii-1) = sigma_i_k(x,ii,k,G,p,a_upper.subvec( (ii-1)*(d*p+1)+1-1, (ii-1)*(d*p+1)+ d*p +1-1),a_lower.subvec( (ii-1)*(d*p+1)+1-1, (ii-1)*(d*p+1)+ d*p +1-1) ) ;
  // };
    return sigma_d;
}