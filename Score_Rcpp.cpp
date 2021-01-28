
//#include <Rcpp.h>
// #include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// // [[Rcpp::depends(RcppParallel)]]
// // [[Rcpp::plugins(openmp)]]

// #include <iostream>
//#include <RcppParallel.h>
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;
using namespace arma;
using namespace std;


// [[Rcpp::export(getH_ik_RCPP)]] //getH_ik
vec H_ik(mat& x, int i, int k, int p, mat Phi, mat eps)
{
  int n = x.n_rows;
  int d = x.n_cols;
  vec ai = Phi.row(i-1).t(); //intercept - index
  
  mat X = x.rows( span(k-p-1,k-1-1)) ;

  vec V = vectorise(flipud(X),1).t();
  vec O = ones(1);
  vec VV = join_cols(O,V);
  
  double y = x(k-1,i-1); //double aV = dot(ai, VV);
  double  e = eps(k-1,i-1);  // y - aV;  //residual
  
  vec  out = - y*VV +  VV*VV.t()*ai - e*VV ;

  return out;
  
}

// [[Rcpp::export(getH_ik_univ)]] //getH_ik
vec H_ik_univ(mat& x, int i, int k, int p, mat Phi, mat eps)
{
  int n = x.n_rows;
  int d = x.n_cols;
  vec ai = Phi.t(); //intercept - index

  mat X = x( span(k-p-1,k-1-1), i-1) ; //depends on self only
  
  vec V = vectorise(flipud(X),1).t(); //time order
  vec O = ones(1); //intercept
  vec VV = join_cols(O,V);
  
  double y = x(k-1,i-1); //double aV = dot(ai, VV);
  double  e = eps(k-1,i-1);  // y - aV;  //residual
  
  vec  out = - y*VV +  VV*VV.t()*ai - e*VV ;
  
  return out;
  
}





// [[Rcpp::export(makeH_k_RCPP)]] //makeH_k
vec H_k(mat& x, int k, int p, mat Phi, mat eps)
{
  int d = x.n_cols;
  vec H = zeros(d+  d*d*p); //accounts for intercept

  for (int ii=1; ii< d+1; ii++) {
    H.subvec((ii-1)*(d*p+1), (ii-1)*(d*p+1)+ d*p  )= H_ik(x,ii,k,p,Phi,eps) ;
  };
  return H;
}


// [[Rcpp::export(makeH_k_univ)]] //makeH_k
vec H_k_univ(mat& x, int k, int p,  arma::field<mat> PhiList, mat eps)
{
  int d = x.n_cols;
  vec H = zeros(d+  d*p); //accounts for intercept
  
  for (int ii=1; ii< d+1; ii++) {
    H.subvec((ii-1)*(p+1), (ii-1)*(p+1)+ (p+1)-1  )= H_ik_univ(x,ii,k,p,PhiList(ii-1),eps) ;
  };
  return H;
}


// [[Rcpp::export(makeH_all_RCPP)]] //makeH_all
mat H_all(mat& x,  int p, int G, mat Phi, mat eps)
{
  int n = x.n_rows;
  int d = x.n_cols;
  int nr = d+ d*d*p;
  //int nc = u-l+1;
  mat H; H.zeros(nr,n); //matrix of H values #accounts for intercept
  for (int t=p+1; t <(n-p); t++ ) {
    H.col(t) = H_k(x, t, p, Phi, eps) ;//-1-1
  };
  return H;

}


// [[Rcpp::export(makeH_all_univ)]] //makeH_all
mat H_all_univ(mat& x,  int p, int G, arma::field<mat> PhiList, mat eps)
{
  int n = x.n_rows;
  int d = x.n_cols;
  int nr = d+ d*p;
  //int nc = u-l+1;
  mat H; H.zeros(nr,n); mat xx;  //matrix of H values #accounts for intercept
  for (int t=p+1; t <(n-p); t++ ) {
    mat xin =x;
    H.col(t) = H_k_univ(xin, t, p, PhiList, eps) ;//-1-1
  };
  return H;
}


// struct Matrix_types {
//   arma::mat m;
//   arma::sp_mat M;
// }; 
// // [[Rcpp::export]] // INTERLEAVE rows of input list
// arma::mat write_rows2(Rcpp::List data,  int nrows, int ncols) {//Rcpp::CharacterVector clss,
//   
//   const int len = data.length();
//   std::vector<Matrix_types> matr(len);
//   std::vector<bool> is_dense(len);
//   arma::mat result(nrows*len, ncols);
//   
//   // populate the structs
//   for (int j = 0; j < len; j++) {
//     //is_dense[j] = (clss[j] == "matrix");
//     //if (is_dense[j]) {
//     matr[j].m = Rcpp::as<arma::mat>(data[j]);
//     //}
//     //else {
//     //  matr[j].M = Rcpp::as<arma::sp_mat>(data[j]);//SPARSE M
//     //}
//   }
//   
//   // populate the result
//   for (int i = 0, k = 0; i < nrows; i++) {
//     for (int j = 0; j < len; j++, k++) {
//       //if (is_dense[j]) {
//       result.row(k) = matr[j].m.row(i);
//       //}
//       //else {
//       //  arma::rowvec r(matr[j].M.row(i));
//       //  result.row(k) = r;
//       //}
//     }
//   }
//   return result;
// }


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


// // [[Rcpp::export(getsigma_i_kLOCAL1_RCPP)]] //getsigma_i_kLOCAL1
// double sigma_i_k(mat x, int i,int k,int G,int p,vec a_upper,vec a_lower)
// {
//   mat x_upper; mat x_lower;
//   mat O = ones(G,1) ;
//   if(p==1) x_upper = join_rows(O, x.rows( k-1,k+G-1-1) ); //upper sample
//   if(p==2) x_upper = join_rows(O, x.rows( k-1,k+G-1-1), x.rows(k-1-1, k+G-2-1 ) ); 
//   if(p==3) x_upper = join_rows(O, x.rows( k-1,k+G-1-1), x.rows(k-1-1, k+G-2-1 ), x.rows(k-2-1, k+G-3-1 ) ); 
//   if(p==4){
//     x_upper = join_rows(x.rows( k-1,k+G-1-1), x.rows(k-1-1, k+G-2-1 ), x.rows(k-2-1, k+G-3-1 ),x.rows(k-3-1, k+G-4-1 ) );
//     x_upper.insert_cols(0,1); 
//   } ;
//   rowvec res_upper =  x( span(k+1-1,k+G-1), i-1).t() - a_upper.t() * x_upper.t(); //upper residuals
//   
//   if(p==1)x_lower = join_rows(O, x.rows( k-G-1,k-1-1) );
//   if(p==2)x_lower  = join_rows(O, x.rows( k-G-1,k-1-1), x.rows(k-G-2, k-2-1));
//   if(p==3)x_lower  = join_rows(O, x.rows( k-G-1,k-1-1), x.rows(k-G-2, k-2-1), x.rows(k-G-3, k-3-1));
//   if(p==4){
//     x_lower = join_rows(x.rows( k-G-1,k-1-1), x.rows(k-G-2, k-2-1), x.rows(k-G-3, k-3-1),x.rows(k-G-4, k-4-1) );
//     x_lower.insert_cols(0,1);
//   };
//   rowvec res_lower =  x( span(k-G+1-1,k-1), i-1).t() - a_lower.t() * x_lower.t(); //lower residuals
//   double sigma_i =  (sum(square(res_upper)) + sum(square(res_lower)) ) /(2*G);
//   return sigma_i;
//   
// }

// // [[Rcpp::export(getsigma_d_kLOCAL1_RCPP)]] //getsigma_d_kLOCAL1
// vec sigma_d_k(mat x, int k,int G,int p,vec a_upper,vec a_lower)
// {
//   int d = x.n_cols;
//   vec sigma_d = zeros(d);
//   for (int ii=1; ii< d+1; ii++) {
//     sigma_d(ii-1) = sigma_i_k(x,ii,k,G,p, a_upper( span((ii-1)*(d*p+1)+1-1, (ii-1)*(d*p+1)+ d*p +1-1)),
//             a_lower( span((ii-1)*(d*p+1)+1-1, (ii-1)*(d*p+1)+ d*p +1-1)) ) ;
//   };
//   return sigma_d;
// }

 
// // SIGMA ESTIMATORS ------------------------

// [[Rcpp::export(get_DiagH_RCPP)]] //get_DiagH
sp_mat DiagH(mat x, int k, int G,int p, mat h_all)//sp_mat
{
  int n = x.n_rows;
  int d = x.n_cols;
  field<mat> H_list(d);
  mat H_u_i; mat H_l_i; mat H_out;
  vec Hbar_u; vec Hbar_l;
  vec evals; mat evecs;
  for (int ii=1; ii< d+1; ii++) {
    H_u_i = h_all.rows( ((ii-1)*(d*p+1)+1-1), ((ii-1)*(d*p+1)+ d*p +1-1) ).cols((k-1),(k+G-1-1) );
    H_l_i = h_all.rows( ((ii-1)*(d*p+1)+1-1), ((ii-1)*(d*p+1)+ d*p +1-1) ).cols((k-G-1),(k-1-1) );
    Hbar_u = mean(H_u_i,1) ;
    H_u_i.each_col() -= Hbar_u ;//H_u_centred
    Hbar_l = mean(H_l_i,1) ;
    H_l_i.each_col() -= Hbar_l; //H_l_centred
    H_out = H_l_i * H_l_i.t() + H_u_i * H_u_i.t();
    //eigen decomposition
    eig_sym(evals,evecs,H_out);
    H_out = evecs * diagmat(pow(evals,-0.5))*  evecs.t();
    H_list(ii-1) = H_out;
  }
  sp_mat Sig_ =  sqrt(2*G) * blockDiag(H_list); //coerce into block diagonal form
  return Sig_;
}

// 

// [[Rcpp::export(get_FullH_RCPP)]] //get_FullH
mat FullH(mat x, int k, int G, mat h_all)
{
  int n = x.n_rows;
  int d = x.n_cols;
  //field<mat> H_list(d);
  mat H_u = h_all.cols((k-1),(k+G-1-1) );
  mat H_l = h_all.cols((k-G-1),(k-1-1) );
  mat H_out;
  vec Hbar_u; vec Hbar_l;
  vec evals; mat evecs;

  Hbar_u = mean(H_u,1) ;
  H_u.each_col() -= Hbar_u ;//H_u_centred
  Hbar_l = mean(H_l,1) ;
  H_l.each_col() -= Hbar_l; //H_l_centred
  H_out = H_l * H_l.t() + H_u * H_u.t();
  //eigen decomposition
  eig_sym(evals,evecs,H_out);
  H_out = evecs * diagmat(pow(evals/(2*G),-0.5))*  evecs.t();


  //DEAL WITH OVERFLOW
  // int n = x.n_rows;
  // int d = x.n_cols;
  // //field<mat> H_list(d);
  // mat H_; Mat<long double> H_out;
  // vec Hbar_u; vec Hbar_l;
  // vec evals; mat evecs; std::vector<long double> evld;
  //
  // Hbar_u = mean(H_u,1) ;
  // H_u.each_col() -= Hbar_u ;//H_u_centred
  // Hbar_l = mean(H_l,1) ;
  // H_l.each_col() -= Hbar_l; //H_l_centred
  // H_ = H_l * H_l.t() + H_u * H_u.t();
  // //eigen decomposition
  // eig_sym(evals,evecs,H_);
  // evld = conv_to<std::vector<long double>>::from(evals);
  // //Col<long double> pp = conv_to<Col<long double>>::from(powl(evld/(2*G),-0.5) );
  // Mat<long double> D = diagmat()evld;
  // H_out = evecs * (D)*  evecs.t();
  return H_out; // RETURNS NaNs
}

// [[Rcpp::export(getA_RCPP)]] //getA
vec getA(mat x, int k, int G,int p, mat eps, mat h_all)
{
  int d = x.n_cols;
  mat r = h_all.cols( (k+1-1),(k+G-1) ); 
  mat l = h_all.cols( (k-G+1-1), (k-1) ); //left/right of window
   
  vec A = sum(r,1) - sum(l,1) ;
  return A;
}

// sigma_i options---------------------------------------------

// [[Rcpp::export(getsigma_iGlobal_RCPP)]] //getsigma_iGlobal
double getsigma_iGlobal(mat eps, int p, int i){
  int n = eps.n_rows;
  vec out = square(eps.col(i-1).rows(p,n-1)) ;
  return mean(out);//mean(out);
}
// [[Rcpp::export(getsigma_dGlobal_RCPP)]] //getsigma_dGlobal
mat getsigma_dGlobal(mat eps, int p){
  // int n = eps.n_rows;
  // int d = eps.n_cols;
  // vec out(d);
  // for(int i = 1; i<d+1; i++){
  //   out(i-1) = getsigma_iGlobal(eps,p,i);
  // }
  // return conv_to<mat>::from(out);//mean(out);
  return cov(eps);
}


// // [[Rcpp::export(getsigma_iLocal_RCPP)]] //getsigma_iLocal
// double getsigma_iLocal(mat eps, int i, int k, int G){
//   double el = var( eps.rows( (k-G+1-1),(k-1)).col(i-1)) ; 
//   double eu = var( eps.rows( (k+1-1),(k+G-1)).col(i-1) );
//   return el + eu;
// }
// [[Rcpp::export(getsigma_dLocal_RCPP)]] //getsigma_dLocal 
mat getsigma_dLocal(mat eps, int k, int p, int G){ //CURRENTLY WORKS getsigma_dLocal_RCPP( as.matrix(cbind(rnorm(1000),rnorm(1000,0,2),rnorm(1000,0,2), rnorm(1000))), 500, 1, 450 )
  mat sigma_d = (cov(eps.rows(k,k+G-1), 1) + cov(eps.rows(k-G,k-1),1)) /2 ; // should this be /(2)
  //mat upper = eps.rows(k,k+G-1); mat lower = eps.rows(k-G,k-1);
  //mat upper_c = upper - mean(upper,0); mat lower_c = lower- mean(lower,0);
 // mat sigma_d = (upper_c.t() * upper_c + lower_c.t() * lower_c) / (2*G);
  return sigma_d;
}



// [[Rcpp::export(get_DiagC_RCPP)]] //get_DiagC
mat DiagC(mat x, int p, mat sigma_d, int k, int G)
{
  int n = x.n_rows;
  int d = x.n_cols;
  mat xk;vec evals; mat evecs; vec evalS; mat evecS;
  mat O = ones(2*G,1) ;
  if(p==1) xk = join_rows(O, x.rows( k-G-1,k+G-1-1) );
  if(p==2) xk = join_rows(O, x.rows( k-G-1,k+G-1-1) , x.rows( k-G-1-1,k+G-2-1));
  if(p==3) xk = join_rows(O, x.rows( k-G-1,k+G-1-1) , x.rows( k-G-1-1,k+G-2-1),x.rows( k-G-2-1,k+G-3-1));
  if(p==4) {
    xk = join_rows(x.rows( k-G-1,k+G-1-1) , x.rows( k-G-1-1,k+G-2-1),x.rows( k-G-2-1,k+G-3-1), x.rows( k-G-3-1,k+G-4-1));
    xk.insert_cols(0,1);
  };
  if(p==5) {
    xk = join_rows(x.rows( k-G-1-1,k+G-2-1),x.rows( k-G-2-1,k+G-3-1), x.rows( k-G-3-1,k+G-4-1),x.rows( k-G-4-1,k+G-5-1));
    xk.insert_cols(0, join_rows(O,x.rows( k-G-1,k+G-1-1)));
  };
  if(p==6) {
    xk = join_rows(x.rows( k-G-2-1,k+G-3-1), x.rows( k-G-3-1,k+G-4-1),x.rows( k-G-4-1,k+G-5-1),x.rows( k-G-5-1,k+G-6-1));
    xk.insert_cols(0, join_rows(O,x.rows( k-G-1,k+G-1-1), x.rows( k-G-1-1,k+G-2-1) ));
  };
  mat C =  xk.t() * xk /(2*G -1); //works as intended
  //eigen decomposition
  eig_sym(evals,evecs,C);
  mat C_ = evecs * diagmat( 1/(evals) )*  evecs.t(); //now returns inverse Sigma
  eig_sym(evalS,evecS,sigma_d);
  mat S_ = evecS * diagmat( 1/(evalS) )*  evecS.t();

  //mat evecKron = kron(evecs, evecS);
  
  //eig_sym(evals,evecs,kron(sigma_d,C));
  //mat out = evecs * diagmat( pow(evals, -0.5) )*  evecs.t();
  //mat out = evecKron * kron(diagmat( pow(evals, -0.5) ), diagmat( pow(evalS, -0.5) )) * evecKron.t();
  mat out = kron(S_, C_);
  return out;
}

// [[Rcpp::export(get_DiagC_univ)]] //get_DiagC
mat DiagC_univ(mat x, int p, mat sigma_d, int k, int G)
{
  int n = x.n_rows;
  int d = x.n_cols;
  mat xk(2*G,d*p+d, fill::ones) ;
  //mat O = ones(2*G,1) ;
  for(int t=0; t<2*G ; t++){
    for(int ii=0; ii<d; ii++){
      xk(t, span( (ii)*(p+1) + 1, (ii)*(p+1) +p ) ) = x( span(k- G + t,k- G + t+p-1),  ii ).t()  ;
    }
  }
  
  mat C =  xk.t() * xk /(2*G -1); //works as intended
  mat S = repelem(sigma_d, p+1, p+1);
  mat out = pinv(S % C);
  return out;
}


// [[Rcpp::export(get_Tkn_RCPP)]] //get_Tkn
double Tkn(mat x, int k, int p,  int G, mat Phi, mat eps, mat h_all , String estim, String var_estim, mat sgd, bool univariate = 0)
{
  int n = x.n_rows;
  int d = x.n_cols;
  double out;
  mat A = getA(x,k,G,p,eps,h_all);
  //double v;
  //Sigma estimator options------
  if(estim == "DiagC"){
    //if(var_estim == "Local"){sgd = conv_to<vec>::from(sigma_d.row(k));} else if (var_estim == "Global") {sgd = conv_to<vec>::from(sigma_d);}
    if (var_estim == "Local") {sgd = getsigma_dLocal(eps, k, p, G);}
    mat Sig_;
    if(univariate){
      Sig_ = DiagC_univ(x,p,sgd,k,G) ;
    } else {
      Sig_ = DiagC(x,p,sgd,k,G) ;
    }
     
    mat prod = A.t() * Sig_ * A;
    out = sqrt( prod(0,0) ) /sqrt(16*G) ;//out = norm(Sig_ * A)  WHY is it scaling like this? undoing estimator scaling?
  } else{
    if(estim == "FullH") {
      mat Sig_ = FullH(x,k,G,h_all);  //FullH estimator
      out = norm(Sig_ * A) / (sqrt(8*G)); //2*
    }
  }
  //------------------------------
  return out;
}


// [[Rcpp::export(get_Tkn_RCPP)]] //get_Tkn
double Tkn_bootstrap(mat x, int k, int p,  int G, mat Phi, mat eps, mat h_all , cube DCcube)
{
  int n = x.n_rows;
  int d = x.n_cols;
  double out;
  mat A = getA(x,k,G,p,eps,h_all);

  //Sigma estimator options------
    mat Sig_ = DCcube.slice(k);
    mat prod = A.t() * Sig_ * A;
    out = sqrt( prod(0,0) ) /sqrt(16*G) ;//out = norm(Sig_ * A)  WHY is it scaling like this? undoing estimator scaling?
  return out;
}




// double floop(int k) {
//   mat x; int p; int G; String estim;
//   return Wkn(x,p,k,G,estim);
// } //wrapper
// 

// [[Rcpp::export(get_T_RCPP)]] //get_T
vec T(mat x, int p, int G, mat Phi, mat eps,arma::field<mat> PhiList, String estim = "DiagC",  String var_estim = "Local", bool univariate = 0)
{
  int n = x.n_rows;
  mat h_all;
  if(univariate){
    h_all  = H_all_univ(x, p, G, PhiList, eps);
  } else {
    h_all  = H_all(x, p, G, Phi, eps);
  }
  
 // mat sigma_d;
  vec out(n, fill::zeros);
  // if(var_estim == "Local"){
  //   sigma_d = getsigma_dLocal(eps,p,G);
  // } else if(var_estim == "Global"){ 
  //     sigma_d = getsigma_dGlobal(eps,p);
  // }
  mat sgd;
  if (var_estim == "Global") sgd = cov(eps.rows(p+1,n-1) );
  for (int k=G+ p+1; k< n-G-p; k++) {
    out(k) = Tkn(x,k,p,G,Phi,eps,h_all ,estim, var_estim, sgd, univariate);
  }
  return out;
}


// [[Rcpp::export(get_T_multiplier)]] //get_T
vec T_multiplier(mat x, int p, int G, mat Phi, mat eps, mat h_all, cube DCcube)
{
  int n = x.n_rows;
  vec out(n, fill::zeros);
  for (int k=G+ p+1; k< n-G-p; k++) {
    out(k) = Tkn_bootstrap( x,  k,  p,   G,  Phi,  eps,  h_all ,  DCcube);
  }
  return out;
}


vec cps(vec Wn, double D_n, int G, double nu = 0.25)
{
  int n = Wn.size(); vec out;
  //rshift <- c(Tn[-1],0); lshift <- c(0,Tn[-n]); 
  vec rshift = shift(Wn,-1); vec lshift = shift(Wn,1); 
  uvec over = find(Wn >D_n); //indices are greater than D_n?
  uvec lunder = find(lshift < D_n);
  uvec v = intersect(over, lunder); //lowers
  
  uvec runder = find(rshift < D_n);
  uvec w = intersect(over, runder); //uppers
  uvec nu_remove = find(w-v >= nu*G); //(epsilon) test for distance between 
  uvec v_nu = v.elem(nu_remove); uvec w_nu = w.elem(nu_remove); //w_nu.insert_rows(w_nu.size(),1); w_nu(w_nu.size()-1)=n;
  int q = nu_remove.size(); //number of CPs
  if(q>0){
    out = zeros(q);
    for (int ii=0; ii< q; ii++) {
      out(ii) = v_nu(ii) + Wn( span(v_nu(ii)-1,w_nu(ii)-1)).index_max() ;
    };
  };   
  return out;
}

// [[Rcpp::export(test_Score_RCPP)]] //Score
List test_Score(mat x, int p, int G, mat Phi, mat eps, double alpha =0.05, String estim = "DiagC"){
  int n = x.n_rows;
  int d = x.n_cols;
  vec cp; //double nu=1/4;
  //Test setup----------------------------
  double c_alpha = -log(log( pow((1-alpha),-0.5)) ); //critical value
  double a = sqrt(2*log(n/G)); //test transform multipliers
  double g = lgamma(d*(d*p+1)/2); double l23 = log(2)-log(3);
  double b = 2*log(n/G) + (d*(d*p+1)/2) * log(log(n/G)) -g - l23   ;
  // D_n =  ;//threshold
  double D_n = max((b + c_alpha)/a, sqrt(2*log(n)) + c_alpha/sqrt(2*log(n)) ); //##ASYMPTOTIC correction
  int Reject = 0; //
  //Run test-----------------------------
  vec Tn = T(x,p,G, Phi, eps, arma::field<arma::mat>(0), estim); //evaluate statistic at each time k
  double test_stat = Tn.max();
  if(test_stat > D_n){ //compare test stat with threshold
    Reject = 1; //test true
    cp = cps(Tn,D_n,G);
    if( cp.size()==0 ) Reject = 0 ;//doesn't pass eps-test
  } ;
  //Plot------------------------------------
  //       plot(Wn) # plot test statistic
  //         abline(h = D_n, col = "blue") #add threshold
  //         if(Reject==TRUE) abline(v = cps, col = "red")  #if rejecting H0, add estimated cps
  //           pl <- recordPlot()
  // #plot( a*Tn - b); abline(h=c_alpha, col="blue") #rescaled plot
  //Output------------------------------------
  List out = List::create(Named("Reject") = Reject,  _["ChangePoints"] = cp);//, _["D_n"]=D_n, _["Wn"] = Wn,);
  return out ;
}

// [[Rcpp::export(MFA_Score)]] 
List MFA(mat x, int p, vec Gset, mat Phi, mat eps, String estim = "DiagC", double alpha = 0.05){
  Gset = sort(Gset); //into ascending order
  int Glen = Gset.size();
  NumericVector cps;
  bool Reject = FALSE;
  NumericVector v(Glen); List tests(Glen);
  List t;
  for(int ii =0; ii < Glen; ii++){
    t =  test_Score(x, p, Gset[ii], Phi, eps, alpha, estim);
    tests[ii] = t;
    if(t["Reject"]){ Reject = TRUE;}
  }
  if(Reject){
    cps = as<List>(tests[0])["ChangePoints"];
    if(Glen > 1){
      for(int ii=1;  ii<Glen; ii++){
        NumericVector K = as<List>(tests[ii])["ChangePoints"];
        for(int jj=0; jj < K.size(); jj++) {
          if(min(abs(K[jj] - cps)) > Gset[ii]) {cps.push_back(K[jj]);}
        }
      }
    }
  } //list(Reject, cps, q= length(cps))
  return(List::create(Named("Reject") = Reject, _["ChangePoints"]=cps, _["q"] = cps.size() ));
}



// [[Rcpp::export]]
arma::uvec sample_index(const int &size, const int &n_tilde){
  // return uvec containing size-many samples from 1:n_tilde
  arma::uvec sequence = arma::linspace<arma::uvec>(0, n_tilde-1, n_tilde);
  arma::uvec out = Rcpp::RcppArmadillo::sample(sequence, size, true);
  return out;
}

// [[Rcpp::export(bootstrap)]]
vec bootstrap(mat x, int p, int G, field<mat> bsPhi, mat bsResid, int n_tilde, int M, String estim, String var_estim, bool univariate = 1){
  int n = x.n_rows;
  int d = x.n_cols;
  uvec s;  uvec i;
  vec stat_m(n); vec max_m(M);
  mat epsi;
  
  for(int m = 0; m < M; m++){
    s = sample_index(n, n_tilde); //randi(n, distr_param(0, n_tilde-1));
    // stat_m = zeros(n);
    // for(int ii = 0; ii < d; ii++ ){
    //   i = ii; //uvec conversion
       epsi = bsResid.rows(s); //bsResid.submat(s, i);
    //   mat Phi = bsPhi(ii);
    //   stat_m = stat_m + T(x.col(ii), p, G, Phi, epsi, arma::field<mat>(0), estim, var_estim);
    // }
  stat_m = T(x, p, G, bsPhi(0), epsi, bsPhi, estim, var_estim, univariate)  ;
  max_m(m) = max(stat_m);
  }

  return max_m;
}

// [[Rcpp::export(multiplier_bootstrap)]]
vec multiplier_bootstrap(mat x, int p, int G, field<mat> PhiList, mat eps, vec cps, int L, int M, String estim, String var_estim){
  int n = x.n_rows;
  int d = x.n_cols;
  vec stat_m(n); vec max_m(M); 
  vec perturb; vec perturb_new;
  int K = n/L; int remainder = n - K*L;
  vec scaled_cps =  cps  * ( (double) L/n ) ; 
  vec lshift = scaled_cps - (double) 1.0; vec rshift = scaled_cps + (double) 1.0;
  scaled_cps = join_cols(lshift, scaled_cps, rshift); //adjacent blocks to 0
  uvec u_cps = conv_to<uvec>::from(scaled_cps);

  cube DCcube(d*p+d, d*p+d, n);
  for(int k=G+ p+1; k< n-G-p; k++) {
    mat sgd = getsigma_dLocal(eps, k, p, G);
    DCcube.slice(k) = DiagC_univ( x,  p,  sgd,  k,  G);
  }
  
  mat h_all = H_all_univ(x,p,G,PhiList,eps);
  mat h_temp;
  for(int m = 0; m < M; m++){
    perturb = randn(L);
    perturb(u_cps) = zeros(u_cps.n_elem);
    perturb_new = repelem(perturb, K,1);
    if(remainder>0) perturb_new.insert_rows(L*K, remainder ); //pad
    h_temp = h_all.each_row() % perturb_new.t();
    stat_m = T_multiplier(x, p, G, PhiList(0),  eps, h_temp, DCcube)  ;
    max_m(m) = max(stat_m);
  }
  
  return max_m;
}
