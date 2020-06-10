
//#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(openmp)]]

#include <iostream>
#include <RcppParallel.h>
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


// [[Rcpp::export(makeH_k_RCPP)]] //makeH_k
vec H_k(mat& x, int k, int p, mat Phi, mat eps)
{
  int d = x.n_cols;
  vec H = zeros(d+  d*d*p); //accounts for intercept

  for (int ii=1; ii< d+1; ii++) {
    H.subvec((ii-1)*(d*p+1)+1-1, (ii-1)*(d*p+1)+ d*p +1-1 )= H_ik(x,ii,k,p,Phi,eps) ;
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
  int n = eps.n_rows;
  int d = eps.n_cols;
  vec out(d);
  for(int i = 1; i<d+1; i++){
    out(i-1) = getsigma_iGlobal(eps,p,i);
  }
  return conv_to<mat>::from(out);//mean(out);
}


// [[Rcpp::export(getsigma_iLocal_RCPP)]] //getsigma_iLocal
double getsigma_iLocal(mat eps, int i, int k, int G){
  double el = var( eps.rows( (k-G+1-1),(k-1)).col(i-1)) ; 
  double eu = var( eps.rows( (k+1-1),(k+G-1)).col(i-1) );
  return el + eu;
}
// [[Rcpp::export(getsigma_dLocal_RCPP)]] //getsigma_dLocal
mat getsigma_dLocal(mat eps, int p, int G){
  int n = eps.n_rows;
  int d = eps.n_cols;
  mat sigma_d(n,d);
  sigma_d.head_rows(p) = zeros(p,d); //fill NAs
  
    for(int t = (G+p+1); t<(n-G-1); t++) {
      for(int i=1; i<d+1; i++){
        sigma_d.row(t-1).col(i-1) = getsigma_iLocal(eps,i,t,G);
      };
    } ;
    return sigma_d;
}



// [[Rcpp::export(get_DiagC_RCPP)]] //get_DiagC
sp_mat DiagC(mat x, int p, vec sigma_d, int k, int G)
{
  int n = x.n_rows;
  int d = x.n_cols;
  mat xk;vec evals; mat evecs;
  mat O = ones(2*G,1) ;
  if(p==1) xk = join_rows(O, x.rows( k-G-1,k+G-1-1) );
  if(p==2) xk = join_rows(O, x.rows( k-G-1,k+G-1-1) , x.rows( k-G-1-1,k+G-2-1));
  if(p==3) xk = join_rows(O, x.rows( k-G-1,k+G-1-1) , x.rows( k-G-1-1,k+G-2-1),x.rows( k-G-2-1,k+G-3-1));
  if(p==4) {
    xk = join_rows(x.rows( k-G-1,k+G-1-1) , x.rows( k-G-1-1,k+G-2-1),x.rows( k-G-2-1,k+G-3-1), x.rows( k-G-3-1,k+G-4-1));
    xk.insert_cols(0,1);
  };
  mat C =  xk.t() * xk /(2*G);
  //eigen decomposition
  eig_sym(evals,evecs,C);
  mat C_ = evecs * diagmat( pow(evals, -0.5) )*  evecs.t();
  field<mat> Clist(d);
  for (int ii=0; ii< d; ii++) {
    Clist(ii) =   C_ /sqrt(2*sigma_d(ii));
  };
  sp_mat out = blockDiag(Clist); //coerce into block diagonal form

  return out;
}


// [[Rcpp::export(get_Tkn_RCPP)]] //get_Tkn
double Tkn(mat x, int k, int p,  int G, mat Phi, mat eps, mat h_all, mat sigma_d, String estim, String var_estim)
{
  int n = x.n_rows;
  int d = x.n_cols;
  double out;
  mat A = getA(x,k,G,p,eps,h_all);
  //double v;
  //Sigma estimator options------
  if(estim == "DiagC"){
    vec sgd;
    if(var_estim == "Local"){sgd = conv_to<vec>::from(sigma_d.row(k));} else if (var_estim == "Global") {sgd = conv_to<vec>::from(sigma_d);}
    sp_mat Sig_ = DiagC(x,p,sgd,k,G) ;
    out = pow(2*G, -0.5) * norm(Sig_ * A, "fro");
  } else{
    if(estim == "DiagH") {
      sp_mat Sig_ = DiagH(x,k,G,p,h_all); //DiagH estimator for Sigma
      out = pow(2*G, -0.5) * norm(Sig_ * A, "fro");
    }
    if(estim == "FullH") {
      mat Sig_ = FullH(x,k,G,h_all);  //FullH estimator
      out = pow(2*G, -0.5) * norm(Sig_ * A, "fro");
    }

  }
  //------------------------------
  return out;
}
// 
// 
// double floop(int k) {
//   mat x; int p; int G; String estim;
//   return Wkn(x,p,k,G,estim);
// } //wrapper
// 

// [[Rcpp::export(get_T_RCPP)]] //get_T
vec T(mat x, int p, int G, mat Phi, mat eps, String estim = "DiagC",  String var_estim = "Local")
{
  int n = x.n_rows;
  mat h_all = H_all(x, p, G, Phi, eps);
  mat sigma_d; vec out(n, fill::zeros);
  if(var_estim == "Local"){
    sigma_d = getsigma_dLocal(eps,p,G);
  } else if(var_estim == "Global"){ 
      sigma_d = getsigma_dGlobal(eps,p);
  }

  for (int k=G+ p+1; k< n-G-p; k++) {
    out(k) = Tkn(x,k,p,G,Phi,eps,h_all,sigma_d,estim, var_estim);
  }
  return out;
}


// 
// //SIMULATION -------------------------------------
// 
// 
// // [[Rcpp::export(get_cps_RCPP)]] //get_cps
// vec cps(vec Wn, double D_n, int G, double nu = 0.25)
// {
//   int n = Wn.size(); vec out;
//   //rshift <- c(Tn[-1],0); lshift <- c(0,Tn[-n]); 
//   vec rshift = shift(Wn,-1); vec lshift = shift(Wn,1); 
//   uvec over = find(Wn >D_n); //indices are greater than D_n?
//   uvec lunder = find(lshift < D_n);
//   uvec v = intersect(over, lunder); //lowers
//   
//   uvec runder = find(rshift < D_n);
//   uvec w = intersect(over, runder); //uppers
//   uvec nu_remove = find(w-v >= nu*G); //(epsilon) test for distance between 
//   uvec v_nu = v.elem(nu_remove); uvec w_nu = w.elem(nu_remove); //w_nu.insert_rows(w_nu.size(),1); w_nu(w_nu.size()-1)=n;
//   int q = nu_remove.size(); //number of CPs
//   if(q>0){
//     out = zeros(q);
//     for (int ii=0; ii< q; ii++) {
//       out(ii) = v_nu(ii) + Wn( span(v_nu(ii)-1,w_nu(ii)-1)).index_max() ;
//     };
//   };   
//   return out;
// }
// 
// // [[Rcpp::export(test_Wald_RCPP)]] //test_Wald
// List test_Wald(mat x, int p, int G, double alpha =0.05, String estim = "DiagC"){
//   int n = x.n_rows;
//   int d = x.n_cols;
//   vec cp; //double nu=1/4;
//   //Test setup----------------------------
//   double c_alpha = -log(log( pow((1-alpha),-0.5)) ); //critical value
//   double a = sqrt(2*log(n/G)); //test transform multipliers
//   double g = lgamma(d*(d*p+1)/2); double l23 = log(2)-log(3);
//   double b = 2*log(n/G) + (d*(d*p+1)/2) * log(log(n/G)) -g - l23   ;
//   // D_n =  ;//threshold
//   double D_n = max((b + c_alpha)/a, sqrt(2*log(n)) + c_alpha/sqrt(2*log(n)) ); //##ASYMPTOTIC correction
//   int Reject = 0; //
//   //Run test-----------------------------
//   vec Wn = W(x,p,G,estim); //evaluate statistic at each time k
//   double test_stat = Wn.max();
//   if(test_stat > D_n){ //compare test stat with threshold
//     Reject = 1; //test true
//     cp = cps(Wn,D_n,G);
//     if( cp.size()==0 ) Reject = 0 ;//doesn't pass eps-test
//   } ; 
//   //Plot------------------------------------
//   //       plot(Wn) # plot test statistic
//   //         abline(h = D_n, col = "blue") #add threshold
//   //         if(Reject==TRUE) abline(v = cps, col = "red")  #if rejecting H0, add estimated cps
//   //           pl <- recordPlot()
//   // #plot( a*Tn - b); abline(h=c_alpha, col="blue") #rescaled plot
//   //Output------------------------------------
//   List out = List::create(Named("Reject") = Reject,  _["ChangePoints"] = cp);//, _["D_n"]=D_n, _["Wn"] = Wn,);
//   return out ;
// }
// 
// 
// // [[Rcpp::export(sim_data_RCPP)]] 
// mat sim_data(List pars, int n=1000, int d=5, double sd = 0.2){
//   mat errors = randn(n,d)*sd;
//   mat simdata = errors;
//   int p = pars.length();
//   //simdata.rows(0,p-1) = errors.rows(0,p-1);
//   for(int r =p; r < n; r++){
//     //simdata.row(r)   errors[row,]
//     for (int ii=0; ii<p;  ii++) {
//       mat mp = pars[ii];
//       vec pp = mp * simdata.row(r-ii).t();
//       simdata.row(r) = simdata.row(r) +  pp.t() ; 
//     };
//   };
//   return simdata;
// }
// 
// 
// 
// // [[Rcpp::export(var_simulate_RCPP)]] 
// NumericVector var_sim(List pars, int reps =100, int p=2, int G=200, double alpha =0.05, String estim = "DiagC",int ncores =1){
//   vec cp={500,1000,1500};
//   NumericVector out(reps) ;
//   RcppParallel::RVector<double> wo(out);
//   //RcppParallel::RVector<double> wx(x);
//   
// #if defined(_OPENMP)
// #pragma omp parallel for num_threads(ncores)
// #endif
//   for(int repl = 0; repl < reps; repl++ ){
//     List p1 = List::create(pars[0], pars[1]);List p2 = List::create(pars[2], pars[3]);List p3 = List::create(pars[4], pars[5]);
//     mat r1 = sim_data(p1, cp(0), 5);mat r2 = sim_data(p2, cp(1)-cp(0), 5);mat r3 = sim_data(p3, cp(2)-cp(1), 5); //in-regime data
//     mat r = join_cols(r1,r2,r3); //full data
//     List t = test_Wald(r, p, G, alpha,estim);
//     wo[repl] = t[0];
//   };
//   //Output------------------------------------
//   // List out = List::create(Named("Reject") = Reject, _["Wn"] = Wn, _["ChangePoints"] = cp, _["D_n"]=D_n);
//   return out ;
// }


