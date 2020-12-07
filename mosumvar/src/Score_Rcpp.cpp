
//#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;


// [[Rcpp::export(getH_ik_RCPP)]] //getH_ik
arma::vec H_ik(arma::mat& x, int i, int k, int p, arma::mat Phi, arma::mat eps)
{
  int n = x.n_rows;
  int d = x.n_cols;
  arma::vec ai = Phi.row(i-1).t(); //intercept - index

  arma::mat X = x.rows( span(k-p-1,k-1-1)) ;

  arma::vec V = arma::vectorise(flipud(X),1).t();
  arma::vec O = ones(1);
  arma::vec VV = join_cols(O,V);

  double y = x(k-1,i-1); //double aV = dot(ai, VV);
  double  e = eps(k-1,i-1);  // y - aV;  //residual

  arma::vec  out = - y*VV +  VV*VV.t()*ai - e*VV ;

  return out;

}


// [[Rcpp::export(makeH_k_RCPP)]] //makeH_k
arma::vec H_k(arma::mat& x, int k, int p, arma::mat Phi, arma::mat eps)
{
  int d = x.n_cols;
  arma::vec H = zeros(d+  d*d*p); //accounts for intercept

  for (int ii=1; ii< d+1; ii++) {
    H.subvec((ii-1)*(d*p+1)+1-1, (ii-1)*(d*p+1)+ d*p +1-1 )= H_ik(x,ii,k,p,Phi,eps) ;
  };
  return H;
}


// [[Rcpp::export(makeH_all_RCPP)]] //makeH_all
arma::mat H_all(arma::mat& x,  int p, int G, arma::mat Phi, arma::mat eps)
{
  int n = x.n_rows;
  int d = x.n_cols;
  int nr = d+ d*d*p;
  //int nc = u-l+1;
  arma::mat H; H.zeros(nr,n); //arma::matrix of H values #accounts for intercept
  for (int t=p+1; t <(n-p); t++ ) {
    H.col(t) = H_k(x, t, p, Phi, eps) ;//-1-1
  };
  return H;

}

// struct arma::matrix_types {
//   arma::arma::mat m;
//   arma::sp_arma::mat M;
// };
// // [[Rcpp::export]] // INTERLEAVE rows of input list
// arma::arma::mat write_rows2(Rcpp::List data,  int nrows, int ncols) {//Rcpp::Characterarma::vector clss,
//
//   const int len = data.length();
//   std::arma::vector<arma::matrix_types> arma::matr(len);
//   std::arma::vector<bool> is_dense(len);
//   arma::arma::mat result(nrows*len, ncols);
//
//   // populate the structs
//   for (int j = 0; j < len; j++) {
//     //is_dense[j] = (clss[j] == "arma::matrix");
//     //if (is_dense[j]) {
//     arma::matr[j].m = Rcpp::as<arma::arma::mat>(data[j]);
//     //}
//     //else {
//     //  arma::matr[j].M = Rcpp::as<arma::sp_arma::mat>(data[j]);//SPARSE M
//     //}
//   }
//
//   // populate the result
//   for (int i = 0, k = 0; i < nrows; i++) {
//     for (int j = 0; j < len; j++, k++) {
//       //if (is_dense[j]) {
//       result.row(k) = arma::matr[j].m.row(i);
//       //}
//       //else {
//       //  arma::rowarma::vec r(arma::matr[j].M.row(i));
//       //  result.row(k) = r;
//       //}
//     }
//   }
//   return result;
// }


// // [[Rcpp::export]]
arma::sp_mat blockDiag( arma::field<arma::mat>& xlist ) {

  //xlist: list of arma::matrices

  unsigned int n = xlist.n_rows ;
  int dimen = 0 ;
  arma::vec dimvec(n) ;

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
// double sigma_i_k(arma::mat x, int i,int k,int G,int p,arma::vec a_upper,arma::vec a_lower)
// {
//   arma::mat x_upper; arma::mat x_lower;
//   arma::mat O = ones(G,1) ;
//   if(p==1) x_upper = join_rows(O, x.rows( k-1,k+G-1-1) ); //upper sample
//   if(p==2) x_upper = join_rows(O, x.rows( k-1,k+G-1-1), x.rows(k-1-1, k+G-2-1 ) );
//   if(p==3) x_upper = join_rows(O, x.rows( k-1,k+G-1-1), x.rows(k-1-1, k+G-2-1 ), x.rows(k-2-1, k+G-3-1 ) );
//   if(p==4){
//     x_upper = join_rows(x.rows( k-1,k+G-1-1), x.rows(k-1-1, k+G-2-1 ), x.rows(k-2-1, k+G-3-1 ),x.rows(k-3-1, k+G-4-1 ) );
//     x_upper.insert_cols(0,1);
//   } ;
//   rowarma::vec res_upper =  x( span(k+1-1,k+G-1), i-1).t() - a_upper.t() * x_upper.t(); //upper residuals
//
//   if(p==1)x_lower = join_rows(O, x.rows( k-G-1,k-1-1) );
//   if(p==2)x_lower  = join_rows(O, x.rows( k-G-1,k-1-1), x.rows(k-G-2, k-2-1));
//   if(p==3)x_lower  = join_rows(O, x.rows( k-G-1,k-1-1), x.rows(k-G-2, k-2-1), x.rows(k-G-3, k-3-1));
//   if(p==4){
//     x_lower = join_rows(x.rows( k-G-1,k-1-1), x.rows(k-G-2, k-2-1), x.rows(k-G-3, k-3-1),x.rows(k-G-4, k-4-1) );
//     x_lower.insert_cols(0,1);
//   };
//   rowarma::vec res_lower =  x( span(k-G+1-1,k-1), i-1).t() - a_lower.t() * x_lower.t(); //lower residuals
//   double sigma_i =  (sum(square(res_upper)) + sum(square(res_lower)) ) /(2*G);
//   return sigma_i;
//
// }

// // [[Rcpp::export(getsigma_d_kLOCAL1_RCPP)]] //getsigma_d_kLOCAL1
// arma::vec sigma_d_k(arma::mat x, int k,int G,int p,arma::vec a_upper,arma::vec a_lower)
// {
//   int d = x.n_cols;
//   arma::vec sigma_d = zeros(d);
//   for (int ii=1; ii< d+1; ii++) {
//     sigma_d(ii-1) = sigma_i_k(x,ii,k,G,p, a_upper( span((ii-1)*(d*p+1)+1-1, (ii-1)*(d*p+1)+ d*p +1-1)),
//             a_lower( span((ii-1)*(d*p+1)+1-1, (ii-1)*(d*p+1)+ d*p +1-1)) ) ;
//   };
//   return sigma_d;
// }


// // SIGMA ESTIarma::matORS ------------------------

// [[Rcpp::export(get_DiagH_RCPP)]] //get_DiagH
arma::sp_mat DiagH(arma::mat x, int k, int G,int p, arma::mat h_all)//sp_arma::mat
{
  int n = x.n_rows;
  int d = x.n_cols;
  field<arma::mat> H_list(d);
  arma::mat H_u_i; arma::mat H_l_i; arma::mat H_out;
  arma::vec Hbar_u; arma::vec Hbar_l;
  arma::vec evals; arma::mat evecs;
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
  arma::sp_mat Sig_ =  sqrt(2*G) * blockDiag(H_list); //coerce into block diagonal form
  return Sig_;
}

//

// [[Rcpp::export(get_FullH_RCPP)]] //get_FullH
arma::mat FullH(arma::mat x, int k, int G, arma::mat h_all)
{
  int n = x.n_rows;
  int d = x.n_cols;
  //field<arma::mat> H_list(d);
  arma::mat H_u = h_all.cols((k-1),(k+G-1-1) );
  arma::mat H_l = h_all.cols((k-G-1),(k-1-1) );
  arma::mat H_out;
  arma::vec Hbar_u; arma::vec Hbar_l;
  arma::vec evals; arma::mat evecs;

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
  // //field<arma::mat> H_list(d);
  // arma::mat H_; arma::mat<long double> H_out;
  // arma::vec Hbar_u; arma::vec Hbar_l;
  // arma::vec evals; arma::mat earma::vecs; std::arma::vector<long double> evld;
  //
  // Hbar_u = mean(H_u,1) ;
  // H_u.each_col() -= Hbar_u ;//H_u_centred
  // Hbar_l = mean(H_l,1) ;
  // H_l.each_col() -= Hbar_l; //H_l_centred
  // H_ = H_l * H_l.t() + H_u * H_u.t();
  // //eigen decomposition
  // eig_sym(evals,earma::vecs,H_);
  // evld = conv_to<std::arma::vector<long double>>::from(evals);
  // //Col<long double> pp = conv_to<Col<long double>>::from(powl(evld/(2*G),-0.5) );
  // arma::mat<long double> D = diagarma::mat()evld;
  // H_out = earma::vecs * (D)*  earma::vecs.t();
  return H_out; // RETURNS NaNs
}

// [[Rcpp::export(getA_RCPP)]] //getA
arma::vec getA(arma::mat x, int k, int G,int p, arma::mat eps, arma::mat h_all)
{
  int d = x.n_cols;
  arma::mat r = h_all.cols( (k+1-1),(k+G-1) );
  arma::mat l = h_all.cols( (k-G+1-1), (k-1) ); //left/right of window

  arma::vec A = sum(r,1) - sum(l,1) ;
  return A;
}

// sigma_i options---------------------------------------------

// [[Rcpp::export(getsigma_iGlobal_RCPP)]] //getsigma_iGlobal
double getsigma_iGlobal(arma::mat eps, int p, int i){
  int n = eps.n_rows;
  arma::vec out = square(eps.col(i-1).rows(p,n-1)) ;
  return mean(out);//mean(out);
}
// [[Rcpp::export(getsigma_dGlobal_RCPP)]] //getsigma_dGlobal
arma::mat getsigma_dGlobal(arma::mat eps, int p){
  // int n = eps.n_rows;
  // int d = eps.n_cols;
  // arma::vec out(d);
  // for(int i = 1; i<d+1; i++){
  //   out(i-1) = getsigma_iGlobal(eps,p,i);
  // }
  // return conv_to<arma::mat>::from(out);//mean(out);
  return cov(eps);
}


// // [[Rcpp::export(getsigma_iLocal_RCPP)]] //getsigma_iLocal
// double getsigma_iLocal(arma::mat eps, int i, int k, int G){
//   double el = var( eps.rows( (k-G+1-1),(k-1)).col(i-1)) ;
//   double eu = var( eps.rows( (k+1-1),(k+G-1)).col(i-1) );
//   return el + eu;
// }
// [[Rcpp::export(getsigma_dLocal_RCPP)]] //getsigma_dLocal
arma::mat getsigma_dLocal(arma::mat eps, int k, int p, int G){ //CURRENTLY WORKS getsigma_dLocal_RCPP( as.arma::matrix(cbind(rnorm(1000),rnorm(1000,0,2),rnorm(1000,0,2), rnorm(1000))), 500, 1, 450 )
  arma::mat sigma_d = (cov(eps.rows(k,k+G-1), 1) + cov(eps.rows(k-G,k-1),1)) /2 ; // should this be /(2)
  //arma::mat upper = eps.rows(k,k+G-1); arma::mat lower = eps.rows(k-G,k-1);
  //arma::mat upper_c = upper - mean(upper,0); arma::mat lower_c = lower- mean(lower,0);
 // arma::mat sigma_d = (upper_c.t() * upper_c + lower_c.t() * lower_c) / (2*G);
  return sigma_d;
}



// [[Rcpp::export(get_DiagC_RCPP)]] //get_DiagC
arma::mat DiagC(arma::mat x, int p, arma::mat sigma_d, int k, int G)
{
  int n = x.n_rows;
  int d = x.n_cols;
  arma::mat xk;arma::vec evals; arma::mat evecs; arma::vec evalS; arma::mat evecS;
  arma::mat O = ones(2*G,1) ;
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
  arma::mat C =  xk.t() * xk /(2*G -1); //works as intended
  //eigen decomposition
  eig_sym(evals,evecs,C);
  arma::mat C_ = evecs * diagmat( 1/(evals) )*  evecs.t(); //now returns inverse Sigma
  eig_sym(evalS,evecS,sigma_d);
  arma::mat S_ = evecS * diagmat( 1/(evalS) )*  evecS.t();

  //arma::mat earma::vecKron = kron(earma::vecs, earma::vecS);

  //eig_sym(evals,earma::vecs,kron(sigma_d,C));
  //arma::mat out = earma::vecs * diagarma::mat( pow(evals, -0.5) )*  earma::vecs.t();
  //arma::mat out = earma::vecKron * kron(diagarma::mat( pow(evals, -0.5) ), diagarma::mat( pow(evalS, -0.5) )) * earma::vecKron.t();
  arma::mat out = kron(S_, C_);
  return out;
}


// [[Rcpp::export(get_Tkn_RCPP)]] //get_Tkn
double Tkn(arma::mat x, int k, int p,  int G, arma::mat Phi, arma::mat eps, arma::mat h_all , String estim, String var_estim)
{
  int n = x.n_rows;
  int d = x.n_cols;
  double out;
  arma::mat A = getA(x,k,G,p,eps,h_all);
  //double v;
  //Sigma estiarma::mator options------
  if(estim == "DiagC"){
    arma::mat sgd;
    if (var_estim == "Global") sgd = cov(eps.rows(p+1,n-1) );
    //if(var_estim == "Local"){sgd = conv_to<arma::vec>::from(sigma_d.row(k));} else if (var_estim == "Global") {sgd = conv_to<arma::vec>::from(sigma_d);}
    if (var_estim == "Local") {
      sgd = getsigma_dLocal(eps, k, p, G);}
      arma::mat Sig_ = DiagC(x,p,sgd,k,G) ;
      arma::mat prod = A.t() * Sig_ * A;
      out = sqrt( prod(0,0) ) /sqrt(16*G) ;//out = norm(Sig_ * A)  WHY is it scaling like this? undoing estiarma::mator scaling?
    } else{
   // if(estim == "DiagH") {
   //   sp_arma::mat Sig_ = DiagH(x,k,G,p,h_all); //DiagH estiarma::mator for Sigma
   //   out = pow(2*G, -0.5) * norm(Sig_ * A, "fro");
   // }
    if(estim == "FullH") {
      arma::mat Sig_ = FullH(x,k,G,h_all);  //FullH estiarma::mator
      out = norm(Sig_ * A) / (sqrt(8*G)); //2*
    }

  }
  //------------------------------
  return out;
}
//
//
// double floop(int k) {
//   arma::mat x; int p; int G; String estim;
//   return Wkn(x,p,k,G,estim);
// } //wrapper
//

// [[Rcpp::export(get_T_RCPP)]] //get_T
arma::vec T(arma::mat x, int p, int G, arma::mat Phi, arma::mat eps, String estim = "DiagC",  String var_estim = "Local")
{
  int n = x.n_rows;
  arma::mat h_all = H_all(x, p, G, Phi, eps);
 // arma::mat sigma_d;
  arma::vec out(n, fill::zeros);
  // if(var_estim == "Local"){
  //   sigma_d = getsigma_dLocal(eps,p,G);
  // } else if(var_estim == "Global"){
  //     sigma_d = getsigma_dGlobal(eps,p);
  // }

  for (int k=G+ p+1; k< n-G-p; k++) {
    out(k) = Tkn(x,k,p,G,Phi,eps,h_all ,estim, var_estim);
  }
  return out;
}



arma::vec cps(arma::vec Wn, double D_n, int G, double nu = 0.25)
{
  int n = Wn.size(); arma::vec out;
  //rshift <- c(Tn[-1],0); lshift <- c(0,Tn[-n]);
  arma::vec rshift = shift(Wn,-1); arma::vec lshift = shift(Wn,1);
  arma::uvec over = find(Wn >D_n); //indices are greater than D_n?
  arma::uvec lunder = find(lshift < D_n);
  arma::uvec v = intersect(over, lunder); //lowers

  arma::uvec runder = find(rshift < D_n);
  arma::uvec w = intersect(over, runder); //uppers
  arma::uvec nu_remove = find(w-v >= nu*G); //(epsilon) test for distance between
  arma::uvec v_nu = v.elem(nu_remove); arma::uvec w_nu = w.elem(nu_remove); //w_nu.insert_rows(w_nu.size(),1); w_nu(w_nu.size()-1)=n;
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
List test_Score(arma::mat x, int p, int G, arma::mat Phi, arma::mat eps, double alpha =0.05, String estim = "DiagC"){
  int n = x.n_rows;
  int d = x.n_cols;
  arma::vec cp; //double nu=1/4;
  //Test setup----------------------------
  double c_alpha = -log(log( pow((1-alpha),-0.5)) ); //critical value
  double a = sqrt(2*log(n/G)); //test transform multipliers
  double g = lgamma(d*(d*p+1)/2); double l23 = log(2)-log(3);
  double b = 2*log(n/G) + (d*(d*p+1)/2) * log(log(n/G)) -g - l23   ;
  // D_n =  ;//threshold
  double D_n = max((b + c_alpha)/a, sqrt(2*log(n)) + c_alpha/sqrt(2*log(n)) ); //##ASYMPTOTIC correction
  int Reject = 0; //
  //Run test-----------------------------
  arma::vec Tn = T(x,p,G, Phi, eps, estim); //evaluate statistic at each time k
  double test_stat = Tn.max();
  if(test_stat > D_n){ //compare test stat with threshold
    Reject = 1; //test true
    cp = cps(Tn,D_n,G);
    if( cp.size()==0 ) Reject = 0 ;//doesn't pass eps-test
  } ;
  //Plot------------------------------------
  //       plot(Wn) # plot test statistic
  //         abline(h = D_n, col = "blue") #add threshold
  //         if(Reject==TRUE) abline(v = cps, col = "red")  #if rejecting H0, add estiarma::mated cps
  //           pl <- recordPlot()
  // #plot( a*Tn - b); abline(h=c_alpha, col="blue") #rescaled plot
  //Output------------------------------------
  List out = List::create(Named("Reject") = Reject,  _["ChangePoints"] = cp);//, _["D_n"]=D_n, _["Wn"] = Wn,);
  return out ;
}

// [[Rcpp::export(MFA_Score)]]
List MFA_Score(arma::mat x, int p, arma::vec Gset, arma::mat Phi, arma::mat eps, String estim = "DiagC", double alpha = 0.05){
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




