#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

/* Given an open cover x and a function v, it computes the pushforward of v to
 * the open sets of x given by the average. It also performs perm random permutations
 * of the values of v and computes the corresponding pushforwards.
 * 
 * Adapted from Camara Lab's RayleighSelection R package
 */

// [[Rcpp::export]]
List pushCpp(arma::vec v, List x, SEXP perm) {
  List xlist(x);
  arma::vec xv(v);
  arma::vec co = xv;
  int n = xlist.size();
  int perms = as<int >(perm);
  arma::mat res_vertices = arma::zeros(perms+1, n);
  
  for (int i=0; i<perms+1; i++) {
    int kk = 0;
    for (int j=0; j<n; j++) {
      
      arma::ivec o1 = xlist[j];
      int u = o1.size();
      
      for (int k=0; k<u; k++) {
        res_vertices(i,j) += co(o1[k]-1);
      }
      
      res_vertices(i,j) /= u;
  
    }
    
    co = arma::shuffle(co);
  };
  
  return List::create(Named("vertices") = res_vertices);
}
