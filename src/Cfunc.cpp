#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

//' Natural spline method
//' @param x: Attribute vector
//' @param knots: Vector of knots (default: all unique values of x)
//' @return: Natural Spline matrix N, each column represents a basis function, each row corresponds to a value of x
//' @export
// [[Rcpp::export]]
NumericMatrix NaturalSplineC(NumericVector& x) {
  NumericVector knots = unique(x);
  sort(knots.begin(), knots.end());
  NumericMatrix d(x.size(), knots.size() - 1);
  for (int i = 0; i < x.size(); i++) {
    for (int k = 0; k < knots.size() - 1; k++) {
      double nume1 = max(pow(x[i]-knots[k], 3), 0.0);
      double nume2 = max(pow(x[i]-knots[knots.size() - 1], 3), 0.0);
      double denom = knots[knots.size() - 1] - knots[k];
      d(i, k) = (nume1 + nume2) / denom;
    }
  }
  
  NumericMatrix N(x.size(), knots.size());
  NumericMatrix::Column c1 = N.column(0);
  c1 = c1 + 1;
  NumericMatrix::Column c2 = N.column(1);
  c2 = x;
  for (int k = 0; k < knots.size() - 2; k++) {
    NumericMatrix::Column Nc = N.column(k + 2);
    NumericMatrix::Column dc1 = d.column(k);
    NumericMatrix::Column dc2 = d.column(knots.size() - 2);
    Nc = dc1 - dc2;
  }
  
  return N;
}

// [[Rcpp::export]]
double calCoefDiffC(int p, NumericMatrix& X, NumericVector& y, double alpha_hat, 
                    NumericMatrix& f_hat, NumericVector& lambda, List& Coef, List& Coef_old, double Coef_diff) {
  
  Function SmoothingSpline_UnivariateC("SmoothingSpline_UnivariateC");
  Function norm("norm");
  for(int j = 0; j < p; j++) {
    NumericMatrix selected_f(f_hat.nrow(), f_hat.ncol()-1);
    for (int i = 0; i < j; i++) {
      selected_f(_, i) = f_hat(_, i);
    }
    for (int i = j+1; i < f_hat.ncol(); i++) {
      selected_f(_, i-1) = f_hat(_, i);
    }
    NumericVector yj = y - alpha_hat - rowSums(selected_f);
    NumericVector xj = X(_, j);
    Coef[j] = SmoothingSpline_UnivariateC(xj, yj, lambda[j]);
    f_hat(_, j) = NaturalSplineC(xj) * Coef[j];
    f_hat(_, j) = f_hat(_, j) - mean(f_hat(_, j));
    NumericVector coefj(Coef[j]);
    NumericVector nm( norm(coefj - Coef_old[j], "2") );
    Coef_diff = Coef_diff + nm[0];
  }
  return Coef_diff;
}


//' kernel matrix method
//' @param x: Attribute vector
//' @return: Kernel Matrix K
//' @export
// [[Rcpp::export]]
NumericMatrix KernelMatrixC(NumericVector& x) {
  NumericVector knots = unique(x);
  sort(knots.begin(), knots.end());
  int N = knots.size();
  NumericMatrix K(N, N);
  for (int i = 2; i < N; i++) {
    for (int j = i; j < N; j++) {
      double term1 = (12*(pow(knots[N-1], 3)-pow(knots[j-2], 3)) - 
                      18*(knots[i-2]+knots[j-2])*(pow(knots[N-1], 2)-pow(knots[j-2], 2)) +
                      36*knots[i-2]*knots[j-2]*(knots[N-1]-knots[j-2])) / 
                      ((knots[N-1]-knots[i-2])*(knots[N-1]-knots[j-2]));
      
      double term2 = (12*(pow(knots[N-1], 3)-pow(knots[N-2], 3)) - 
                      18*(knots[j-2]+knots[N-2])*(pow(knots[N-1], 2)-pow(knots[N-2], 2)) +
                      36*knots[j-2]*knots[N-2]*(knots[N-1]-knots[N-2])) / 
                      ((knots[N-1]-knots[N-2])*(knots[N-1]-knots[j-2]));
      
      double term3 = (12*(pow(knots[N-1], 3)-pow(knots[N-2], 3)) - 
                      18*(knots[i-2]+knots[N-2])*(pow(knots[N-1], 2)-pow(knots[N-2], 2)) +
                      36*knots[i-2]*knots[N-2]*(knots[N-1]-knots[N-2])) / 
                      ((knots[N-1]-knots[i-2])*(knots[N-1]-knots[N-2]));
      
      double term4 = (12*(pow(knots[N-1], 3)-pow(knots[N-2], 3)) -
                      36*knots[N-2]*(pow(knots[N-1], 2)-pow(knots[N-2], 2))+
                      36*pow(knots[N-2], 2)*(knots[N-1]-knots[N-2])) /
                        (pow((knots[N-1]-knots[N-2]), 2));
      
      K(i,j) = term1 - term2 - term3 + term4;
    }
  }
  
  NumericMatrix trans_K = transpose(K);
  NumericMatrix out(N, N);
  for (int i = 0; i < N; ++i) {
    NumericMatrix::Row r1 = K.row(i);
    NumericMatrix::Row r2 = trans_K.row(i);
    NumericMatrix::Row oRow = out.row(i);
    oRow = r1 + r2;
    out(i, i) -= K(i, i);
  }
  
  return out;
}