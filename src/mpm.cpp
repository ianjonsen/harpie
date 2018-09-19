#include <TMB.hpp>

using namespace density;

template<class Type>
Type objective_function<Type>::operator() ()
  {
    
    DATA_MATRIX(x);                   // locations
    DATA_VECTOR(d);                   // dive metric
    DATA_INTEGER(A);                  // number of animals
    DATA_VECTOR(dt);                  // time difference between dive stats
    DATA_FACTOR(idx);                // cumsum of number of locations for each animal

    PARAMETER_VECTOR(lg);		          // Autocorrelation parameter for locations (link scale)
    PARAMETER_VECTOR(ll);		          // Autocorrelation parameter for dives (link scale)
    PARAMETER_VECTOR(log_sigma);	    // Innovation variance (log scale)
    PARAMETER(log_sigma_g);           // logistic scale parameter of rw on lg (log scale)
    PARAMETER(log_sigma_l);


    // Backtransform parameters from link scale
    vector<Type> gamma = Type(1.0) / (Type(1.0) + exp(-lg));
    vector<Type> lambda = Type(1.0) / (Type(1.0) + exp(-ll));
    vector<Type> sigma = exp(log_sigma);
    Type sigma_g = exp(log_sigma_g);
    Type sigma_l = exp(log_sigma_l);


    // 2x2 covariance matrix for innovations
    matrix<Type> cov(2,2);
    cov(0,0) = sigma(0) * sigma(0);
    cov(0,1) = 0.0;
    cov(1,0) = 0.0;
    cov(1,1) = sigma(1) * sigma(1);

    Type jnll = 0.0;
    vector<Type> mu(2);

    MVNORM_t<Type> nll_dens(cov);   // Multivariate Normal density
    int i,j;
         
    for(i = 0; i < A; ++i) {
      for(j = (idx(i)+1); j < idx(i+1); ++j) {
        jnll -= dnorm(lg(j), lg(j-1), sigma_g, TRUE);  // RW on logit(gamma)
      }
      for(j = 1; j < dt.size(); ++j) {
        jnll -= dnorm(ll(j), ll(j-1), dt(j) * sigma_l, TRUE);  // RW on logit(lambda)
      }
      
      for(j = (idx(i)+2); j < idx(i+1); ++j){
        mu = x.row(j) - x.row(j-1) - gamma(j-1) * (x.row(j-1) - x.row(j-2));  // first diff RW on locations
        jnll += nll_dens(mu);
      }
      for(j = 1; j < dt.size(); ++j) {
        jnll -= dnorm(d(j), lambda(j-1) * d(j-1), dt(j) * sigma(2), TRUE);
      }
    }

    ADREPORT(sigma_g);
    ADREPORT(sigma_l);
    ADREPORT(sigma);

    return jnll;
}
  

