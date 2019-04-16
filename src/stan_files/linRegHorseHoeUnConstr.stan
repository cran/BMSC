data {
  int N; //the number of observations
  int K; //the number of columns in the model matrix
  int K2; //the number of columns in the model matrix 2 (non-regularized coefficients)
  int K1; //Intercept yes/no?
  real y[N]; //the response
  vector[N] yUncertainty; //sd of uncertainties in Y
  matrix[N, K] X; //the model matrix
  matrix[N, K2] X2; //the model matrix
}
parameters {
  vector[K1] beta0; //the intercept
  vector[K2] beta2; //the must include variables
  real logsigma; //the standard deviation
  vector[K-K1] zvalues;
  real <lower = 0> r1_global ;
  real <lower = 0> r2_global ;
  vector <lower = 0>[K-K1] r1_local ;
  vector <lower = 0>[K-K1] r2_local ;
}
transformed parameters {
  real <lower = 0> tau ; // global shrinkage parameter
  vector <lower = 0>[K-K1] lambda ; // local shrinkage parameters
  vector[K-K1] beta; //the regression parameters
  vector[K+K2] betaAll;
  vector[N] sigma;
  sigma = sqrt(square(exp(logsigma)) + square(yUncertainty)) ;
  lambda = r1_local .* sqrt ( r2_local );
  tau = r1_global * sqrt ( r2_global ) * sqrt(square(exp(logsigma)) + mean(square(yUncertainty)));
  beta = zvalues .* lambda * tau ;
  betaAll = append_row(append_row(beta0, beta), beta2);
}
model {
  logsigma ~ student_t(3, 0, 1);
  // half - t priors for lambdas
  zvalues ~ normal(0 , 1);
  r1_local ~ normal(0.0 , 1.0);
  r2_local ~ inv_gamma(0.5, 0.5);
  // half - t prior for tau
  r1_global ~ normal(0.0 , 1);
  r2_global ~ inv_gamma(0.5, 0.5);
  if(K1 > 0){
    beta0 ~ student_t(1, 0, 5);
  }
  if(K2 > 0){
  beta2 ~ student_t(3, 0, 5); //prior for the slopes
}
  y ~ normal(append_col(X, X2) * betaAll, sigma);
}
generated quantities{
  vector[N] log_lik;
  for (n in 1:N) log_lik[n] = normal_lpdf(y[n] | append_col(X[n, ], X2[n, ]) * betaAll, sigma[n]);
}
