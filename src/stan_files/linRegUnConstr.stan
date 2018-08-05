data {
  int N; //the number of observations
  int K; //the number of columns in the model matrix
  int K1; //Intercept yes/no?
  vector[N] y; //the response
  vector[N] yUncertainty; //sd of uncertainties in Y
  matrix[N,K] X; //the model matrix
}
parameters {
  vector[K1] beta0; //the intercept
  vector[K - K1] beta; //the regression parameters
  real logsigma; //the standard deviation
}
transformed parameters{
  vector[N] sigma;
  sigma = sqrt(square(exp(logsigma)) + square(yUncertainty)) ;
}
model {  
  logsigma ~ student_t(3, 0, 1);
  if(K1 > 0){
    beta0 ~ student_t(1, 0, 10); //prior for the intercept
  }
    beta ~ student_t(3, 0, 5); //prior for the slopes
  y ~ normal(X * append_row(beta0, beta), sigma);
}
generated quantities{
  vector[N] log_lik;
  vector[N] sigmasq;
  real rsq;
  for (n in 1:N) sigmasq[n] = square(sigma[n]);
  rsq = 1 - mean(sigmasq) / variance(y); 
  for (n in 1:N) log_lik[n] = normal_lpdf(y[n] | X[n, ] * append_row(beta0, beta), sigma[n]);
}
