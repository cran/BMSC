data {
  int N; //the number of observations
  int K; //the number of columns in the model matrix
  int K1; //Intercept yes/no?
  vector[N] y; //the response
  vector<lower = 0>[N] yUncertainty; //sd of uncertainties in Y
  matrix<lower = 0>[N,K] xUncertaintyMatrix; // sd of uncertainties in X
  matrix[N,K] X; //the model matrix
}
parameters {
  vector[K1] beta0; //the intercept
  simplex[K - K1] beta; //the regression parameters
  real logsigma; //the standard deviation
  vector[K] mu_x; // prior location
  vector[K] log_sigma_x; // prior scale
  matrix[N,K] XTRUE; //the model matrix of True values
}
transformed parameters{
  vector[N] sigma;
  vector[K] sigma_x; // prior scale
  vector[K] betaAll;
  sigma = sqrt(square(exp(logsigma)) + square(yUncertainty)) ;
  sigma_x = exp(log_sigma_x);
  betaAll = append_row(beta0, beta);
}
model {  
  logsigma~ student_t(3, 0, 1);
  if(K1 > 0){
    beta0 ~ student_t(1, 0, 10); //prior for the intercept
  }
  beta ~ student_t(3, 0, 5); //prior for the slopes
  mu_x ~ student_t(1, 0, 10);
  log_sigma_x ~ student_t(3, 0, 5);
  for (i in 1:K){
    XTRUE[,i] ~ normal(mu_x[i], sigma_x[i]);
    X[,i] ~ normal(XTRUE[,i], xUncertaintyMatrix[,i] + 1E-4);
  }
  y ~ normal(XTRUE * betaAll, sigma);
}
generated quantities{
  vector[N] log_lik;
  vector[N] sigmasq;
  real rsq;
  for (n in 1:N) sigmasq[n] = square(sigma[n]);
  rsq = 1 - mean(sigmasq) / variance(y); 
  for (n in 1:N) log_lik[n] = normal_lpdf(y[n] | XTRUE[n, ] * betaAll, sigma[n]);
}
