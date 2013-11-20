data {
  int<lower=0> N;
  int<lower=0> N1;
  int<lower=0> N2;
  int<lower=0> N3;
  int<lower=0> N4;
  real<lower=0> y1[N1];
  real<lower=0> y2[N2];
  real<lower=0> y3[N3];
  real<lower=0> y4[N4];
  int<lower=0> obs1[N1];
  int<lower=0> obs2[N2];
  int<lower=0> obs3[N3];
  int<lower=0> obs4[N4]; 
  real<lower=0> b[4];
  real<lower=0> varsigma[4];
  real<lower=0> sigma[4];
  real eta[4];
  int<lower=0> P;
  real<lower=0> nu;
}
parameters {
  real<lower=0, upper=1> i0;
  real<lower=0, upper=1> x[2,N];
  real<lower=0> beta;
  real<lower=0> gamma;
}
model {
  i0 ~ normal(0.002, 0.0005) T[0,1];
  beta ~ lognormal(-1.3296, 0.1183);
  gamma ~ lognormal(-2.1764, 0.1055);
  x[1,1] ~ normal(1 - i0 - pow(beta*(i0)*(1-i0),nu), sqrt(beta) / P) T[0,1];
  x[2,1] ~ normal(i0*(1 - gamma) + (1-i0) - x[1,1], sqrt(gamma) / P) T[0,1-x[1,1]];
  for(j in 1:(N-1)) {
    x[1,j+1] ~ normal(x[1,j] - pow(beta*x[2,j]*x[1,j],nu), sqrt(beta) / P) T[0,1];
    x[2,j+1] ~ normal(x[2,j]*(1 - gamma) + x[1,j] - x[1,j+1], sqrt(gamma) / P) T[0,1-x[1,j+1]];
  }
  for(j in 1:N1) {
    y1[j] ~ lognormal(pow(b[1]*x[1,obs1[j]], varsigma[1]) + eta[1], sigma[1]);
  }
  for(j in 1:N2) {
    y2[j] ~ lognormal(pow(b[2]*x[1,obs2[j]], varsigma[2]) + eta[2], sigma[2]);
  }
  for(j in 1:N3) {
    y3[j] ~ lognormal(pow(b[3]*x[1,obs3[j]], varsigma[3]) + eta[3], sigma[3]);
  }
  for(j in 1:N4) {
    y4[j] ~ lognormal(pow(b[4]*x[1,obs4[j]], varsigma[4]) + eta[4], sigma[4]);
  }
}
