// generated with brms 2.17.0
functions {
}

data {
  int<lower=1> N;  // total number of observations
  vector[N] DOY;  // response variable
  vector[N] Year; //Explanatory variable
  int<lower=1> nYear; //Number of years
  }

transformed data {
}

parameters {
  //real Intercept;  // temporary intercept for centered predictors
  vector[nYear] b_int_raw; //standardized random effects on mean
  real<lower=0> sigma_int //Variance for random effect on mean
  
  real<lower=0> sigma;  // dispersion parameter
  real alpha;  // skewness parameter
}
transformed parameters {
  real lprior = 0;  // prior contributions to the log posterior
  lprior += normal_lpdf(Intercept | 110, 5);
  lprior += cauchy_lpdf(sigma | 0, 5);
  lprior += normal_lpdf(alpha | 0, 4);
  
  lprior += cauchy_lpdf(sigma_int | 0, 2);
  lprior += normal_lpdf(b1 | 0, sigma_int);

}
model {
  // likelihood including constants
  
    // initialize linear predictor term
    vector[N] mu = Intercept + rep_vector(0.0, N);
    // parameters used to transform the skew-normal distribution
    real delta;  // transformed alpha parameter
    real omega;  // scale parameter
    
    // use efficient skew-normal parameterization
    delta = alpha / sqrt(1 + alpha^2);
    omega = sigma / sqrt(1 - sqrt(2 / pi())^2 * delta^2);
    
    for (n in 1:N) {
      mu[n] = mu[n] - omega * delta * sqrt(2 / pi());
    }
    target += skew_normal_lpdf(DOY | mu, omega, alpha);

  // priors including constants
  target += lprior;
}

generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept;
}
