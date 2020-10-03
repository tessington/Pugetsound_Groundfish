data {
  int ndata;
  int nbeta;
  int ngamma;
  int nyear;
  int ntheta;
  int npresent;

  matrix[ndata, nbeta] x; // these are the intercept, depth and day effects
  matrix[ndata, nyear] t; // these are the years
  matrix[ndata, ngamma] b; // these are the  basin effects
  matrix[ndata, ntheta] r; // these are the region effects
  real y[ndata]; // catch rate
  int pa[ndata]; // vector of 0 and 1 for absent and present
  int present_index[npresent]; // which of the observations have non zero catch
  int basin_index[ntheta]; // says which basin each of the regions are in

}

parameters {
  vector[nbeta] beta_p; //the regression parameters depth and year on n
  vector[nbeta] beta_y; //the regression parameters depth and year on w
  vector[ngamma] gamma_p_raw; // basin effects
  vector[ngamma] gamma_y_raw; // basin effects
  vector[ntheta] theta_p_raw; // region effects
  vector[ntheta] theta_y_raw; // region effects
  vector[nyear] eps_p_raw;
  vector[nyear] eps_y_raw;
  real logitrho_p; // autocorrelation term
  real logitrho_y; // autocorrelation term
  real <lower = 0> sigma_psi_p; // standard deviation of year effect
  real <lower = 0> sigma_psi_y; // standard deviation of year effect
  real <lower = 0, upper = 100> shape; //the gamma shape parameter
  real <lower = 0> sigma_basin_p;
  real <lower = 0> sigma_basin_y;
  vector<lower =0>[ngamma] sigma_region_p;  // 
  vector<lower =0>[ngamma] sigma_region_y; //  variance of regional effects for by basins
}

transformed parameters {
real rho_p;
real rho_y;
vector[nyear] psi_p_raw;
vector[nyear] psi_y_raw;
vector[nyear] psi_p;
vector[nyear] psi_y;
vector[ndata] logit_phat;
vector[ndata] log_yhat;
vector[ndata] phat;
vector[ngamma] gamma_pa;
vector[ngamma] gamma_y;
vector[ntheta] theta_p;
vector[ntheta] theta_y;
vector[ndata] yhat;
real sigma_eps_p;
real sigma_eps_y;
vector[ndata] log_lik_pa;
vector[npresent] log_lik_y;

 rho_p = inv_logit(logitrho_p);
 rho_y = inv_logit(logitrho_y);
 sigma_eps_p = sigma_psi_p * (1 - rho_p)^0.5;
  sigma_eps_y = sigma_psi_y * (1 - rho_y)^0.5;
 
 psi_p_raw[1] = sigma_psi_p * eps_p_raw[1];
 psi_y_raw[1] = sigma_psi_y * eps_y_raw[1];

for (i in 2:nyear){
  psi_p_raw[i] = rho_p * psi_p_raw[i-1] + eps_p_raw[i] *sigma_eps_p;
  psi_y_raw[i] = rho_y * psi_y_raw[i-1] + eps_y_raw[i] *sigma_eps_y;
}
// center so that psi has mean of 0
for (i in 1:nyear) {
  psi_p[i] = psi_p_raw[i] - mean(psi_p_raw);
  psi_y[i] = psi_y_raw[i] - mean(psi_y_raw);
}
// noncentered basin effects
gamma_pa = gamma_p_raw* sigma_basin_p;
gamma_y = gamma_y_raw* sigma_basin_y;

// noncentered region effects
for (i in 1:ntheta) {
  theta_p[i] = theta_p_raw[i]* sigma_region_p[basin_index[i]];
   theta_y[i] = theta_y_raw[i]* sigma_region_y[basin_index[i]];
}

// get predicted values
logit_phat = x * beta_p + t * psi_p + b * gamma_pa + r * theta_p;
log_yhat = x * beta_y  + t * psi_y + b * gamma_y + r * theta_y;
phat = inv_logit(logit_phat);
yhat = exp(log_yhat);

for (i in 1:ndata) log_lik_pa[i] = bernoulli_lpmf(pa[i]|phat[i]);
for (i in 1:npresent) log_lik_y[i] = gamma_lpdf(y[present_index[i]]|shape, shape / yhat[present_index[i]]);

}
model {

sigma_basin_p ~ cauchy(0,2.5);
sigma_basin_y ~ cauchy(0,2.5);
logitrho_p ~ normal(0,1.75);
logitrho_y ~ normal(0,1.75);

// loop through betas to assign priors
for (i in 1:nbeta){
  beta_p[i] ~ cauchy(0,2.5);//prior for the slopes following Gelman 2008
  beta_y[i] ~ cauchy(0,2.5);//prior for the slopes following Gelman 2008
}

for (i in 1:nyear){
  eps_p_raw[i] ~ normal(0, 1);
  eps_y_raw[i] ~ normal(0, 1);
}

// loop through basin gammas and cauchy_scales
for (i in 1:ngamma){
  gamma_p_raw[i] ~ normal(0,1);// random effects for basin
  gamma_y_raw[i] ~ normal(0,1);// random effects for basin
  sigma_region_p[i] ~ cauchy(0,1.5); // scale variance by region
  sigma_region_y[i] ~ cauchy(0,1.5);; // scale variance by region


}

// loop through region gammas and cauchy_scales
for (i in 1:ntheta){
  theta_p_raw[i] ~ normal(0,1);
  theta_y_raw[i] ~ normal(0,1);
}


  // calculate likelihood of presence / absence

  for (i in 1:ndata) pa[i] ~ bernoulli(phat[i]);
  // calculate likelihood for non 0 catches
 // for (i in 1:npresent)  y[present_index[i]] ~ gamma(alpha_g[present_index[i]], beta_g[present_index[i]]);
for (i in 1:npresent)  y[present_index[i]] ~ gamma(shape, shape / yhat[present_index[i]]); 

}

generated quantities {

}
