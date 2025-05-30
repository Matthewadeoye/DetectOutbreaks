functions{

  //Stationary distribution
   vector TPM(real G12, real G21){
    matrix[2, 2] m;
    m[1, 1] = 1 - G12;
    m[1, 2] = G12;
    m[2, 1] = G21;
    m[2, 2] = 1 - G21;

    matrix[2, 2] mT = transpose(m);

    complex_vector[2] E_values = eigenvalues(mT); vector[2] NE_values = get_real(E_values);
    complex_matrix[2, 2] E_vectors = eigenvectors(mT);

   int index;
   index = 1;
   if(abs(NE_values[2] - 1) < 1e-8) {
    index = 2;
  }
   complex_vector[2] stationary_distribution = E_vectors[, index];
   vector[2] Nstationary_distribution = get_real(stationary_distribution);
   Nstationary_distribution /= sum(Nstationary_distribution);
   return(Nstationary_distribution);
  }

  //colsums function
  array[] int colSums(array[,] int M){
    int ncol = dims(M)[2];
    array[ncol] int sums;

    for(i in 1:ncol){
      sums[i] = sum(M[,i]);
    }
    return(sums);
  }

//Forward filtering
array[] matrix Stanforwardfilter(array[,] int y, vector r, vector s, vector u, matrix gamma, matrix e_it, vector B, int Model, matrix z_it, matrix z_it2) {

  int ndept = dims(y)[1];
  int time = dims(y)[2];
  int nstate = rows(gamma);

  matrix[time, nstate] Alphas;
  matrix[time, nstate] Betas;
  matrix[time, nstate] alpha;
  matrix[time, nstate] beta;

  array[ndept] matrix[time, nstate] Allforwardprobs;

  //Model1 or Model2 or Model 4 or Model 5
  if(Model == 1 || Model == 2 || Model == 4 || Model == 5 || Model == 7){
  vector[nstate] init_density = TPM(gamma[1, 2], gamma[2, 1]);

  for (i in 1:ndept) {

    // Initialization of the first time step for each department

    alpha[1, 1] = log(init_density[1]) + poisson_lpmf(y[i, 1] | e_it[i, 1] * exp(r[1] + s[1] + u[i]));
    alpha[1, 2] = 0;
    beta[1, 1] = 0;
    beta[1, 2] = log(init_density[2]) + poisson_lpmf(y[i, 1] | e_it[i, 1] * exp(r[1] + s[1] + u[i] + (B[1] * z_it[i, 1])));

    Alphas[1, 1] = alpha[1, 1];
    Alphas[1, 2] = beta[1, 2];

    // Dynamic programming loop for the remaining time steps
    for (t in 2:time) {

      alpha[t, 1] = Alphas[t-1, 1] + log(gamma[1, 1]) + poisson_lpmf(y[i, t] | e_it[i, t] * exp(r[t] + s[t] + u[i]));
      alpha[t, 2] = Alphas[t-1, 2] + log(gamma[2, 1]) + poisson_lpmf(y[i, t] | e_it[i, t] * exp(r[t] + s[t] + u[i]));
      beta[t, 1] =  Alphas[t-1, 1] + log(gamma[1, 2]) + poisson_lpmf(y[i, t] | e_it[i, t] * exp(r[t] + s[t] + u[i] + (B[1] * z_it[i, t])));
      beta[t, 2] =  Alphas[t-1, 2] + log(gamma[2, 2]) + poisson_lpmf(y[i, t] | e_it[i, t] * exp(r[t] + s[t] + u[i] + (B[1] * z_it[i, t])));
      Alphas[t, 1] = log_sum_exp(alpha[t,]);
      Alphas[t, 2] = log_sum_exp(beta[t,]);
    }

    Allforwardprobs[i] = Alphas;
  }
  return(Allforwardprobs);
 }
 else //Model3 or Model 6
  if(Model == 3 || Model == 6){
    vector[nstate] init_density = TPM(gamma[1, 2], gamma[2, 1]);

  for (i in 1:ndept) {

    // Initialization of the first time step for each department

    alpha[1, 1] = log(init_density[1]) + poisson_lpmf(y[i, 1] | e_it[i, 1] * exp(r[1] + s[1] + u[i]));
    alpha[1, 2] = 0;
    beta[1, 1] = 0;
    beta[1, 2] = log(init_density[2]) + poisson_lpmf(y[i, 1] | e_it[i, 1] * exp(r[1] + s[1] + u[i] + (B[1] * z_it[i, 1]) + (B[2] * z_it2[i, 1])));

    Alphas[1, 1] = alpha[1, 1];
    Alphas[1, 2] = beta[1, 2];

    // Dynamic programming loop for the remaining time steps
    for (t in 2:time) {

      alpha[t, 1] = Alphas[t-1, 1] + log(gamma[1, 1]) + poisson_lpmf(y[i, t] | e_it[i, t] * exp(r[t] + s[t] + u[i]));
      alpha[t, 2] = Alphas[t-1, 2] + log(gamma[2, 1]) + poisson_lpmf(y[i, t] | e_it[i, t] * exp(r[t] + s[t] + u[i]));
      beta[t, 1] =  Alphas[t-1, 1] + log(gamma[1, 2]) + poisson_lpmf(y[i, t] | e_it[i, t] * exp(r[t] + s[t] + u[i] + (B[1] * z_it[i, t]) + (B[2] * z_it2[i, t])));
      beta[t, 2] =  Alphas[t-1, 2] + log(gamma[2, 2]) + poisson_lpmf(y[i, t] | e_it[i, t] * exp(r[t] + s[t] + u[i] + (B[1] * z_it[i, t]) + (B[2] * z_it2[i, t])));
      Alphas[t, 1] = log_sum_exp(alpha[t,]);
      Alphas[t, 2] = log_sum_exp(beta[t,]);
    }

    Allforwardprobs[i] = Alphas;
  }
  return(Allforwardprobs);
 }

 //if Model specification outside expected range
 for(i in 1: ndept){
   matrix[time, nstate] dummy;
   for(t in 1:time){
     dummy[i, t] = 0;
   }
   Allforwardprobs[i] = dummy;
 }
 return Allforwardprobs;
}

//Backward sweep
array[] matrix Stanbackwardsweep(array[,] int y, vector r, vector s, vector u, matrix gamma, matrix e_it, vector B, int Model, matrix z_it, matrix z_it2) {

  int ndept = dims(y)[1];
  int time = dims(y)[2];
  int nstate = rows(gamma);

  matrix[time, nstate] Alphas;
  matrix[time, nstate] alpha;
  matrix[time, nstate] beta;

  array[ndept] matrix[time, nstate] Allbackwardprobs;

  array[time-1] int ind;
  for (i in 1:(time-1)) {
    ind[i] = i;
}

  //Model1 or Model2 or Model 4 or Model 5
  if(Model == 1 || Model == 2 || Model == 4 || Model == 5 || Model == 7){

  for (i in 1:ndept) {

    // Initialization of the last time step for each department

    alpha[time, 1] = 0;
    alpha[time, 2] = 0;
    beta[time, 1] = 0;
    beta[time, 2] = 0;

    Alphas[time, 1] = alpha[time, 1];
    Alphas[time, 2] = beta[time, 2];

    // Dynamic programming loop for the remaining time steps
    for (t in sort_indices_desc(ind)) {

      alpha[t, 1] = Alphas[t+1, 1] + log(gamma[1, 1]) + poisson_lpmf(y[i, t+1] | e_it[i, t+1] * exp(r[t+1] + s[t+1] + u[i]));
      alpha[t, 2] = Alphas[t+1, 2] + log(gamma[1, 2]) + poisson_lpmf(y[i, t+1] | e_it[i, t+1] * exp(r[t+1] + s[t+1] + u[i]));
      beta[t, 1] =  Alphas[t+1, 1] + log(gamma[2, 1]) + poisson_lpmf(y[i, t+1] | e_it[i, t+1] * exp(r[t+1] + s[t+1] + u[i] + (B[1] * z_it[i, t+1])));
      beta[t, 2] =  Alphas[t+1, 2] + log(gamma[2, 2]) + poisson_lpmf(y[i, t+1] | e_it[i, t+1] * exp(r[t+1] + s[t+1] + u[i] + (B[1] * z_it[i, t+1])));
      Alphas[t, 1] = log_sum_exp(alpha[t, ]);
      Alphas[t, 2] = log_sum_exp(beta[t, ]);
    }
    Allbackwardprobs[i] = Alphas;
  }
  return(Allbackwardprobs);
 }
 //Model 3 or Model 6
 else if(Model == 3 || Model == 6){
  for (i in 1:ndept) {

    // Initialization of the last time step for each department

    alpha[time, 1] = 0;
    alpha[time, 2] = 0;
    beta[time, 1] = 0;
    beta[time, 2] = 0;

    Alphas[time, 1] = alpha[time, 1];
    Alphas[time, 2] = beta[time, 2];

    // Dynamic programming loop for the remaining time steps
    for (t in sort_indices_desc(ind)) {

      alpha[t, 1] = Alphas[t+1, 1] + log(gamma[1, 1]) + poisson_lpmf(y[i, t+1] | e_it[i, t+1] * exp(r[t+1] + s[t+1] + u[i]));
      alpha[t, 2] = Alphas[t+1, 2] + log(gamma[1, 2]) + poisson_lpmf(y[i, t+1] | e_it[i, t+1] * exp(r[t+1] + s[t+1] + u[i]));
      beta[t, 1] =  Alphas[t+1, 1] + log(gamma[2, 1]) + poisson_lpmf(y[i, t+1] | e_it[i, t+1] * exp(r[t+1] + s[t+1] + u[i] + (B[1] * z_it[i, t+1]) + (B[2] * z_it2[i, t+1])));
      beta[t, 2] =  Alphas[t+1, 2] + log(gamma[2, 2]) + poisson_lpmf(y[i, t+1] | e_it[i, t+1] * exp(r[t+1] + s[t+1] + u[i] + (B[1] * z_it[i, t+1]) + (B[2] * z_it2[i, t+1])));
      Alphas[t, 1] = log_sum_exp(alpha[t, ]);
      Alphas[t, 2] = log_sum_exp(beta[t, ]);
    }
    Allbackwardprobs[i] = Alphas;
  }
  return(Allbackwardprobs);
 }

 //if Model specification outside expected range
 for(i in 1: ndept){
   matrix[time, nstate] dummy;
   for(t in 1:time){
     dummy[i, t] = 0;
   }
   Allbackwardprobs[i] = dummy;
 }
 return Allbackwardprobs;
}

//Local decoding
matrix StanDecoding(array[,] int y, vector r, vector s, vector u, matrix gamma, matrix e_it, vector B, int Model, matrix z_it, matrix z_it2) {
  int ndept = dims(y)[1];
  int time = dims(y)[2];
  int nstate = rows(gamma);

if(Model == 1 || Model == 2 || Model == 3 || Model == 4 || Model == 5 || Model == 6 || Model == 7){
  array[ndept] matrix[time, nstate] Allforwardprobs;
  Allforwardprobs = Stanforwardfilter(y, r, s, u, gamma, e_it, B, Model, z_it, z_it2);
  array[ndept] matrix[time, nstate] Allbackwardprobs;
  Allbackwardprobs = Stanbackwardsweep(y, r, s, u, gamma, e_it, B, Model, z_it, z_it2);
  matrix[ndept, time] Res;
  real P1;
  real P2;
  matrix[time, nstate] P1P2;

  for(i in 1:ndept){
    for(j in 1:time){
      P1 = Allforwardprobs[i][j,1] + Allbackwardprobs[i][j,1];
      P2 = Allforwardprobs[i][j,2] + Allbackwardprobs[i][j,2];
      P1P2[j, 1] = P1;
      P1P2[j, 2] = P2;
      Res[i,j]= exp(P2 - log_sum_exp(P1P2[j, ]));

    }
  }
  return(Res);
 }
 //If Model specification is zero or out of range.
  matrix[ndept, time] dummy;
 for(i in 1:ndept){
   for(t in 1:time){
     dummy[i, t] = 0;
   }
 }
  return dummy;
}
  // Intrinsic GMRF density
  real IGMRF1_lpdf(vector uconstrained, real kappa_u, matrix R, int rankdef) {
    int n = rows(R);
    return (((n - rankdef) / 2.0) * (log(kappa_u) - log(2.0 * pi())) - (kappa_u / 2.0) * quad_form(R, uconstrained));
  }

  // Random walk density
  real randomwalk2_lpdf(vector r, real kappa_r) {
    int time = dims(r)[1];
    real res = 0;
    for (i in 3:time) {
      res += (r[i-2] - (2 * r[i-1]) + r[i])^2;
    }
    return (((time - 2) / 2.0) * log(kappa_r) - (kappa_r / 2.0) * res);
  }

  // Seasonal components' density
    real seasonalComp_lpdf(vector s, real kappa_s) {
    int time = dims(s)[1];
    real res = 0;
    for (i in 12:time) {
    res += (sum(s[(i-11):(i-0)]))^2;
  }
    return (((time - 11) / 2.0) * log(kappa_s) - (kappa_s / 2.0) * res);
}

  //Gamma matrix function

    matrix G(real G12, real G21) {
    matrix[2, 2] m;
    m[1, 1] = 1 - G12;
    m[1, 2] = G12;
    m[2, 1] = G21;
    m[2, 2] = 1 - G21;
    return m;
  }

  //loglikelihood via forward filtering

  real Stan_Loglikelihood(array[,] int y, vector r, vector s, vector u, matrix gamma, matrix e_it, vector B, int Model, matrix z_it, matrix z_it2) {

  int ndept = dims(y)[1];
  int time = dims(y)[2];
  int nstate = rows(gamma);

  matrix[ndept, nstate] logforwardprobs;
  matrix[time, nstate] Alphas;
  matrix[time, nstate] Betas;
  matrix[time, nstate] alpha;
  matrix[time, nstate] beta;

  if(Model == 0){
  //Model0
    matrix[ndept, time] allLoglikelihood;

  for (i in 1:ndept) {
    for (t in 1:time) {
      if(y[i, t] == -1){
        allLoglikelihood[i, t] = 0;
      }else{
         allLoglikelihood[i, t] = poisson_lpmf(y[i, t] | e_it[i, t] * exp(r[t] + s[t] + u[i]));
      }
   }
}
      real fullLogLikelihood = sum(allLoglikelihood);
      return fullLogLikelihood;
  }

 //Model1 or Model2 or Model4 or Model5 or Model7
  else if(Model == 1 || Model == 2 || Model == 4 || Model == 5 || Model == 7){
    vector[nstate] init_density = TPM(gamma[1, 2], gamma[2, 1]);

  for (i in 1:ndept) {

  // Initialization of the first time step for each department
    if(y[i, 1] == -1){
    alpha[1, 1] = log(init_density[1]);
    alpha[1, 2] = 0;
    beta[1, 1] = 0;
    beta[1, 2] = log(init_density[2]);
  }else{
    alpha[1, 1] = log(init_density[1]) + poisson_lpmf(y[i, 1] | e_it[i, 1] * exp(r[1] + s[1] + u[i]));
    alpha[1, 2] = 0;
    beta[1, 1] = 0;
    beta[1, 2] = log(init_density[2]) + poisson_lpmf(y[i, 1] | e_it[i, 1] * exp(r[1] + s[1] + u[i] + (B[1] * z_it[i, 1])));
  }
    Alphas[1, 1] = alpha[1, 1];
    Alphas[1, 2] = beta[1, 2];

  // Dynamic programming loop for the remaining time steps
      for (t in 2:time) {
        if(y[i, t] == -1){
         alpha[t, 1] = Alphas[t-1, 1] + log(gamma[1, 1]);
         alpha[t, 2] = Alphas[t-1, 2] + log(gamma[2, 1]);
         beta[t, 1] =  Alphas[t-1, 1] + log(gamma[1, 2]);
         beta[t, 2] =  Alphas[t-1, 2] + log(gamma[2, 2]);
         Alphas[t, 1] = log_sum_exp(alpha[t,]);
         Alphas[t, 2] = log_sum_exp(beta[t,]);
        }else{
         alpha[t, 1] = Alphas[t-1, 1] + log(gamma[1, 1]) + poisson_lpmf(y[i, t] | e_it[i, t] * exp(r[t] + s[t] + u[i]));
         alpha[t, 2] = Alphas[t-1, 2] + log(gamma[2, 1]) + poisson_lpmf(y[i, t] | e_it[i, t] * exp(r[t] + s[t] + u[i]));
         beta[t, 1] =  Alphas[t-1, 1] + log(gamma[1, 2]) + poisson_lpmf(y[i, t] | e_it[i, t] * exp(r[t] + s[t] + u[i] + (B[1] * z_it[i, t])));
         beta[t, 2] =  Alphas[t-1, 2] + log(gamma[2, 2]) + poisson_lpmf(y[i, t] | e_it[i, t] * exp(r[t] + s[t] + u[i] + (B[1] * z_it[i, t])));
         Alphas[t, 1] = log_sum_exp(alpha[t,]);
         Alphas[t, 2] = log_sum_exp(beta[t,]);
        }
}

      // save the final forward probabilities

      logforwardprobs[i,1] = log_sum_exp(alpha[time,]);
      logforwardprobs[i,2] = log_sum_exp(beta[time,]);

    }

      vector[ndept] rowlogsumexp;
      for (i in 1:ndept) {
      rowlogsumexp[i] = log_sum_exp(logforwardprobs[i,]);
    }
      real fullLogLikelihood = sum(rowlogsumexp);
      return fullLogLikelihood;
  }

  //Model3 or Model 6
  else if(Model == 3 || Model == 6){
      vector[nstate] init_density = TPM(gamma[1, 2], gamma[2, 1]);
  for (i in 1:ndept) {

  // Initialization of the first time step for each department
    if(y[i, 1] == -1){
    alpha[1, 1] = log(init_density[1]);
    alpha[1, 2] = 0;
    beta[1, 1] = 0;
    beta[1, 2] = log(init_density[2]);
  }else{
    alpha[1, 1] = log(init_density[1]) + poisson_lpmf(y[i, 1] | e_it[i, 1] * exp(r[1] + s[1] + u[i]));
    alpha[1, 2] = 0;
    beta[1, 1] = 0;
    beta[1, 2] = log(init_density[2]) + poisson_lpmf(y[i, 1] | e_it[i, 1] * exp(r[1] + s[1] + u[i] + (B[1] * z_it[i, 1]) + (B[2] * z_it2[i, 1])));
}
    Alphas[1, 1] = alpha[1, 1];
    Alphas[1, 2] = beta[1, 2];

  // Dynamic programming loop for the remaining time steps
      for (t in 2:time) {
         if(y[i, t] == -1){
         alpha[t, 1] = Alphas[t-1, 1] + log(gamma[1, 1]);
         alpha[t, 2] = Alphas[t-1, 2] + log(gamma[2, 1]);
         beta[t, 1] =  Alphas[t-1, 1] + log(gamma[1, 2]);
         beta[t, 2] =  Alphas[t-1, 2] + log(gamma[2, 2]);
         Alphas[t, 1] = log_sum_exp(alpha[t,]);
         Alphas[t, 2] = log_sum_exp(beta[t,]);
        }else{
         alpha[t, 1] = Alphas[t-1, 1] + log(gamma[1, 1]) + poisson_lpmf(y[i, t] | e_it[i, t] * exp(r[t] + s[t] + u[i]));
         alpha[t, 2] = Alphas[t-1, 2] + log(gamma[2, 1]) + poisson_lpmf(y[i, t] | e_it[i, t] * exp(r[t] + s[t] + u[i]));
         beta[t, 1] =  Alphas[t-1, 1] + log(gamma[1, 2]) + poisson_lpmf(y[i, t] | e_it[i, t] * exp(r[t] + s[t] + u[i] + (B[1] * z_it[i, t]) + (B[2] * z_it2[i, t])));
         beta[t, 2] =  Alphas[t-1, 2] + log(gamma[2, 2]) + poisson_lpmf(y[i, t] | e_it[i, t] * exp(r[t] + s[t] + u[i] + (B[1] * z_it[i, t]) + (B[2] * z_it2[i, t])));
         Alphas[t, 1] = log_sum_exp(alpha[t,]);
         Alphas[t, 2] = log_sum_exp(beta[t,]);
        }
}
      // save the final forward probabilities

      logforwardprobs[i,1] = log_sum_exp(alpha[time,]);
      logforwardprobs[i,2] = log_sum_exp(beta[time,]);

    }

      vector[ndept] rowlogsumexp;
      for (i in 1:ndept) {
      rowlogsumexp[i] = log_sum_exp(logforwardprobs[i,]);
    }
      real fullLogLikelihood = sum(rowlogsumexp);
      return fullLogLikelihood;
    }
   return 0;
  }
}

data {
  int<lower=1> ndept;                 // Number of departments
  int<lower=1> time;                  // Time
  int<lower=1> nstate;                // Number of states
  int<lower=1> rankdef;               // Rank deficiency of structure matrix (R)
  array[ndept, time] int y;           // data matrix
  matrix[ndept, time] e_it;           // initial Susceptibles
  matrix[ndept, ndept] R;             // Structure / Precision matrix (IGMRF1/IGMRF2)
  int<lower=0, upper=7> Model;        // Model's functional form
  int<lower=0, upper=2> npar;         // Number of autoregressive parameters
  matrix[ndept, time] z_it;           //designMatrix1
  matrix[ndept, time] z_it2;          //designMatrix2
}

parameters {
  real<lower=0, upper=1> G12;            // transition to hyperendemic
  real<lower=0, upper=1> G21;            // transition to endemic
  real<lower=0> kappa_u;                 // spatial precision parameter
  real<lower=0> kappa_r;                 // trend precision parameter
  real<lower=0> kappa_s;                 // seasonal precision parameter
  vector[ndept-1] u;                     // Spatial components
  vector[time] r;                        // Trend components
  vector[time] s;                        // Seasonal components
  vector<lower=0>[npar] B;               // autoregressive parameters
}

transformed parameters {
  real sumC = sum(u[1:(ndept-1)]);
  vector[ndept] uconstrained;
  uconstrained = append_row(u[1:(ndept-1)], -sumC);
}

model {
  // Priors
  G12 ~ beta(2, 2);
  G21 ~ beta(2, 2);
  kappa_u ~ gamma(1, 0.01);
  kappa_r ~ gamma(1, 0.0001);
  kappa_s ~ gamma(1, 0.001);

  if(Model==1 || Model==2 || Model==4 || Model==5 || Model==7){
    B[1] ~ gamma(2, 2);
  }

  if(Model==3 || Model==6){
    B[1] ~ gamma(2, 2);
    B[2] ~ gamma(2, 2);
  }

  uconstrained ~ IGMRF1(kappa_u, R, rankdef);
  r ~ randomwalk2(kappa_r);
  s ~ seasonalComp(kappa_s);

  // Likelihood
  target += Stan_Loglikelihood(y, r, s, uconstrained, G(G12, G21), e_it, B, Model, z_it, z_it2);
}

generated quantities{
  real log_lik = Stan_Loglikelihood(y, r, s, uconstrained, G(G12, G21), e_it, B, Model, z_it, z_it2);
  real state1_stationary_dist = TPM(G12, G21)[2];
//  array[ndept, time] int predY_it;

//if(Model == 0){
//for(i in 1:ndept){
//  for(t in 1:time){
//    real expected_it = e_it[i, t] * exp(r[t] + s[t] + uconstrained[i]);
//    predY_it[i, t] = poisson_rng(expected_it);
//  }
// }
//}
//else if(Model == 1 || Model == 2 || Model == 4 || Model == 5 || Model == 7){
//  matrix[ndept, time] x_it = StanDecoding(y, r, s, uconstrained, G(G12, G21), e_it, B, Model, z_it, z_it2);
//
//for(i in 1:ndept){
//  for(t in 1:time){
//    int indicator = binomial_rng(1, x_it[i, t]);
//    real expected_it = e_it[i, t] * exp(r[t] + s[t] + uconstrained[i] + indicator * z_it[i, t] * B[1]);
//    predY_it[i, t] = poisson_rng(expected_it);
//  }
// }
//}
//else if(Model == 3 || Model == 6){
//  matrix[ndept, time] x_it = StanDecoding(y, r, s, uconstrained, G(G12, G21), e_it, B, Model, z_it, z_it2);

//for(i in 1:ndept){
//  for(t in 1:time){
//    int indicator = binomial_rng(1, x_it[i, t]);
//    real expected_it = e_it[i, t] * exp(r[t] + s[t] + uconstrained[i] + indicator * z_it[i, t] * B[1] + indicator *  z_it2[i, t] * B[2]);
//    predY_it[i, t] = poisson_rng(expected_it);
//  }
// }
//}

//Overall predicted cases
//  array[time] int colsumY = colSums(predY_it);
}
