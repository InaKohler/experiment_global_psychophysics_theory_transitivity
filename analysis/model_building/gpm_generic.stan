// gpm_generic.stan
//
// Global psychophysics model (Heller, 2021) for cross-modal magnitude
// production with three modalities:
//   1 = loudness (l), 2 = brightness (b), 3 = vibration strength (s)
//
// Reference parameters are estimated as potentially role-dependent:
//   rho_<f>to<g>   : reference of modality f as standard, target g
//   rho_<g>from<f> : reference of modality g as target, standard f
// Standard references are truncated below the lowest standard of their
// modality. alpha_b is fixed to 1 for identifiability. The weighting of the
// production factor is omega_1 * p^omega.
//
// Intensities (std, tgt) are in dB: dB SPL (loudness), dB Lambert
// (brightness), 40 * log10(d / 0.001578822 mm) (vibration).
// The variances are estimated depending on sigindex input.
// With onlyprior = 1 the model samples from the priors only.

functions {
  real db_lambert(real std_b) {return 10 * log10(std_b * pi() / (10 ^ -6));}
  real db_spl(real std_l) {return 20 * log10(std_l / (2 * (10 ^ -5)));}
  real db_inv_lambert(real db_b) {return (10 ^ (db_b / 10)) * (10 ^ -6) / pi();}
  real db_inv_spl(real db_l) {return (10 ^ (db_l / 20)) * 2 * (10 ^ -5);}
  real db_disp(real x_s) { return 40 * log10(x_s / 0.001578822); }
  real db_inv_disp(real db_s) { return 0.001578822 * 10^(db_s / 40); }
  real weigh_fun(real p, real omega_1, real omega) {
      return omega_1 * pow(p, omega);
  }
  real db_convert(real x, int modality) {
    if (modality == 1) return db_spl(x);
    else if (modality == 2) return db_lambert(x);
    else if (modality == 3) return db_disp(x);
    else reject("Unknown modality: ", modality, ". Please choose either 1 (loudness), 2 (brightness), 3 (vibration strength)");
  }
  real db_inverse(real x, int modality) {
    if (modality == 1) return db_inv_spl(x);
    else if (modality == 2) return db_inv_lambert(x);
    else if (modality == 3) return db_inv_disp(x);
    else reject("Unknown modality: ", modality, ". Please choose either 1 (loudness), 2 (brightness), 3 (vibration strength)");
  }
  
  real gpm_predict(real std, real alpha_std, real alpha_tgt,
                      real beta_std, real beta_tgt, 
                      real rho_std, real rho_tgt,
                      real omega_1, int p, real omega,
                      int std_modality, int tgt_modality) {
    return db_convert(
      pow(
        inv(alpha_tgt) *
          (weigh_fun(p, omega_1, omega) * alpha_std * 
             pow(db_inverse(std, std_modality), beta_std) -
           weigh_fun(p, omega_1, omega) * alpha_std * 
             pow(db_inverse(rho_std, std_modality), beta_std) +
           alpha_tgt * pow(db_inverse(rho_tgt, tgt_modality), beta_tgt)
          ),
        inv(beta_tgt))
    , tgt_modality);
  }
}

data {
  int<lower=1> ntotal;
  int<lower=1> ncond;
  int<lower=1> nsig;
  array[ncond] int<lower=1, upper=3> std_modality;
  array[ncond] int<lower=1, upper=3> tgt_modality;
  array[ncond] int<lower=1, upper=3> p;
  vector[ncond] std;                   // no <lower=1>: dB values can be negative
  array[ntotal] int<lower=1> idx;
  array[ntotal] int<lower=1> sigidx;
  array[ntotal] real tgt;              // no <lower=0>: dB values can be negative
  int<lower=0, upper=1> onlyprior;
  
  // lowest standards per modality — used to truncate rho_std parameters
  real lowestsoundstd;
  real lowestlightstd;
  real lowestvibrationstd;
  
  // priors
  real<lower=0> alpha_l_logmu;
  real<lower=0> alpha_l_logsigma;
  real<lower=0> alpha_s_logmu;
  real<lower=0> alpha_s_logsigma;
  real          beta_l_logmu;
  real<lower=0> beta_l_logsigma;
  real          beta_b_logmu;
  real<lower=0> beta_b_logsigma;
  real          beta_s_logmu;
  real<lower=0> beta_s_logsigma;
  real<lower=0> rho_ltol_mu;
  real<lower=0> rho_ltol_sigma;
  real<lower=0> rho_ltob_mu;
  real<lower=0> rho_ltob_sigma;
  real<lower=0> rho_ltos_mu;
  real<lower=0> rho_ltos_sigma;
  real<lower=0> rho_lfroml_mu;
  real<lower=0> rho_lfroml_sigma;
  real<lower=0> rho_lfromb_mu;
  real<lower=0> rho_lfromb_sigma;
  real<lower=0> rho_lfroms_mu;
  real<lower=0> rho_lfroms_sigma;
  real<lower=0> rho_btob_mu;
  real<lower=0> rho_btob_sigma;
  real<lower=0> rho_btol_mu;
  real<lower=0> rho_btol_sigma;
  real<lower=0> rho_btos_mu;
  real<lower=0> rho_btos_sigma;
  real<lower=0> rho_bfromb_mu;
  real<lower=0> rho_bfromb_sigma;
  real<lower=0> rho_bfroml_mu;
  real<lower=0> rho_bfroml_sigma;
  real<lower=0> rho_bfroms_mu;
  real<lower=0> rho_bfroms_sigma;
  real<lower=0> rho_stos_mu;
  real<lower=0> rho_stos_sigma;
  real<lower=0> rho_stol_mu;
  real<lower=0> rho_stol_sigma;
  real<lower=0> rho_stob_mu;
  real<lower=0> rho_stob_sigma;
  real<lower=0> rho_sfroms_mu;
  real<lower=0> rho_sfroms_sigma;
  real<lower=0> rho_sfroml_mu;
  real<lower=0> rho_sfroml_sigma;
  real<lower=0> rho_sfromb_mu;
  real<lower=0> rho_sfromb_sigma;
  real          omega_1_logmu;
  real<lower=0> omega_1_logsigma;
  real          omega_logmu;
  real<lower=0> omega_logsigma;
  real<lower=0> sig_mu;
  real<lower=0> sig_sigma;
}

transformed data {
  real alpha_b = 1.0;
}

parameters {
  real<lower=0> alpha_l;
  real<lower=0> alpha_s;
  real<lower=0> beta_l;
  real<lower=0> beta_b;
  real<lower=0> beta_s;
  real<lower=0> omega_1;
  real<lower=0> omega;
  vector<lower=0>[nsig] sig;
  // rho_std parameters truncated to be below the lowest standard of their modality
  real<upper=lowestsoundstd>     rho_ltob;
  real<upper=lowestsoundstd>     rho_ltos;
  real<upper=lowestsoundstd>     rho_ltol;
  real<upper=lowestlightstd>     rho_btol;
  real<upper=lowestlightstd>     rho_btos;
  real<upper=lowestlightstd>     rho_btob;
  real<upper=lowestvibrationstd> rho_stol;
  real<upper=lowestvibrationstd> rho_stob;
  real<upper=lowestvibrationstd> rho_stos;
  // rho_tgt parameters are unconstrained
  real rho_lfromb;
  real rho_lfroms;
  real rho_lfroml;
  real rho_bfroml;
  real rho_bfroms;
  real rho_bfromb;
  real rho_sfroml;
  real rho_sfromb;
  real rho_sfroms;
}
                      
transformed parameters {
  array[ncond] real mu; 
  for (i in 1:ncond) {
    real alpha_std_i;
    real alpha_tgt_i;
    real beta_std_i;
    real beta_tgt_i;
    real rho_std_i;
    real rho_tgt_i;
    if (std_modality[i] == 1){
      alpha_std_i = alpha_l;
      beta_std_i  = beta_l;
      if (tgt_modality[i] == 1){
        rho_std_i   = rho_ltol;
      }else if(tgt_modality[i] == 2){
        rho_std_i   = rho_ltob;
      }else if(tgt_modality[i] == 3){
        rho_std_i   = rho_ltos;
      }
    } else if (std_modality[i] == 2){
      alpha_std_i = alpha_b;
      beta_std_i  = beta_b;
      if (tgt_modality[i] == 1){
        rho_std_i   = rho_btol;
      }else if(tgt_modality[i] == 2){
        rho_std_i   = rho_btob;
      }else if(tgt_modality[i] == 3){
        rho_std_i   = rho_btos;
      }
    } else if (std_modality[i] == 3){
      alpha_std_i = alpha_s;
      beta_std_i  = beta_s;
      if (tgt_modality[i] == 1){
        rho_std_i   = rho_stol;
      }else if(tgt_modality[i] == 2){
        rho_std_i   = rho_stob;
      }else if(tgt_modality[i] == 3){
        rho_std_i   = rho_stos;
      }
    }
    if (tgt_modality[i] == 1){
      alpha_tgt_i = alpha_l;
      beta_tgt_i  = beta_l;
      if (std_modality[i] == 1){
        rho_tgt_i   = rho_lfroml;
      }else if(std_modality[i] == 2){
        rho_tgt_i   = rho_lfromb;
      }else if(std_modality[i] == 3){
        rho_tgt_i   = rho_lfroms;
      }
    } else if (tgt_modality[i] == 2){
      alpha_tgt_i = alpha_b;
      beta_tgt_i  = beta_b;
      if (std_modality[i] == 1){
        rho_tgt_i   = rho_bfroml;
      }else if(std_modality[i] == 2){
        rho_tgt_i   = rho_bfromb;
      }else if(std_modality[i] == 3){
        rho_tgt_i   = rho_bfroms;
      }
    } else if (tgt_modality[i] == 3){
      alpha_tgt_i = alpha_s;
      beta_tgt_i  = beta_s;
      if (std_modality[i] == 1){
        rho_tgt_i   = rho_sfroml;
      }else if(std_modality[i] == 2){
        rho_tgt_i   = rho_sfromb;
      }else if(std_modality[i] == 3){
        rho_tgt_i   = rho_sfroms;
      }
    }
    mu[i] = gpm_predict(std[i],
                        alpha_std_i, alpha_tgt_i,
                        beta_std_i,  beta_tgt_i,
                        rho_std_i,   rho_tgt_i,
                        omega_1, p[i], omega,
                        std_modality[i], tgt_modality[i]);
  }
}

model {
  target += lognormal_lpdf(alpha_l | alpha_l_logmu, alpha_l_logsigma);
  target += lognormal_lpdf(alpha_s | alpha_s_logmu, alpha_s_logsigma);
  target += lognormal_lpdf(beta_b | beta_b_logmu, beta_b_logsigma);
  target += lognormal_lpdf(beta_l | beta_l_logmu, beta_l_logsigma);
  target += lognormal_lpdf(beta_s | beta_s_logmu, beta_s_logsigma);
  // truncated rho_std priors: rho must be below the lowest standard of that modality
  target += normal_lpdf(rho_ltol | rho_ltol_mu, rho_ltol_sigma) -
              normal_lcdf(lowestsoundstd | rho_ltol_mu, rho_ltol_sigma);
  target += normal_lpdf(rho_ltob | rho_ltob_mu, rho_ltob_sigma) -
              normal_lcdf(lowestsoundstd | rho_ltob_mu, rho_ltob_sigma);
  target += normal_lpdf(rho_ltos | rho_ltos_mu, rho_ltos_sigma) -
              normal_lcdf(lowestsoundstd | rho_ltos_mu, rho_ltos_sigma);
  target += normal_lpdf(rho_btol | rho_btol_mu, rho_btol_sigma) -
              normal_lcdf(lowestlightstd | rho_btol_mu, rho_btol_sigma);
  target += normal_lpdf(rho_btob | rho_btob_mu, rho_btob_sigma) -
              normal_lcdf(lowestlightstd | rho_btob_mu, rho_btob_sigma);
  target += normal_lpdf(rho_btos | rho_btos_mu, rho_btos_sigma) -
              normal_lcdf(lowestlightstd | rho_btos_mu, rho_btos_sigma);
  target += normal_lpdf(rho_stol | rho_stol_mu, rho_stol_sigma) -
              normal_lcdf(lowestvibrationstd | rho_stol_mu, rho_stol_sigma);
  target += normal_lpdf(rho_stob | rho_stob_mu, rho_stob_sigma) -
              normal_lcdf(lowestvibrationstd | rho_stob_mu, rho_stob_sigma);
  target += normal_lpdf(rho_stos | rho_stos_mu, rho_stos_sigma) -
              normal_lcdf(lowestvibrationstd | rho_stos_mu, rho_stos_sigma);
  // unconstrained rho_tgt priors
  target += normal_lpdf(rho_lfroml | rho_lfroml_mu, rho_lfroml_sigma);
  target += normal_lpdf(rho_lfromb | rho_lfromb_mu, rho_lfromb_sigma);
  target += normal_lpdf(rho_lfroms | rho_lfroms_mu, rho_lfroms_sigma);
  target += normal_lpdf(rho_bfroml | rho_bfroml_mu, rho_bfroml_sigma);
  target += normal_lpdf(rho_bfromb | rho_bfromb_mu, rho_bfromb_sigma);
  target += normal_lpdf(rho_bfroms | rho_bfroms_mu, rho_bfroms_sigma);
  target += normal_lpdf(rho_sfroml | rho_sfroml_mu, rho_sfroml_sigma);
  target += normal_lpdf(rho_sfromb | rho_sfromb_mu, rho_sfromb_sigma);
  target += normal_lpdf(rho_sfroms | rho_sfroms_mu, rho_sfroms_sigma);
  target += lognormal_lpdf(omega_1 | omega_1_logmu, omega_1_logsigma);
  target += lognormal_lpdf(omega | omega_logmu, omega_logsigma);
  target += normal_lpdf(sig | sig_mu, sig_sigma) -
              normal_lccdf(0 | sig_mu, sig_sigma);
  
  if (!onlyprior) {
    target += normal_lpdf(tgt | to_vector(mu)[idx], to_vector(sig[sigidx]));
  }
}

generated quantities {
  array[ntotal] real tgt_pred = normal_rng(
    to_vector(mu)[idx],
    to_vector(sig[sigidx]));
}
