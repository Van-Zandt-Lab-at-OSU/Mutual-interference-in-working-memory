


#####################
#####################
# Model for the verbal version

model1 <- "
functions {
  row_vector compute_matrix(real C,real sigma,real md,real M){
      row_vector[3] result;
      real ai;
      real aj;
      
      ai = (1-C/2)^(md-1);
      aj = (C/2)*(1-C/2)^(md-2);
      result[1] = exp(ai/sigma)/(exp(ai/sigma)+(md-1)*exp(aj/sigma)+(M-md));
      result[2] = (md-1)*exp(aj/sigma)/(exp(ai/sigma)+(md-1)*exp(aj/sigma)+(M-md));          
      result[3] = (M-md)/(exp(ai/sigma)+(md-1)*exp(aj/sigma)+(M-md));  
      
      return(result);
  }
  real IG_lpdf(real t, real v,real b) {
      real density1;
      real density2;
      real density;
      density1 = log(b) - 0.5*(log(2)+log(pi())+3*log(t));
      density2 = -0.5*(v*t-b)^2/t;
      density = density1 + density2;
      return(density);
  }
  
  real mc_lpdf(vector rt,int md,row_vector prevector,matrix randommatrix,
       int[] loc,int[] accuracy,int[] ltype,real mu1,real mu2,real mu3,real tau,
       real v_u,real b,real v_r,real kappa,real sigma1,vector a_u,vector a_r,
       vector q,row_vector compute_1){
       
       
       matrix[3,3] cmat;
       real ai;
       real aj;
       row_vector[3] pre[md];
       real density;
       row_vector[3] q0;
       
       density = 0;
       
       
       cmat = randommatrix;
       cmat[1] = compute_1;
       cmat[2] = [1./8.,(md-1.)/8.,(8.-md)/8.];
       cmat[3] = [1./8.,(md-1.)/8.,(8.-md)/8.];  
       q0 = cmat[2];
       
       for(i in 1:md){
          pre[i] = prevector;
          pre[i,1] = 1;
       }
       
       for(i in 1:9){
             pre[loc[i]] = pre[loc[i]]*cmat;  
             if(rt[i]<=tau){
                density = density + log(a_u[1]) + lognormal_lpdf(rt[i] | mu1,1.) + binomial_lpmf(1 | 1,q[ltype[i]]);
             }
             else{
                density = density + log_sum_exp(log_sum_exp(log(a_u[1]) + lognormal_lpdf(rt[i] | mu1,1.) 
                               + binomial_lpmf(1 | 1,q[ltype[i]]),
                               log(a_u[2]) + IG_lpdf(rt[i]-tau | exp(v_u + kappa*pre[loc[i],ltype[i]]),b) 
                               + binomial_lpmf(1 | 1, pre[loc[i],ltype[i]])),
                               log(a_u[3]) + lognormal_lpdf(rt[i]-tau | mu3,1.)
                                + binomial_lpmf(1 | 1,q[ltype[i]]));             
             }
       }
       
       for(i in 1:md){
          pre[loc[9+i]] = pre[loc[9+i]]*cmat;
             if(rt[9+i]<=tau){
                density = density + log(a_r[1]) + lognormal_lpdf(rt[9+i] | mu2,sigma1) + binomial_lpmf(1 | 1,q[ltype[i]]);
             }
             else{
                density = density + log_sum_exp(log_sum_exp(log(a_r[1]) + lognormal_lpdf(rt[9+i] | mu2,sigma1)
                                + binomial_lpmf(1 | 1,q[ltype[9+i]]),
                               log(a_r[2]) + IG_lpdf(rt[9+i]-tau | exp(v_r + kappa*pre[loc[9+i],ltype[9+i]]),b)
                                + binomial_lpmf(1 | 1, pre[loc[9+i],ltype[9+i]])),
                               log(a_r[3]) + lognormal_lpdf(rt[9+i]-tau | mu3,1.)
                                + binomial_lpmf(1 | 1,q[ltype[9+i]]));             
             }
       }

      return density;
  }
}

data {
  int<lower=1> N;
  int<lower=1> Nsub;
  int<lower=1> Ngroup;
  int<lower=1> md[N];    // memory demand
  int<lower=0> id[N];    // subject belonging
  int<lower=0> gid[Nsub];
  int accuracy[N,14];       // accuracy of response
  int loc[N,14];
  int ltype[N,14];
  int task[N];
  int md_id[N];
  vector[14] rt[N];
  matrix[3,3] randommatrix;
  row_vector[3] prevector;
  real tau_mean[Nsub];
}


parameters {

  real C_0[2];
  real C_c[Ngroup,2];
  real C_raw[Nsub,2];
  real sd_C_raw[Ngroup,2];
  
  real v_u_0[2];
  real v_u_c[Ngroup,2];
  real v_u[Nsub,2];
  real sd_v_u_raw[Ngroup,2];

  real b_0[2];
  real b_c[Ngroup,2];
  real b_raw[Nsub,2];
  real sd_b_raw[Ngroup,2];
  
  simplex[3] a_u[Nsub,2];

  real v_r_0[2];
  real v_r_c[Ngroup,2];
  real v_r[Nsub,2];
  real sd_v_r_raw[Ngroup,2];

  real kappa_0[2];
  real kappa_c[Ngroup,2];
  real kappa[Nsub,2];
  real sd_kappa_raw[Ngroup,2];
  
  simplex[3] a_r[Nsub,2];
  
  real<lower=0,upper=0.5> tau[Nsub];
  simplex[3] q[Nsub,2];
  
  real sigma_0[2];
  real sigma_c[Ngroup,2];
  real sigma_raw[Nsub,2];
  real sd_sigma_raw[Ngroup,2];

  real<upper=-1> mu1[Nsub];
  real<upper=-0.5> mu2[Nsub];
  real<lower=1.5>  mu3[Nsub];
  
  real sdd_mu2_raw[Nsub];

}

transformed parameters {
  real<lower=0,upper=1> C[Nsub,2];
  real<lower=0> b[Nsub,2];
  
  real<lower=0> sd_C[Ngroup,2];
  real<lower=0> sd_v_u[Ngroup,2];
  real<lower=0> sd_b[Ngroup,2];

  real<lower=0> sd_v_r[Ngroup,2];
  real<lower=0> sd_kappa[Ngroup,2];
  
  real<lower=0,upper=1> sigma[Nsub,2];
  real<lower=0> sd_sigma[Ngroup,2];
  real<lower=0,upper=1> sdd_mu2[Nsub];
  
  for(n in 1:Nsub){
     for(l in 1:2){
       C[n,l] = inv_logit(C_raw[n,l]);
       b[n,l] = exp(b_raw[n,l]);         
       sigma[n,l] = inv_logit(sigma_raw[n,l]);
     }
     sdd_mu2[n] = inv_logit(sdd_mu2_raw[n]);
  }
  
  for(n in 1:Ngroup){
     for(l in 1:2){
       sd_C[n,l] = exp(sd_C_raw[n,l]);
       sd_v_u[n,l] = exp(sd_v_u_raw[n,l]);
       sd_b[n,l] = exp(sd_b_raw[n,l]);
         
       sd_v_r[n,l] = exp(sd_v_r_raw[n,l]);
       sd_kappa[n,l] = exp(sd_kappa_raw[n,l]);
       sd_sigma[n,l] = exp(sd_sigma_raw[n,l]);
     }
  }
}

model {
   int K;
   int M;
   row_vector[3] cmat[Nsub,2,2];
   real ai;
   real aj;

   for(l in 1:2){
     C_0[l] ~ normal(0.,0.2);
     v_u_0[l] ~ normal(0.,0.2);
     b_0[l] ~ normal(0.,0.2);
   
     v_r_0[l] ~ normal(0.,0.1);
     kappa_0[l] ~ normal(0.,0.2);
     sigma_0[l] ~ normal(0.,0.2);
   }
   
   for(n in 1:Ngroup){
      for(l in 1:2){
         C_c[n,l] ~ normal(C_0[l],0.2);
         v_u_c[n,l] ~ normal(v_u_0[l],0.2);
         kappa_c[n,l] ~ normal(kappa_0[l],0.2);
         
         v_r_c[n,l] ~ normal(v_r_0[l],0.1);
         b_c[n,l] ~ normal(b_0[l],0.2);
         
         sigma_c[n,l] ~ normal(sigma_0[l],0.2); 
         
         sd_C_raw[n,l] ~ normal(-1,0.2);
         sd_v_u_raw[n,l] ~ normal(-1,0.2);
         sd_b_raw[n,l] ~ normal(-1,0.2);
         sd_v_r_raw[n,l] ~ normal(-1,0.2);
         sd_kappa_raw[n,l] ~ normal(-1,0.2);
         sd_sigma_raw[n,l] ~ normal(-1.,0.2);  

      }
   }
   
   for(n in 1:Nsub){
      K = gid[n];
      for(l in 1:2){
         C_raw[n,l] ~ normal(C_c[K,l],sd_C[K,l]);
         v_u[n,l] ~ normal(v_u_c[K,l],sd_v_u[K,l]);
         b_raw[n,l] ~ normal(b_c[K,l],sd_b[K,l]);

         v_r[n,l] ~ normal(v_r_c[K,l],sd_v_r[K,l]);
         kappa[n,l] ~ normal(kappa_c[K,l],sd_kappa[K,l]);
         sigma_raw[n,l] ~ normal(sigma_c[K,l],sd_sigma[K,l]); 
         
      }
      
      sdd_mu2_raw[n] ~ normal(0.,1.);
      
      mu1[n] ~ normal(-2,0.05);
      mu2[n] ~ normal(-1,0.05);
      mu3[n] ~ normal(3,0.05);
      
      tau[n] ~ normal(tau_mean[n],0.001);   
   }
   
   for(n in 1:Nsub){
      cmat[n,1,1] = compute_matrix(C[n,1],sigma[n,1],2,8);
      cmat[n,1,2] = compute_matrix(C[n,2],sigma[n,2],2,8);
      cmat[n,2,1] = compute_matrix(C[n,1],sigma[n,1],4,8);
      cmat[n,2,2] = compute_matrix(C[n,2],sigma[n,2],4,8);
   }
   

   for(n in 1:N){
      M = task[n];
      K = id[n];
      target += mc_lpdf(rt[n] | md[n],prevector,randommatrix,loc[n],accuracy[n],ltype[n],
                   mu1[K],mu2[K],mu3[K],tau[K],v_u[K,M],b[K,M],v_r[K,M],kappa[K,M],
                   sdd_mu2[K],a_u[K,M],a_r[K,M],q[K,M],
                   cmat[K,md_id[n],M]);
   }
}
"

require(rstan)

load("updating_verbal.Rdata")


set.seed(117)  # seeds used are 111,117,177

dat$tau_mean = rep(0.1,dat$Nsub)
dat$tau_mean[142] = 0.15
dat$tau_mean[178] = 0.15
dat$md_id = ifelse(dat$md==2,1,2)

fit <- stan(model_code=model1,data=dat,chains=1,iter=2500,warmup=500,
            control = list(adapt_delta = 0.8,max_treedepth=10))

save(fit,file="test_1.Rdata")






