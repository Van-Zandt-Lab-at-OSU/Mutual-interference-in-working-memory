


# Hierarchical Bayesian model for binned data

model1 <- "
functions {
  row_vector compute_matrix(real C,real sigma,real md,real M,real r,real up_t){
      row_vector[3] result;
      real ai;
      real aj;
      
      ai = (1-C/2)^(md-1)*(1-exp(-r*up_t));
      aj = (C/2)*(1-C/2)^(md-2)*(1-exp(-r*up_t));
      result[1] = exp(ai/sigma)/(exp(ai/sigma)+(md-1)*exp(aj/sigma)+(9-md));
      result[2] = (md-1)*exp(aj/sigma)/(exp(ai/sigma)+(md-1)*exp(aj/sigma)+(9-md));          
      result[3] = (9-md)/(exp(ai/sigma)+(md-1)*exp(aj/sigma)+(9-md));  
      
      return(result);
  }
  real IG_lpdf(real t, real v,real b) {
      real density1;
      real density2;
      real density;
      density1 = log(b) - 0.5*(log(2) +log(pi())+3*log(t));
      density2 = -0.5*(v*t-b)^2/t;
      density = density1 + density2;
      return(density);
  }
  real mc_lpdf(vector rt,int md,int mm,row_vector prevector,matrix randommatrix,
       int[] accuracy,int[] ltype,real mu1,real mu2,real tau,
       real C,real sigma,real v,real kappa,real b,real sigma1,
       vector a,vector q,real r,row_vector compute_1){
       
       
       matrix[3,3] cmat;
       real ai;
       real aj;
       row_vector[3] pre;
       real density;
       
       density = 0;
       
       
       cmat = randommatrix;
       cmat[1] = compute_1;
       cmat[2] = [1./9.,(md-1.)/9.,(9.-md)/9.];
       cmat[3] = [1./9.,(md-1.)/9.,(9.-md)/9.];    
       
        pre = prevector;
        pre[1] = 1.;
        
       
       for(i in 1:mm){
          pre = pre*cmat; 
       }
       
       for(i in 1:md){
             if(rt[i]<=tau){
                density = density + log(a[1]) + lognormal_lpdf(rt[i] | mu1,sigma1) + binomial_lpmf(1 | 1,q[ltype[i]]);
             }
             else{
                density = density + log_sum_exp(log_sum_exp(log(a[1]) + lognormal_lpdf(rt[i] | mu1,sigma1)
                                + binomial_lpmf(1 | 1,q[ltype[i]]),
                               log(a[2]) + IG_lpdf(rt[i]-tau | exp(v+pre[ltype[i]]*kappa),b)
                                + binomial_lpmf(1 | 1, pre[ltype[i]])),
                               log(a[3]) + lognormal_lpdf(rt[i]-tau | mu2,1.)
                                + binomial_lpmf(1 | 1,q[ltype[i]]));             
             }
       }

      return density;
  }
}

data {
  int<lower=1> N;
  int<lower=1> Nsub;      // participant number
  int<lower=1> Ngroup;    // group number
  int N_upt;
  real up_t_list[N_upt];
  int up_t[N];
  int gid[Nsub];
  int id[N];
  int<lower=1> md[N];    // memory demand
  int mm[N];
  int accuracy[N,6];       // accuracy of response
  int ltype[N,6];
  vector[6] rt[N];
  matrix[3,3] randommatrix;
  row_vector[3] prevector;
}


parameters {

  real C_0;               
  real C_c[Ngroup];
  real C_raw[Nsub];
  real sd_C_raw[Ngroup];
  
  real r_0;               
  real r_c[Ngroup];
  real r_raw[Nsub];
  real sd_r_raw[Ngroup];
  
  real sigma_0;               
  real sigma_c[Ngroup];
  real sigma_raw[Nsub];
  real sd_sigma_raw[Ngroup];
  
  real v_0;
  real v_c[Ngroup];
  real v[Nsub];
  real sd_v_raw[Ngroup];
  
  real kappa_0;
  real kappa_c[Ngroup];
  real kappa[Nsub];
  real sd_kappa_raw[Ngroup];
  
  real b_0;
  real b_c[Ngroup];
  real b_raw[Nsub];
  real sd_b_raw[Ngroup];

  real<upper=-0.5> mu1[Nsub];
  real<lower=3> mu2[Nsub];
  
  real sdd_mu1_raw[Nsub];
  
  real<lower=0.1,upper=0.5> tau[Nsub];
  simplex[3] a[Nsub];
  simplex[3] q[Nsub];
}

transformed parameters {
  real<lower=0> b[Nsub];
  real<lower=0> sigma[Nsub];
  real<lower=0,upper=1> C[Nsub];
  real<lower=0> r[Nsub];
  real<lower=0,upper=1> sdd_mu1[Nsub];
  
  real<lower=0> sd_C[Ngroup];
  real<lower=0> sd_sigma[Ngroup];
  real<lower=0> sd_v[Ngroup];
  real<lower=0> sd_kappa[Ngroup];
  real<lower=0> sd_b[Ngroup];
  real<lower=0> sd_r[Ngroup];
  
  for(i in 1:Nsub){
     b[i] = exp(b_raw[i]);
     sigma[i] = exp(sigma_raw[i]);
     C[i] = inv_logit(C_raw[i]);
     r[i] = exp(r_raw[i]);
     sdd_mu1[i] = inv_logit(sdd_mu1_raw[i]);
  }
  
  for(i in 1:Ngroup){
    sd_C[i] = exp(sd_C_raw[i]);
    sd_sigma[i] = exp(sd_sigma_raw[i]); 
    sd_v[i] = exp(sd_v_raw[i]);
    sd_kappa[i] = exp(sd_kappa_raw[i]);
    sd_b[i] = exp(sd_b_raw[i]);
    sd_r[i] = exp(sd_r_raw[i]);    
  }
}

model {

  int K;
  row_vector[3] cmat[Nsub,N_upt,4];

  v_0 ~ normal(0.,0.2);
  b_0 ~ normal(0.,0.2);
  kappa_0 ~ normal(0.,0.2);
  sigma_0 ~ normal(-1,0.2);
  C_0 ~ normal(0.,0.2);
  r_0 ~ normal(0.,0.2);
     
  for(i in 1:Ngroup){
    v_c[i] ~ normal(v_0,0.2);
    b_c[i] ~ normal(b_0,0.2);
    kappa_c[i] ~ normal(kappa_0,0.2);
    sigma_c[i] ~ normal(sigma_0,0.2);
    C_c[i] ~ normal(C_0,0.2);
    r_c[i] ~ normal(r_0,0.2);
    
    sd_v_raw[i] ~ normal(-1,0.2);
    sd_b_raw[i] ~ normal(-1,0.2);
    sd_kappa_raw[i] ~ normal(-1,0.2);
    sd_sigma_raw[i] ~ normal(-1,0.2);
    sd_C_raw[i] ~ normal(-1,0.2);
    sd_r_raw[i] ~ normal(-1,0.2);  
  }

  for(i in 1:Nsub){
    K = gid[i];
    v[i] ~ normal(v_c[K],sd_v[K]);
    b_raw[i] ~ normal(b_c[K],sd_b[K]);
    kappa[i] ~ normal(kappa_c[K],sd_kappa[K]);
    sigma_raw[i] ~ normal(sigma_c[K],sd_sigma[K]);
    C_raw[i] ~ normal(C_c[K],sd_C[K]);
    r_raw[i] ~ normal(r_c[K],sd_r[K]); 
    
    mu1[i] ~ normal(-1,0.05);
    mu2[i] ~ normal(3,0.05);
    sdd_mu1_raw[i] ~ normal(0.,1.);
    tau[i] ~ normal(0.15,0.001);
  }

  for(i in 1:Nsub){
     for(n in 1:N_upt){
        for(j in 1:4){
            cmat[i,n,j] = compute_matrix(C[i],sigma[i],j,9,r[i],up_t_list[n]);
        }
     }
  }

  for(n in 1:N){
     K = id[n];
     target += mc_lpdf(rt[n] | md[n],mm[n],prevector,randommatrix,accuracy[n],ltype[n],
                      mu1[K],mu2[K],tau[K],C[K],sigma[K],v[K],kappa[K],b[K],sdd_mu1[K],
                      a[K],q[K],r[K],cmat[K,up_t[n],md[n]]);

  }
}
"


require(rstan)

set.seed(117)  # seeds used are 117, 177 and 111
  
load("dat_interference.Rdata")

fit <- stan(model_code=model1,data=sub,chains=1,iter=2500,warmup=500,
              control = list(adapt_delta = 0.8,max_treedepth=10))

save(fit,file="result_1.Rdata")
  






























