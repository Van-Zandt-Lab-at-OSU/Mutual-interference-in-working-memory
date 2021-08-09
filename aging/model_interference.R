
##################################
# This model obtains parameters for the sub-cognitive/pre-activation processes

model1 <- "
functions {

  real subcog_lpdf(real[] rt,int md,int mm,int[] accuracy,int[] ltype,
       real alpha,real kappa,row_vector prevector,matrix zeromatrix,real rt_bound,
       int n){
       
       matrix[3,3] cmat;
       row_vector[3] pre;
       real ai;
       real aj;
       real scale_part1;
       real scale_part2;
       real density;

       density = 0;
       
       for(i in 1:md){
            if(rt[i] < 0.6){
                 density = density + weibull_lpdf(rt[i] | kappa,alpha);
            }              
       }
       return(density);
  }
}

data {
  int<lower=1> N;
  int<lower=1> Nsub;
  int<lower=1> Ngroup;
  int Nupt;
  int s_group[Nsub];
  real up_t_list[Nupt];
  int<lower=1> md[N];    // memory demand
  int<lower=1> id[N];    // subject belonging
  int<lower=1> mm[N];
  int accuracy[N,6];       // accuracy of response
  int ltype[N,6];
  real t[N,6];
  int up_t[N];
  real rt_bound[Nsub];
  matrix[3,3] zeromatrix;
  row_vector[3] prevector;
}

parameters {

  real alpha_0;
  real alpha_g[Ngroup];
  real alpha_raw[Nsub];
  
  real kappa_0;
  real kappa_g[Ngroup];
  real kappa_raw[Nsub];
}

transformed parameters {

  real<lower=0> alpha[Nsub];
  real<lower=0> kappa[Nsub];

  for(n in 1:Nsub){
      alpha[n] = exp(alpha_raw[n]);
      kappa[n] = exp(kappa_raw[n]);
  }
}

model {
   int K;
   int M;
   
   alpha_0 ~ normal(0.,1.);
   kappa_0 ~ normal(0.,1.);
   
   for(n in 1:Ngroup){
      alpha_g[n] ~ normal(alpha_0,1.);
      kappa_g[n] ~ normal(kappa_0,1.);
   }
   
   for(n in 1:Nsub){
      K = s_group[n];
      alpha_raw[n] ~ normal(alpha_g[K],1.);
      kappa_raw[n] ~ normal(kappa_g[K],1.);
   } 

   
   for(n in 1:N){
      K = id[n];
      target += subcog_lpdf(t[n] | md[n],mm[n],accuracy[n],ltype[n],alpha[K],
                   kappa[K],prevector,zeromatrix,rt_bound[K],
                   n);
   }
}
"

require(rstan)

set.seed(117)

load("dat_interference_full.Rdata")

dat$keep <- sample(0:1,dat$N,replace = T)


Nsub <- dat$Nsub
rt_bound <- rep(0,Nsub)
dat1 <- data.frame(t=dat$t,id=dat$id)
for(n in 1:Nsub){
  sub <- subset(dat1,id==n)
  t <- c(sub$t.1,sub$t.2,sub$t.3,
         sub$t.4,sub$t.5,sub$t.6)
  t <- t[t>0]
  rt_bound[n] <- quantile(t,probs=0.1)
}



dat$zeromatrix <- matrix(rep(0,9),nrow=3)
dat$prevector <- rep(0,3)
dat$rt_bound <- rt_bound


fit <- stan(model_code=model1,data=dat,chains=1,iter=2000,warmup=1000,
            control = list(adapt_delta = 0.99,max_treedepth=10))

# Obtain parameter estimates to plug in the main model

post <- extract(fit)
kappa <- apply(post$kappa,2,mean)
alpha <- apply(post$alpha,2,mean)

dat$kappa_value <- kappa
dat$alpha_value <- alpha





##################################
# Main model


model2 <- "
functions {

  real mc_lpdf(real[] rt,int md,int mm,int[] accuracy,int[] ltype,
       real alpha,real kappa,vector a,real q,row_vector prevector,matrix zeromatrix,real psi,real rt_bound,
       row_vector compute_1,row_vector compute_2,real kappa_value,real alpha_value){
       
       matrix[3,3] cmat;
       row_vector[3] pre;
       real ai;
       real aj;
       real scale_part1;
       real scale_part2;
       real density;

       density = 0;
       pre = prevector;
       pre[1] = 1;
       
/////////////////

       cmat = zeromatrix;
       cmat[1] = compute_1;  
       cmat[2] = [1./9.,(md-1.)/9.,(9.-md)/9.];
       cmat[3] = [1./9.,(md-1.)/9.,(9.-md)/9.];       
       
       for(i in 1:mm){
             pre = pre*cmat;          
       }
       
       
       cmat[1] = compute_2;  
       
       pre = pre*cmat;
       
       for(i in 1:md){
          scale_part1 = (alpha*pre[ltype[i]])^(1/kappa);
          scale_part2 = (alpha-alpha*pre[ltype[i]])^(1/kappa);
          if(rt[i]>rt_bound){
             density = density + log_sum_exp(log_sum_exp(log(a[1]) +weibull_lpdf(rt[i] | kappa,1/scale_part1)
                    + log(1.00001-weibull_cdf(rt[i],kappa,1/scale_part2)),
                    log(a[2]) + exponential_lpdf(rt[i]-rt_bound | psi) + binomial_lpmf(accuracy[i] | 1,q)),
                    log(a[3]) + weibull_lpdf(rt[i] | kappa_value,alpha_value) + binomial_lpmf(accuracy[i] | 1,q));          
          } 
          else{
             density = density + log_sum_exp(log(a[1]) +weibull_lpdf(rt[i] | kappa,1/scale_part1)
                    + log(1.00001-weibull_cdf(rt[i],kappa,1/scale_part2)),
                    log(a[3]) + weibull_lpdf(rt[i] | kappa_value,alpha_value) + binomial_lpmf(accuracy[i] | 1,q));               
          }
       }
       return(density);
  }
}

data {
  int<lower=1> N;         // trial number
  int<lower=1> Nsub;      // participant number
  int<lower=1> Ngroup;    // group number
  int Nupt;               // number of bins of time limit
  int s_group[Nsub];      // group identifier of 
  real up_t_list[Nupt];   // updating time limit of each bin
  int<lower=1> md[N];     // memory demand
  int<lower=1> id[N];     // subject belonging
  int<lower=1> mm[N];     // number of updates
  int accuracy[N,6];      // accuracy of response
  int ltype[N,6];         // type of response (target, competitor, non-competitor)
  real t[N,6];            // RT
  int up_t[N];            // identifies the bin of the time limit
  real rt_bound[Nsub];    // threshold for supra-cognitive processes (b)
  matrix[3,3] zeromatrix; // matrix of 0s
  row_vector[3] prevector;
  real kappa_value[Nsub]; // Estimated parameters from Model 1
  real alpha_value[Nsub]; // Estimated parameters from Model 1
  int keep[N];            // used in the out-of-sample validation
}

parameters {

  real C_0;               
  real C_g[Ngroup];
  real C_raw[Nsub];
  
  real r_0;               
  real r_g[Ngroup];
  real r_raw[Nsub];
  
  real sigma_0;               
  real sigma_g[Ngroup];
  real sigma_raw[Nsub];
  
  real alpha_0;
  real alpha_g[Ngroup];
  real alpha_raw[Nsub];
  
  real kappa_0;
  real kappa_g[Ngroup];
  real kappa_raw[Nsub];
  
  real q_0;
  real q_g[Ngroup];
  real q_raw[Nsub];
  
  simplex[3] a[Nsub];
  
  real psi_0;
  real psi_g[Ngroup];
  real psi_raw[Nsub];
}

transformed parameters {
  real<lower=0,upper=1> C[Nsub];
  real<lower=0> sigma[Nsub];
  real<lower=0> r[Nsub];
  real<lower=0> alpha[Nsub];
  real<lower=0> kappa[Nsub];
  real<lower=0,upper=1> q[Nsub];
  real<lower=0> psi[Nsub];

  for(n in 1:Nsub){
      C[n] = inv_logit(C_raw[n]);
      sigma[n] = exp(sigma_raw[n]);
      r[n] = exp(r_raw[n]);
      alpha[n] = exp(alpha_raw[n]);
      kappa[n] = exp(kappa_raw[n]);
      q[n] = inv_logit(q_raw[n]);  
      psi[n] = exp(psi_raw[n]);
  }
}

model {
   int K;
   int M;
   
   row_vector[3] compute_once[Nsub,4,Nupt];
   row_vector[3] compute_once_no_time[Nsub,4];
   real ai;
   real aj;

   C_0 ~ normal(0.,1.);
   sigma_0 ~ normal(0.,1.);
   alpha_0 ~ normal(0.,1.);
   kappa_0 ~ normal(0.,1.);
   q_0 ~ normal(0.,1.);
   r_0 ~ normal(0.,1.);
   psi_0 ~ normal(0.,1.);
   
   for(n in 1:Ngroup){
      C_g[n] ~ normal(C_0,1.);
      sigma_g[n] ~ normal(sigma_0,1.);
      alpha_g[n] ~ normal(alpha_0,1.);
      kappa_g[n] ~ normal(kappa_0,1.);
      q_g[n] ~ normal(q_0,1.);
      r_g[n] ~ normal(r_0,1.);
      psi_g[n] ~ normal(psi_0,1.);
   }
   
   for(n in 1:Nsub){
      K = s_group[n];
      C_raw[n] ~ normal(C_g[K],1.);
      sigma_raw[n] ~ normal(sigma_g[K],1.);
      alpha_raw[n] ~ normal(alpha_g[K],1.);
      kappa_raw[n] ~ normal(kappa_g[K],1.);
      q_raw[n] ~ normal(q_g[K],1.);
      r_raw[n] ~ normal(r_g[K],1.);
      psi_raw[n] ~ normal(psi_g[K],1.);
   } 

   for(n in 1:Nsub){
      for(i in 1:4){
         for(j in 1:Nupt){
            if(i>1){
              ai = (1-C[n]/2)^(i-1)*(1-exp(-r[n]*up_t_list[j]));
              aj = (C[n]/2)*(1-C[n]/2)^(i-2)*(1-exp(-r[n]*up_t_list[j]));              
            }
            else{
              ai = (1-exp(-r[n]*up_t_list[j]));
              aj = 0;
            }
            compute_once[n,i,j,1] = exp(ai/sigma[n])/(exp(ai/sigma[n])+(i-1)*exp(aj/sigma[n])+(9-i));
            compute_once[n,i,j,2] = (i-1)*exp(aj/sigma[n])/(exp(ai/sigma[n])+(i-1)*exp(aj/sigma[n])+(9-i));
            compute_once[n,i,j,3] = (9-i)/(exp(ai/sigma[n])+(i-1)*exp(aj/sigma[n])+(9-i));
         }
         if(i>1){
            ai = (1-C[n]/2)^(i-1);
            aj = (C[n]/2)*(1-C[n]/2)^(i-2);              
         }
         else{
            ai = 1;
            aj = 0;
         }
         compute_once_no_time[n,i,1] = exp(ai/sigma[n])/(exp(ai/sigma[n])+(i-1)*exp(aj/sigma[n])+(9-i));
         compute_once_no_time[n,i,2] = (i-1)*exp(aj/sigma[n])/(exp(ai/sigma[n])+(i-1)*exp(aj/sigma[n])+(9-i));
         compute_once_no_time[n,i,3] = (9-i)/(exp(ai/sigma[n])+(i-1)*exp(aj/sigma[n])+(9-i));         
      }
   }
   
   
   for(n in 1:N){
      K = id[n];
      target += mc_lpdf(t[n] | md[n],mm[n],accuracy[n],ltype[n],alpha[K],
                   kappa[K],a[K],q[K],prevector,zeromatrix,psi[K],rt_bound[K],
                   compute_once[K,md[n],up_t[n]],compute_once_no_time[K,md[n]],
                   kappa_value[K],alpha_value[K]);
   }
}
"

require(rstan)

load("dat_interference_full.Rdata")
set.seed(117)

Nsub <- dat$Nsub
rt_bound <- rep(0,Nsub)
dat1 <- data.frame(t=dat$t,id=dat$id)
for(n in 1:Nsub){
  sub <- subset(dat1,id==n)
  t <- c(sub$t.1,sub$t.2,sub$t.3,
         sub$t.4,sub$t.5,sub$t.6)
  t <- t[t>0]
  rt_bound[n] <- max(quantile(t,probs=0.95),2)
}

dat$kappa_value <- kappa
dat$alpha_value <- alpha


dat$zeromatrix <- matrix(rep(0,9),nrow=3)
dat$prevector <- rep(0,3)
dat$rt_bound <- rt_bound


fit <- stan(model_code=model2,data=dat,chains=1,iter=6000,warmup=1000,
            control = list(adapt_delta = 0.8,max_treedepth=10))




















