


##################################
# This model obtains parameters for the sub-cognitive/pre-activation processes



model1 <- "
functions {

  real subcog_lpdf(vector rt,int md,real alpha,real kappa){
       
       real density;
       density = 0;
       
       for(i in 1:9){
             if(rt[i]<0.65){
                density = density + weibull_lpdf(rt[i] | kappa,alpha);      
             }
       }
       
       for(i in 1:md){
             if(rt[9+i]<0.65){
                density = density + weibull_lpdf(rt[9+i] | kappa,alpha);      
             }
       }
       
       return(density);
  }
}

data {
  int<lower=1> N;
  int<lower=1> Nsub;
  int<lower=1> Ngroup;
  int<lower=1> md[N];    // memory demand
  int<lower=0> id[N];    // subject belonging
  int accuracy[N,14];       // accuracy of response
  int pre_accuracy[N,14];
  int pre_num[N,14];
  int answer_num[N,14];
  int ope1[N,14];
  int ope2[N,14];
  int recall[N,14];
  int loc[N,14];
  int swit_rep[N,14];
  int pre_loc[N,5];
  int ltype[N,14];
  int task[N];
  vector[14] rt[N];
  real rt_bound[Nsub];
  matrix[3,3] randommatrix;
  row_vector[3] prevector;
}

parameters {

  real alpha_0;
  real kappa_0;
}

transformed parameters {
  real<lower=0> alpha;
  real<lower=0> kappa;


  for(n in 1:Nsub){

         alpha = exp(alpha_0);
         kappa = exp(kappa_0);

  }
}

model {
   int K;
   int M;

   alpha_0 ~ normal(0.,1.);
   kappa_0 ~ normal(0.,1.);

   
   for(n in 1:N){
      K = id[n];
      M = 2-task[n];
      target += subcog_lpdf(rt[n] | md[n],alpha,kappa);
   }
}
"

require(rstan)

load("dat_vector_197_loc.Rdata")

set.seed(117)
dat$whether_keep <- sample(0:1,dat$N,replace = T)


Nsub <- dat$Nsub

data <- read.table("data_sub.csv",header=T,sep=",")
dat1 <- data.frame(id=data$id,accuracy=data$accuracy,task=data$task,rt=data$rt/1000)
rt_bound <- rep(0,Nsub)
for(i in 1:Nsub){
  sub <- subset(dat1,id==i) 
  rt_bound[i] <- quantile(sub$rt,probs=0.9)
}

dat$rt_bound <- rt_bound
dat$randommatrix <- matrix(rep(0,9),nrow=3)
dat$prevector <- rep(0,3)



fit <- stan(model_code=model1,data=dat,chains=1,iter=5000,warmup=1000,
            control = list(adapt_delta = 0.8,max_treedepth=10))


# Obtain parameter estimates to plug in the main model

post <- extract(fit)

kappa_value <- rep(mean(post$kappa),Nsub)
alpha_value <- rep(mean(post$alpha),Nsub)







################################3
# the actual chains
#################################

model2 <- "
functions {

  real mc_lpdf(vector rt,int md,int[] pre_num,int[] answer_num,int[] ope1,int[] ope2,
       int[] recall,row_vector prevector,matrix randommatrix,int[] loc,int[] pre_loc,int[] accuracy,real C,
       real sigma,real alpha,real kappa,vector a,real psi,real rt_bound,int[] swit_rep,int[] ltype,real q,
       row_vector compute_1,real kappa_value,real alpha_value){
       
       matrix[3,3] cmat;
       real ai;
       real aj;
       row_vector[3] pre[md];
       real scale_part1;
       real scale_part2;
       real density;
       density = 0;
       
       
       cmat = randommatrix;
       cmat[1] = compute_1;
       cmat[2] = [1./9.,(md-1.)/9.,(9.-md)/9.];
       cmat[3] = [1./9.,(md-1.)/9.,(9.-md)/9.];    
       
       for(i in 1:md){
          pre[i] = prevector;
          pre[i,1] = 1;
       }
       
       for(i in 1:9){
             pre[loc[i]] = pre[loc[i]]*cmat;  
             scale_part1 = (alpha*pre[loc[i],ltype[i]])^(1/kappa);
             scale_part2 = (alpha-alpha*pre[loc[i],ltype[i]])^(1/kappa);  
             if(rt[i]>rt_bound){
                density = density + log_sum_exp(log(a[1]+a[2]) + weibull_lpdf(rt[i] | kappa,1/scale_part1)
                    + log(1.00001-weibull_cdf(rt[i],kappa,1/scale_part2)),
                    log(a[3]) + exponential_lpdf(rt[i]-rt_bound | psi) + binomial_lpmf(accuracy[i] | 1,q));          
             }
             else{
                density = density + log(a[1]+a[2]) + weibull_lpdf(rt[i] | kappa,1/scale_part1)
                    + log(1.00001-weibull_cdf(rt[i],kappa,1/scale_part2));          
             }
       }
       
       for(i in 1:md){
          pre[loc[9+i]] = pre[loc[9+i]]*cmat;
          scale_part1 = (alpha*pre[loc[9+i],ltype[9+i]])^(1/kappa);
          scale_part2 = (alpha-alpha*pre[loc[9+i],ltype[9+i]])^(1/kappa);  
          if(rt[9+i]>rt_bound){
              density = density + log_sum_exp(log_sum_exp(log(a[1]) +weibull_lpdf(rt[9+i] | kappa,1/scale_part1)
                    + log(1.00001-weibull_cdf(rt[9+i],kappa,1/scale_part2)),
                    log(a[3]) + exponential_lpdf(rt[9+i]-rt_bound | psi) + binomial_lpmf(accuracy[9+i] | 1,q)),
                    log(a[2]) + weibull_lpdf(rt[9+i] | kappa_value,alpha_value) + binomial_lpmf(accuracy[9+i] | 1,q));
          }
          else{
               density = density + log_sum_exp(log(a[1]) + weibull_lpdf(rt[9+i] | kappa,1/scale_part1)
                    + log(1.00001-weibull_cdf(rt[9+i],kappa,1/scale_part2)),
                    log(a[2]) + weibull_lpdf(rt[9+i] | kappa_value,alpha_value) + binomial_lpmf(accuracy[9+i] | 1,q));
          }
       }
       
       return(density);
  }
}

data {
  int<lower=1> N;        // trial number
  int<lower=1> Nsub;     // participant number
  int<lower=1> Ngroup;   // group number
  int<lower=1> md[N];    // memory demand
  int<lower=0> id[N];    // subject belonging
  int<lower=0> gid[Nsub]; // group belonging
  int accuracy[N,14];       // accuracy of response
  int pre_accuracy[N,14];   // accuracy of previous updating step
  int pre_num[N,14];        // response in the previous step
  int answer_num[N,14];     // response in the current step
  int ope1[N,14];           // type of updating operation: 1 add, 2 substract
  int ope2[N,14];           // amount added or substracted
  int recall[N,14];         // whether the operation is a recall without updating
  int loc[N,14];            // location updated & recalled
  int swit_rep[N,14];       
  int pre_loc[N,5];
  int ltype[N,14];          // response of item types
  int task[N];              // whether it is pre-test (1) or post-test (2)
  vector[14] rt[N];         // RT
  real rt_bound[Nsub];      // RT threshold for supra-cognitive processes
  matrix[3,3] randommatrix; 
  row_vector[3] prevector;
  
//
  int whether_keep[N];
  
  real kappa_value[Nsub];
  real alpha_value[Nsub];
}

parameters {

  real C_0;               // shift
  real C_g[Ngroup];
  real C_raw[Nsub];
  
  real sigma_0;               // direction
  real sigma_g[Ngroup];
  real sigma_raw[Nsub];
  
  real alpha_0;
  real alpha_g[Ngroup];
  real alpha_raw[Nsub];
  
  real kappa_0;
  real kappa_g[Ngroup];
  real kappa_raw[Nsub];
  
  simplex[3] a[Nsub];

  real C1_0;               // shift
  real C1_g[Ngroup];
  real C1_raw[Nsub];
  
  real sigma1_0;               // direction
  real sigma1_g[Ngroup];
  real sigma1_raw[Nsub];
  
  real alpha1_0;
  real alpha1_g[Ngroup];
  real alpha1_raw[Nsub];
  
  real kappa1_0;
  real kappa1_g[Ngroup];
  real kappa1_raw[Nsub];
  
  simplex[3] a1[Nsub];
  
  real q_0;
  real q_g[Ngroup];
  real q_raw[Nsub];

  real psi_0;
}

transformed parameters {
  real<lower=0,upper=1> C[Nsub];
  real<lower=0> sigma[Nsub];
  real<lower=0> alpha[Nsub];
  real<lower=0> kappa[Nsub];
  real<lower=0,upper=1> q[Nsub];

  real<lower=0,upper=1> C1[Nsub];
  real<lower=0> sigma1[Nsub];
  real<lower=0> alpha1[Nsub];
  real<lower=0> kappa1[Nsub];
  real<lower=0,upper=1> q1[Nsub];
  
  real<lower=0> psi;

  for(n in 1:Nsub){
         C[n] = inv_logit(C_raw[n]);
         sigma[n] = exp(sigma_raw[n]);
         alpha[n] = exp(alpha_raw[n]);
         kappa[n] = exp(kappa_raw[n]);
         q[n] = inv_logit(q_raw[n]);  
         
         C1[n] = inv_logit(C1_raw[n]);
         sigma1[n] = exp(sigma1_raw[n]);
         alpha1[n] = exp(alpha1_raw[n]);
         kappa1[n] = exp(kappa1_raw[n]);
         q1[n] = inv_logit(q_raw[n]); 
         
         psi = exp(psi_0);
  }
}

model {
   int K;
   int M;
   row_vector[3] cmat[Nsub,2,2];
   real ai;
   real aj;

   C_0 ~ normal(0.,0.3);
   sigma_0 ~ normal(0.,0.3);
   alpha_0 ~ normal(0.,0.3);
   kappa_0 ~ normal(0.,0.3);
   q_0 ~ normal(0.,0.3);
   
   C1_0 ~ normal(0.,0.3);
   sigma1_0 ~ normal(0.,0.3);
   alpha1_0 ~ normal(0.,0.3);
   kappa1_0 ~ normal(0.,0.3);

   psi_0 ~ normal(0.,0.3);
   
   for(n in 1:Ngroup){
      C_g[n] ~ normal(C_0,0.3);
      sigma_g[n] ~ normal(sigma_0,0.3);
      alpha_g[n] ~ normal(alpha_0,0.3);
      kappa_g[n] ~ normal(kappa_0,0.3);  
      q_g[n] ~ normal(q_0,0.3);
      
      C1_g[n] ~ normal(C1_0,0.3);
      sigma1_g[n] ~ normal(sigma1_0,0.3);
      alpha1_g[n] ~ normal(alpha1_0,0.3);
      kappa1_g[n] ~ normal(kappa1_0,0.3);
   }
   for(n in 1:Nsub){
      K = gid[n];
      C_raw[n] ~ normal(C_g[K],0.3);
      sigma_raw[n] ~ normal(sigma_g[K],0.3);
      alpha_raw[n] ~ normal(alpha_g[K],0.3);
      kappa_raw[n] ~ normal(kappa_g[K],0.3);
      q_raw[n] ~ normal(q_g[K],0.3);
      
      C1_raw[n] ~ normal(C1_g[K],0.3);
      sigma1_raw[n] ~ normal(sigma1_g[K],0.3);
      alpha1_raw[n] ~ normal(alpha1_g[K],0.3);
      kappa1_raw[n] ~ normal(kappa1_g[K],0.3);
   }

   for(n in 1:Nsub){
      for(l in 1:2){
          ai = (1-C[n]/2)^(1+2*l-1);
          aj = (C[n]/2)*(1-C[n]/2)^(1+2*l-2);
          cmat[n,l,1,1] = exp(ai/sigma[n])/(exp(ai/sigma[n])+(1+2*l-1)*exp(aj/sigma[n])+(9-1-2*l));
          cmat[n,l,1,2] = (1+2*l-1)*exp(aj/sigma[n])/(exp(ai/sigma[n])+(1+2*l-1)*exp(aj/sigma[n])+(9-1-2*l));          
          cmat[n,l,1,3] = (9-1-2*l)/(exp(ai/sigma[n])+(1+2*l-1)*exp(aj/sigma[n])+(9-1-2*l));   
          ai = (1-C1[n]/2)^(1+2*l-1);
          aj = (C1[n]/2)*(1-C1[n]/2)^(1+2*l-2);
          cmat[n,l,2,1] = exp(ai/sigma1[n])/(exp(ai/sigma1[n])+(1+2*l-1)*exp(aj/sigma1[n])+(9-1-2*l));
          cmat[n,l,2,2] = (1+2*l-1)*exp(aj/sigma1[n])/(exp(ai/sigma1[n])+(1+2*l-1)*exp(aj/sigma1[n])+(9-1-2*l));          
          cmat[n,l,2,3] = (9-1-2*l)/(exp(ai/sigma1[n])+(1+2*l-1)*exp(aj/sigma1[n])+(9-1-2*l)); 
      }
   }
   
   for(n in 1:N){
      K = id[n];
      M = 2-task[n];
      if(whether_keep[n]==1){
      target += mc_lpdf(rt[n] | md[n],pre_num[n],answer_num[n],ope1[n],ope2[n],
                   recall[n],prevector,randommatrix,loc[n],pre_loc[n],accuracy[n],
                   M*C[K]+(1-M)*C1[K],M*sigma[K]+(1-M)*sigma1[K],M*alpha[K]+(1-M)*alpha1[K],
                   M*kappa[K]+(1-M)*kappa1[K],M*a[K]+(1-M)*a1[K],psi,
                   rt_bound[K],swit_rep[n],ltype[n],M*q[K]+(1-M)*q1[K],cmat[K,(md[n]-1)/2,task[n]],
                   kappa_value[K],alpha_value[K]);
      }
   }
}
"

require(rstan)

load("dat_vector_197_loc.Rdata")

set.seed(117)

Nsub <- dat$Nsub

data <- read.table("data_sub.csv",header=T,sep=",")
dat1 <- data.frame(id=data$id,accuracy=data$accuracy,task=data$task,rt=data$rt/1000)
rt_bound <- rep(0,Nsub)
for(i in 1:Nsub){
  sub <- subset(dat1,id==i) 
  rt_bound[i] <- quantile(sub$rt,probs=0.9)
}

dat$rt_bound <- rt_bound
dat$randommatrix <- matrix(rep(0,9),nrow=3)
dat$prevector <- rep(0,3)

dat$kappa_value <- kappa_value
dat$alpha_value <- alpha_value


initial_values <- list(list(C_raw=rep(0,Nsub),C1_raw=rep(0,Nsub),
                            sigma_raw=rep(0,Nsub),sigma1_raw=rep(0,Nsub),
                            alpha_raw=rep(0,Nsub),alpha1_raw=rep(0,Nsub),
                            kappa_raw=rep(0,Nsub),kappa1_raw=rep(0,Nsub),
                            q_raw=rep(0,Nsub)))


fit <- stan(model_code=model2,data=dat,chains=1,iter=6000,warmup=1000,
            control = list(adapt_delta = 0.8,max_treedepth=10),init=initial_values)


