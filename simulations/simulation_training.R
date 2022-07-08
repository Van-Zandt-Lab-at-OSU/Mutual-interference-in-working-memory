


#####################
#####################
# Simulations for both updating and recall


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
  vector[14] rt[N];
  matrix[3,3] randommatrix;
  row_vector[3] prevector;
  real tau_mean[Nsub];
  int md_id[N];
}


parameters {

  real C_0;
  real C_c[Ngroup];
  real C_raw[Nsub];
  real sd_C_raw[Ngroup];

  real v_u_0;
  real v_u_c[Ngroup];
  real v_u[Nsub];
  real sd_v_u_raw[Ngroup];

  real b_0;
  real b_c[Ngroup];
  real b_raw[Nsub];
  real sd_b_raw[Ngroup];
  
  simplex[3] a_u[Nsub];
  
  real v_r_0;
  real v_r_c[Ngroup];
  real v_r[Nsub];
  real sd_v_r_raw[Ngroup];

  real kappa_0;
  real kappa_c[Ngroup];
  real kappa[Nsub];
  real sd_kappa_raw[Ngroup];
  
  simplex[3] a_r[Nsub];
  
  real<lower=0,upper=0.5> tau[Nsub];
  simplex[3] q[Nsub];
  
  real sigma_0;
  real sigma_c[Ngroup];
  real sigma_raw[Nsub];
  real sd_sigma_raw[Ngroup];

  real<upper=-1> mu1[Nsub];
  real<upper=-0.5> mu2[Nsub];
  real<lower=1.5>  mu3[Nsub];
  
  real sdd_mu2_raw[Nsub];

}

transformed parameters {
  real<lower=0,upper=1> C[Nsub];
  real<lower=0> b[Nsub];
  
  real<lower=0> sd_C[Ngroup];
  real<lower=0> sd_v_u[Ngroup];
  real<lower=0> sd_b[Ngroup];

  real<lower=0> sd_v_r[Ngroup];
  real<lower=0> sd_kappa[Ngroup];
  
  real<lower=0,upper=1> sigma[Nsub];
  real<lower=0> sd_sigma[Ngroup];
  real<lower=0,upper=1> sdd_mu2[Nsub];
  
  for(n in 1:Nsub){
      C[n] = inv_logit(C_raw[n]);
      b[n] = exp(b_raw[n]);         
      sigma[n] = inv_logit(sigma_raw[n]);
      sdd_mu2[n] = inv_logit(sdd_mu2_raw[n]);
  }
  
  for(n in 1:Ngroup){
       sd_C[n] = exp(sd_C_raw[n]);
       sd_v_u[n] = exp(sd_v_u_raw[n]);
       sd_b[n] = exp(sd_b_raw[n]);
         
       sd_v_r[n] = exp(sd_v_r_raw[n]);
       sd_kappa[n] = exp(sd_kappa_raw[n]);
       sd_sigma[n] = exp(sd_sigma_raw[n]);
  }
}

model {
   int K;
   int M;
   row_vector[3] cmat[Nsub,2];
   real ai;
   real aj;

     C_0 ~ normal(0.,0.2);
     v_u_0 ~ normal(0.,0.2);
     b_0 ~ normal(0.,0.2);
   
     v_r_0 ~ normal(0.,0.1);
     kappa_0 ~ normal(0.,0.2);
     sigma_0 ~ normal(0.,0.2);
   
   for(n in 1:Ngroup){
         C_c[n] ~ normal(C_0,0.2);
         v_u_c[n] ~ normal(v_u_0,0.2);
         kappa_c[n] ~ normal(kappa_0,0.2);
         v_r_c[n] ~ normal(v_r_0,0.1);
         b_c[n] ~ normal(b_0,0.2);
         
         sigma_c[n] ~ normal(sigma_0,0.2); 
         
         sd_C_raw[n] ~ normal(-1,0.2);
         sd_v_u_raw[n] ~ normal(-1,0.2);
         sd_b_raw[n] ~ normal(-1,0.2);
         sd_v_r_raw[n] ~ normal(-1,0.2);
         sd_kappa_raw[n] ~ normal(-1,0.2);
         sd_sigma_raw[n] ~ normal(-1.,0.2);  
   }
   
   for(n in 1:Nsub){
      K = gid[n];
         C_raw[n] ~ normal(C_c[K],sd_C[K]);
         v_u[n] ~ normal(v_u_c[K],sd_v_u[K]);
         b_raw[n] ~ normal(b_c[K],sd_b[K]);

         v_r[n] ~ normal(v_r_c[K],sd_v_r[K]);
         kappa[n] ~ normal(kappa_c[K],sd_kappa[K]);
         
         sigma_raw[n] ~ normal(sigma_c[K],sd_sigma[K]); 
      
      sdd_mu2_raw[n] ~ normal(0.,1.);
      
      mu1[n] ~ normal(-2,0.05);
      mu2[n] ~ normal(-1,0.05);
      mu3[n] ~ normal(3,0.05);
      
      tau[n] ~ normal(tau_mean[n],0.001);   
   }
   
   for(n in 1:Nsub){
      cmat[n,1] = compute_matrix(C[n],sigma[n],3,9);
      cmat[n,2] = compute_matrix(C[n],sigma[n],5,9);
   }
   

   for(n in 1:N){
      K = id[n];
      target += mc_lpdf(rt[n] | md[n],prevector,randommatrix,loc[n],accuracy[n],ltype[n],
                   mu1[K],mu2[K],mu3[K],tau[K],v_u[K],b[K],v_r[K],kappa[K],
                   sdd_mu2[K],a_u[K],a_r[K],q[K],cmat[K,md_id[n]]);
   }
}
"



# bar plot

# posterior predictive

# compute and plot the accuracy of each subject


require(rstan)
require(statmod)

compute_matrix <- function(C,sigma,md,M){
  
  cmat <- matrix(rep(0,9),nrow=3)
  ai <- (1-C/2)^(md-1)
  aj <- (C/2)*(1-C/2)^(md-2)
  
  cmat[1,1] = exp(ai/sigma)/(exp(ai/sigma)+(md-1)*exp(aj/sigma)+(M-md))
  cmat[1,2] = (md-1)*exp(aj/sigma)/(exp(ai/sigma)+(md-1)*exp(aj/sigma)+(M-md))
  cmat[1,3] = (M-md)/(exp(ai/sigma)+(md-1)*exp(aj/sigma)+(M-md))  
  cmat[2,] = c(1/M,(md-1)/M,(M-md)/M)
  cmat[3,] = c(1/M,(md-1)/M,(M-md)/M) 
  
  return(cmat)
}

posterior_predictive <- function(md,loc,sigma,C,
                                 mu1,mu2,mu3,tau,v_u,b,v_r,kappa,sdd_mu2,
                                 a_u,a_r,q,M){
  
  pre <- array(0,dim=c(5,3))
  ltype <- rep(0,9+md)
  rt <- rep(0,9+md)
  cmat <- matrix(rep(0,9),nrow=3)
  # compute posterior predictive accuracy
  
  cmat <- compute_matrix(C,sigma,md,M)
  
  for(i in 1:md){
    pre[i,1] = 1
  }
  
  for(i in 1:9){
    pre[loc[i],] = pre[loc[i],] %*% cmat  
    choice <- sum(rmultinom(1,size=1,prob=a_u)*(1:3))
    if(choice==1){
      rt[i] <- rlnorm(1,mu1,1)
      ltype[i] <- sum(rmultinom(1,size=1,prob=q)*(1:3))
    }
    else if(choice==2){
      ltype[i] <- sum(rmultinom(1,size=1,prob=pre[loc[i],])*(1:3))
      rt[i] <- tau + rinvgauss(1,b/exp(v_u + kappa*pre[loc[i],ltype[i]]),b^2)
    }
    else{
      rt[i] <- tau + rlnorm(1,mu3,1)
      ltype[i] <- sum(rmultinom(1,size=1,prob=q)*(1:3))      
    }
  }
  
  cmat <- compute_matrix(C,sigma,md,M)
  
  for(i in 1:md){
    pre[loc[9+i],] = pre[loc[9+i],] %*% cmat
    choice <- sum(rmultinom(1,size=1,prob=a_r)*(1:3))
    if(choice==1){
      rt[9+i] <- rlnorm(1,mu2,sdd_mu2)
      ltype[9+i] <- sum(rmultinom(1,size=1,prob=q)*(1:3))
    }
    else if(choice==2){
      ltype[9+i] <- sum(rmultinom(1,size=1,prob=pre[loc[9+i],])*(1:3))
      rt[9+i] <- tau + rinvgauss(1,b/exp(v_r + kappa*pre[loc[9+i],ltype[9+i]]),b^2)
    }
    else{
      rt[9+i] <- tau + rlnorm(1,mu3,1)
      ltype[9+i] <- sum(rmultinom(1,size=1,prob=q)*(1:3))      
    }
    
  }
  
  result <- list(ltype=ltype,rt=rt)
  return(result)
  
}

obtain_quantiles <- function(param,L){
  result <- array(0,dim=c(3,dim(param)[2:L]))
  if(L==2){
    J <- dim(param)[2]
    for(j in 1:J){
      result[,j] <- quantile(param[,j],probs=c(0.025,0.5,0.975))
    }
  }
  else if(L==3){
    J <- dim(param)[2]
    M <- dim(param)[3]
    for(j in 1:J){
      for(m in 1:M){
        result[,j,m] <- quantile(param[,j,m],probs=c(0.025,0.5,0.975))        
      }
    }
  }
  return(result)
}




###########################################

###########################################


# A set up with only one condition
# 1,4,5,6,7,8,9,10,11,12
# 14,16,17,20,22,23,24,25,29,30

seed <- c(1,4,5,6,7,8,9,10,11,12,14,16,17,20,22,23,24,25,29,30)
name <- c(1,4,5,6,7,8,9,10,11,12,14,16,17,20,22,23,24,25,29,30)
names <- paste("result_long_",name,".Rdata",sep="")
names1 <- paste("data_",name,".Rdata",sep="")

for(m in 1:10){


require(rstan)

set.seed(seed[m])

load("dat_vector_197_loc.Rdata")

dat0 <- dat

sid <- 1:15

Nsub <- length(sid)

# size 1=208
J <- 1

seq <- seq(1,J*32*Nsub*3,2)


new_dat <- list(accuracy=dat$accuracy[seq,],
                loc=dat$loc[seq,],
                md=dat$md[seq],
                ltype=dat$ltype[seq,],
                rt=dat$rt[seq,],
                id=dat$id[seq],
                task=rep(1,length(seq)))


dat <- new_dat

dat$id <- sort(rep(sid,J*16*3))

dat$gid <- c(rep(1,5),rep(2,5),rep(3,5))

ltype <- array(0,dim=dim(dat$ltype))
rt <- array(0,dim=dim(dat$rt))



load("model_result_1.Rdata")

post <- extract(fit)

# obtain the parameters by generating them from prior distributions
N1 <- 5
N2 <- 5
N3 <- 5

generate_parameters <- function(N1,N2,N3,m,sd){
  c1 <- rnorm(N1,m[1],exp(sd[1]))
  c2 <- rnorm(N1,m[2],exp(sd[2]))
  c3 <- rnorm(N1,m[3],exp(sd[3]))
  result <- c(c1,c2,c3)
  return(result)
}


C_c <- apply(post$C_c[,,2],2,mean)
sd_C_raw <- apply(post$sd_C_raw[,,2],2,mean)
C <- generate_parameters(N1,N2,N3,C_c,sd_C_raw)
C <- 1/(1+exp(-C))

sigma_c <- apply(post$sigma_c[,,2],2,mean)
sd_sigma_raw <- apply(post$sd_sigma_raw[,,2],2,mean)
sigma <- generate_parameters(N1,N2,N3,sigma_c,sd_sigma_raw)
sigma <- 1/(1+exp(-sigma))

kappa_c <- apply(post$kappa_u_c[,,2],2,mean)
sd_kappa_raw <- apply(post$sd_kappa_u_raw[,,2],2,mean)
kappa <- generate_parameters(N1,N2,N3,kappa_c,sd_kappa_raw)

v_u_c <- apply(post$v_u_c[,,2],2,mean)
sd_v_u_raw <- apply(post$sd_v_u_raw[,,2],2,mean)
v_u <- generate_parameters(N1,N2,N3,v_u_c,sd_v_u_raw)

v_r_c <- apply(post$v_r_c[,,2],2,mean)
sd_v_r_raw <- apply(post$sd_v_r_raw[,,2],2,mean)
v_r <- generate_parameters(N1,N2,N3,v_r_c,sd_v_r_raw)

b_c <- apply(post$b_c[,,2],2,mean)
sd_b_raw <- apply(post$sd_b_raw[,,2],2,mean)
b <- generate_parameters(N1,N2,N3,b_c,sd_b_raw)
b <- exp(b)

pick_id <- sample(1:196,15)
K <- sample(1:300,1)
a_u <- post$a_u[K,pick_id,2,]
a_r <- post$a_r[K,pick_id,2,]
q <- post$q[K,pick_id,2,]

mu1 <- rep(-2,15)
mu2 <- rep(-1,15)
mu3 <- rep(3,15)
tau <- rep(0.15,15)
sdd_mu2 <- post$sdd_mu2[K,pick_id]

true_param <- list(C=C,C_c=C_c,sd_C_raw=sd_C_raw,
                   sigma=sigma,sigma_c=sigma_c,sd_sigma_raw=sd_sigma_raw,
                   kappa=kappa,kappa_c=kappa_c,sd_kappa_raw=sd_kappa_raw,
                   v_u=v_u,v_u_c=v_u_c,sd_v_u_raw=sd_v_u_raw,
                   v_r=v_r,v_r_c=v_r_c,sd_v_r_raw=sd_v_r_raw,
                   a_u=a_u,a_r=a_r,q=q,
                   b=b,b_c=b_c,sd_b_raw=sd_b_raw,
                   sdd_mu2=sdd_mu2)


for(j in 1:length(dat$md)){
  n <- dat$id[j]
  K <- dat$task[j]

    
  result <- posterior_predictive(dat$md[j],dat$loc[j,],sigma[n],C[n],
                                 mu1[n],mu2[n],mu3[n],tau[n],v_u[n],b[n],v_r[n],kappa[n],
                                 sdd_mu2[n],a_u[n,],a_r[n,],q[n,],9)
  
  ltype[j,1:(9+dat$md[j])] <- result$ltype
  rt[j,1:(9+dat$md[j])] <- result$rt
  
}

dat$ltype <- ltype
dat$rt <- rt
dat$accuracy <- ifelse(ltype==1,1,0)

dat$N <- length(dat$id)
dat$Nsub <- 15
dat$Ngroup <- 3
dat$gid <- c(rep(1,5),rep(2,5),rep(3,5))
dat$sid <- sid




data.sets <- list(dat=dat,true_param=true_param)

save(data.sets,file=names1[m])

###########################


load(names1[m])

dat <- data.sets$dat


Nsub <- dat$Nsub


dat$tau_mean <- rep(0.15,15)
dat$md_id <- ifelse(dat$md==3,1,2)
dat$randommatrix <- matrix(rep(0,9),nrow=3)
dat$prevector <- rep(0,3)



fit <- stan(model_code=model1,data=dat,chains=1,iter=1500,warmup=500,
            control = list(adapt_delta = 0.8,max_treedepth=10))


post <- extract(fit)


C <- obtain_quantiles(post$C,2)
sigma <- obtain_quantiles(post$sigma,2)
kappa <- obtain_quantiles(post$kappa,2)
v_u <- obtain_quantiles(post$v_u,2)
v_r <- obtain_quantiles(post$v_r,2)
b <- obtain_quantiles(post$b,2)
a_u <- obtain_quantiles(post$a_u,3)
a_r <- obtain_quantiles(post$a_r,3)
q <- obtain_quantiles(post$q,3)
sdd_mu2 <- obtain_quantiles(post$sdd_mu2,2)
mu1 <- obtain_quantiles(post$mu1,2)
mu2 <- obtain_quantiles(post$mu2,2)
mu3 <- obtain_quantiles(post$mu3,2)
tau <- obtain_quantiles(post$tau,2)


result <- list(C=C,sigma=sigma,kappa=kappa,v_u=v_u,v_r=v_r,
               b=b,a_u=a_u,a_r=a_r,q=q,sdd_mu2=sdd_mu2,
               mu1=mu1,mu2=mu2,mu3=mu3,tau=tau)

data.results <- list(result=result,param=data.sets$true_param,seed=seed[m])

save(data.results,file=names[m])

print(m)

}
