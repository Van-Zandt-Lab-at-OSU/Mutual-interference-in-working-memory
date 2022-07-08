


#####################
#####################
# one subject

# different distributions 
# 23 0.2  46 0.2 74 0.2


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
  real<lower=1.5> mu2[Nsub];
  
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

# bar plot

# posterior predictive

# compute and plot the accuracy of each subject


require(rstan)
require(statmod)

compute_matrix <- function(C,sigma,md,M,r,up_t){
  
  cmat <- matrix(rep(0,9),nrow=3)
  ai <- (1-C/2)^(md-1)*(1-exp(-r*up_t))
  aj <- ifelse(md>1,(C/2)*(1-C/2)^(md-2)*(1-exp(-r*up_t)),0)
  
  cmat[1,1] = exp(ai/sigma)/(exp(ai/sigma)+(md-1)*exp(aj/sigma)+(M-md))
  cmat[1,2] = (md-1)*exp(aj/sigma)/(exp(ai/sigma)+(md-1)*exp(aj/sigma)+(M-md))
  cmat[1,3] = (M-md)/(exp(ai/sigma)+(md-1)*exp(aj/sigma)+(M-md))  
  cmat[2,] = c(1/M,(md-1)/M,(M-md)/M)
  cmat[3,] = c(1/M,(md-1)/M,(M-md)/M) 
  
  return(cmat)
}

posterior_predictive <- function(md,mm,sigma,C,
                                 mu1,mu2,tau,v,b,kappa,sdd_mu1,a,q,M,r,up_t){
  
  pre <- rep(0,3)
  ltype <- rep(0,md)
  rt <- rep(0,md)
  cmat <- matrix(rep(0,9),nrow=3)
  # compute posterior predictive accuracy
  
  cmat <- compute_matrix(C,sigma,md,M,r,up_t)
  
  pre[1] = 1
  
  for(i in 1:mm){
    pre = pre %*% cmat
  }
  
  
  for(i in 1:md){
    choice <- sum(rmultinom(1,size=1,prob=a)*(1:3))
    if(choice==1){
      rt[i] <- rlnorm(1,mu1,sdd_mu1)
      if(md==1){
        ltype[i] <- sum(rmultinom(1,size=1,prob=c(q[1],q[3]))*c(1,3))
      }
      else{
        ltype[i] <- sum(rmultinom(1,size=1,prob=q)*(1:3))
      }
    }
    else if(choice==2){
      ltype[i] <- sum(rmultinom(1,size=1,prob=pre)*(1:3))
      rt[i] <- tau + rinvgauss(1,b/exp(v + kappa*pre[ltype[i]]),b^2)
    }
    else{
      rt[i] <- tau + rlnorm(1,mu2,1)
      if(md==1){
        ltype[i] <- sum(rmultinom(1,size=1,prob=c(q[1],q[3]))*c(1,3))
      }
      else{
        ltype[i] <- sum(rmultinom(1,size=1,prob=q)*(1:3))
      }    
    }
    
  }
  
  result <- list(ltype=ltype,rt=rt)
  return(result)
  
}

generate_parameters <- function(N1,N2,m,sd){
  c1 <- rnorm(N1,m[1],exp(sd[1]))
  c2 <- rnorm(N2,m[2],exp(sd[2]))
  result <- c(c1,c2)
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

seed <- 1:50
name <- 1:50
names <- paste("result_",name,".Rdata",sep="")
names1 <- paste("dat_",name,".Rdata",sep="")

for(m in 29:30){


require(rstan)

set.seed(seed[m])

load("dat_interference.Rdata")

dat0 <- dat

sid <- 1:20

Nsub <- length(sid)

# size 1=208
J <- 80



seq0 <- sort(sample(469:1170,J,rep=F))

seq <- rep(0,length(seq0)*Nsub)
for(n in 1:Nsub){
  seq[(J*n-J+1):(J*n)] = seq0 + J*(n-1)
}



new_dat <- list(accuracy=dat$accuracy[seq,],
                md=dat$md[seq],
                mm=dat$mm[seq],
                ltype=dat$ltype[seq,],
                rt=dat$rt[seq,],
                id=dat$id[seq],
                up_t_list=dat$up_t_list,
                up_t=dat$up_t[seq],
                PT=dat$PT[seq])


dat <- new_dat

dat$id <- sort(rep(sid,J))

dat$gid <- c(rep(1,10),rep(2,10))

ltype <- array(0,dim=dim(dat$ltype))
rt <- array(0,dim=dim(dat$rt))



load("result_upd_1.Rdata")

post <- extract(fit)

# obtain the parameters by generating them from prior distributions
N1 <- 10
N2 <- 10

C_c <- apply(post$C_c,2,mean)
sd_C_raw <- apply(post$sd_C_raw,2,mean)
C <- generate_parameters(N1,N2,C_c,sd_C_raw)
C <- 1/(1+exp(-C))

sigma_c <- apply(post$sigma_c,2,mean)
sd_sigma_raw <- apply(post$sd_sigma_raw,2,mean)
sigma <- generate_parameters(N1,N2,sigma_c,sd_sigma_raw)
sigma <- 1/(1+exp(-sigma))

kappa_c <- apply(post$kappa_c,2,mean)
sd_kappa_raw <- apply(post$sd_kappa_raw,2,mean)
kappa <- generate_parameters(N1,N2,kappa_c,sd_kappa_raw)

v_c <- apply(post$v_c,2,mean)
sd_v_raw <- apply(post$sd_v_raw,2,mean)
v <- generate_parameters(N1,N2,v_c,sd_v_raw)

b_c <- apply(post$b_c,2,mean)
sd_b_raw <- apply(post$sd_b_raw,2,mean)
b <- generate_parameters(N1,N2,b_c,sd_b_raw)
b <- exp(b)

r_c <- apply(post$r_c,2,mean)
sd_r_raw <- apply(post$sd_r_raw,2,mean)
r <- generate_parameters(N1,N2,r_c,sd_r_raw)
r <- exp(r)

pick_id <- sample(1:33,Nsub)
K <- sample(1:300,1)
a <- post$a[K,pick_id,]
q <- post$q[K,pick_id,]

mu1 <- rep(-1,Nsub)
mu2 <- rep(3,Nsub)
tau <- rep(0.15,Nsub)
sdd_mu1 <- post$sdd_mu1[K,pick_id]

true_param <- list(C=C,C_c=C_c,sd_C_raw=sd_C_raw,
                   sigma=sigma,sigma_c=sigma_c,sd_sigma_raw=sd_sigma_raw,
                   kappa=kappa,kappa_c=kappa_c,sd_kappa_raw=sd_kappa_raw,
                   v=v,v_c=v_c,sd_v_raw=sd_v_raw,
                   r=r,r_c=r_c,sd_r_raw=sd_r_raw,
                   a=a,q=q,
                   b=b,b_c=b_c,sd_b_raw=sd_b_raw,
                   sdd_mu1=sdd_mu1)


for(j in 1:length(dat$md)){
  n <- dat$id[j]

    
  result <- posterior_predictive(dat$md[j],dat$mm[j],sigma[n],C[n],
                                 mu1[n],mu2[n],tau[n],v[n],b[n],kappa[n],sdd_mu1[n],
                                 a[n,],q[n,],9,r[n],dat$PT[j])
  
  ltype[j,1:(dat$md[j])] <- result$ltype
  rt[j,1:(dat$md[j])] <- result$rt
  
}

dat$ltype <- ltype
dat$rt <- rt
dat$accuracy <- ifelse(ltype==1,1,0)

dat$N <- length(dat$id)
dat$Nsub <- Nsub
dat$Ngroup <- 2
dat$gid <- c(rep(1,10),rep(2,10))
dat$sid <- sid




data.sets <- list(dat=dat,true_param=true_param)

save(data.sets,file=names1[m])

###########################


load(names1[m])

dat <- data.sets$dat


Nsub <- dat$Nsub

dat$N_upt <- 4
dat$tau_mean <- rep(0.15,Nsub)
dat$randommatrix <- matrix(rep(0,9),nrow=3)
dat$prevector <- rep(0,3)


init <- list(list(C_raw=1/(1+exp(-C)),
                  r_raw=log(r),
                  kappa_raw=kappa,
                  v_raw=v,
                  b_raw=log(b),
                  sigma_raw=log(sigma),
                  a=a,
                  q=q,
                  tau=tau,
                  mu1=mu1,
                  mu2=mu2,
                  sdd_mu1_raw=1/(1+exp(-sdd_mu1))))

fit <- stan(model_code=model1,data=dat,chains=1,iter=1500,warmup=500,
            control = list(adapt_delta = 0.8,max_treedepth=10),init=init)


post <- extract(fit)


C <- obtain_quantiles(post$C,2)
sigma <- obtain_quantiles(post$sigma,2)
kappa <- obtain_quantiles(post$kappa,2)
v <- obtain_quantiles(post$v,2)
b <- obtain_quantiles(post$b,2)
a <- obtain_quantiles(post$a,3)
r <- obtain_quantiles(post$r,2)
q <- obtain_quantiles(post$q,3)
sdd_mu1 <- obtain_quantiles(post$sdd_mu1,2)
mu1 <- obtain_quantiles(post$mu1,2)
mu2 <- obtain_quantiles(post$mu2,2)
tau <- obtain_quantiles(post$tau,2)


result <- list(C=C,sigma=sigma,kappa=kappa,v=v,r=r,
               b=b,a=a,q=q,sdd_mu1=sdd_mu1,
               mu1=mu1,mu2=mu2,tau=tau)

data.results <- list(result=result,param=data.sets$true_param,seed=seed[m])

save(data.results,file=names[m])

print(m)

}
