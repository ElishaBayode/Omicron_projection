negbin.loglik <- function (params,times, test_prop) {
  x <-   as.data.frame(deSolve::ode(y=init,time=times,func= sveirs,
                                    parms= params)) 
  prediction <-  (x$Er+ x$Erv + x$Erw
                  + x$Em+ x$Emv + x$Emw)*test_prop 
  sum(dnbinom(x=dat_omic$value,
              mu=params["p"]*prediction,size=1/params["theta"],
              log=TRUE))
}

func_loglik <- function (par,test_prop, dat_omic) {
  #params <- c(parameters, beta_r= exp(par[1]),beta_m=exp(par[2]), p=expit(par[3]),theta=exp(par[4]))
  params <- c(parameters, beta_m=(par[1]), p=(par[2]),theta=0.1)
  times <- 1:nrow(dat_omic)
  -negbin.loglik(params, times =  times, test_prop)
}


negbin.loglik_2 <- function (params) {
  x <-   as.data.frame(deSolve::ode(y=init,time=times,func= sveirs,
                                    parms= params)) 
  prediction <-  (x$Er+ x$Erv + x$Erw
                  + x$Em+ x$Emv + x$Emw)
  sum(dnbinom(x=dat_omic$value,
              mu=params["p"]*prediction,size=1/params["theta"],
              log=TRUE))
}

func_loglik_2 <- function (par,dat_omic) {
  params <- c(parameters, beta_r= exp(par[1]),beta_m=exp(par[2]), p=expit(par[3]),theta=exp(par[4]))
  -negbin.loglik_2(params)
}





