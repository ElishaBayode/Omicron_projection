negbin.loglik <- function (params,times, test_prop) {
  x <-   as.data.frame(deSolve::ode(y=init,time=times,func= sveirs,
                                    parms= params)) 
  prediction <-  params["p"]*params["sigma"]*(x$Er+ x$Erv + x$Erw
                  + x$Em+ x$Emv + x$Emw)*test_prop[1:nrow(x)]
  sum(dnbinom(x=dat_omic$value,
              mu=prediction,size=1/params["theta"],
              log=TRUE))
}


plot.loglik.info <- function(params,times, test_prop, init) { 
  x <-   as.data.frame(deSolve::ode(y=init,time=times,func= sveirs,
                                    parms= params)) 
  prediction <-  params["p"]*params["sigma"]*(x$Er+ x$Erv + x$Erw
                  + x$Em+ x$Emv + x$Emw)*test_prop[1:nrow(x)]
  tmp = data.frame(time =x$time[1:nrow(dat_omic)], 
                   model = prediction[1:nrow(dat_omic)],
                   data = dat_omic$value)
  return(ggplot(data = tmp, aes(x=time, y=model))+geom_line() + geom_point(aes(x=time, y=data), alpha=0.5))
#  sum(dnbinom(x=dat_omic$value,
#              mu=params["p"]*prediction,size=1/params["theta"],
#              log=TRUE))
}

func_loglik <- function (par,test_prop,dat_omic,parameters,init) {
  # 'par' contains the parameters to be fit. Replace their default value in 'parameters'
  if(!all(names(par) %in% names(parameters))) stop('Names of pars to fit do not match names in parameters object')
  parameters[names(par)] <- par
  times <- 1:nrow(dat_omic)
  -negbin.loglik(params = parameters, times =  times, test_prop)
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

func_loglik_2 <- function (par,dat_omic,parameters) {
  # 'par' contains the parameters to be fit. Replace their default value in 'parameters'
  if(!all(names(par) %in% names(parameters))) stop('Names of pars to fit do not match names in parameters object')
  parameters[names(par)] <- par
  -negbin.loglik_2(parameters)
}





