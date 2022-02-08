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



# Penalized maximum likelihood: 
penalized_negbin.loglik <- function (params,times, test_prop, pen.weight, penalties) {
  x <-   as.data.frame(deSolve::ode(y=init,time=times,func= sveirs,
                                    parms= params)) 
  prediction <-  params["p"]*params["sigma"]*(x$Er+ x$Erv + x$Erw
                                              + x$Em+ x$Emv + x$Emw)*test_prop[1:nrow(x)]
  ll <- sum(dnbinom(x=dat_omic$value,
              mu=prediction,size=1/params["theta"],
              log=TRUE))
  # Calculate squared difference between known_prop = 0.5 and proportion resident strain on date_known_prop = Dec 12
  pen <- (((x$Er+x$Erv+x$Erw)/(x$Er+x$Erv+x$Erw + x$Em+x$Emv+x$Emw))[penalties$date_known_prop] - penalties$known_prop)^2
  # Calculate squared difference between known_growth = 0.2 and mutant relative growth rate during period_known_growth
  pen2 <- (get_growth_rate(x, startoffset = penalties$period_known_growth[1], 
           duration = penalties$period_known_growth[2]-penalties$period_known_growth[1])$mutrate - penalties$known_growth)^2
  return(ll - (pen+pen2)*pen.weight/2)
}
func_penloglik <- function (par,test_prop,dat_omic,parameters, pen.size, penalties) {
  
  # 'par' contains the parameters to be fit. Replace their default value in 'parameters'
  if(!all(names(par) %in% names(parameters))) stop('Names of pars to fit do not match names in parameters object')
  parameters[names(par)] <- par
  
  # Convert dates to numbered days
  times <- 1:nrow(dat_omic)
  penalties$date_known_prop <- dat_omic$day[dat_omic$date==penalties$date_known_prop]
  penalties$period_known_growth <- c(dat_omic$day[dat_omic$date==period_known_growth[1]],dat_omic$day[dat_omic$date==period_known_growth[2]])
  
  # Penalty weight calculated as relative to negloglh at initial values
  pen.weight <- -pen.size*negbin.loglik(parameters,times,test_prop) 
  
  -penalized_negbin.loglik(params = parameters, times =  times, test_prop, pen.weight=pen.weight, penalties = penalties)
}




