negbin.loglik_ckeck <- function (params) {
  x <-   out_BC
  prediction <-  (x$Er+ x$Erv + x$Erw
                  + x$Em+ x$Emv + x$Emw)*test_prop
  sum(dnbinom(x=dat_omic$value,
              mu=params["p"]*prediction,size=1/params["theta"],
              log=TRUE))
}

f_loglik_check <- function (par) {
  params <- c(S_0=init[[1]],Er_0=init[[2]],Em_0=init[[3]],Ir_0=init[[4]],
              Im_0=init[[5]],R_0=init[[6]],V_0=init[[7]],Erv_0=init[[8]], 
              Emv_0=init[[9]],Irv_0=init[[10]],Imv_0=init[[11]],Rv_0=init[[12]],
              W_0=init[[13]],Erw_0=init[[14]],Emw_0=init[[15]],Irw_0=init[[16]],
              Imw_0=init[[17]],Rw_0=init[[18]] ,parameters,
              beta_r=exp(par[1]),p=expit(par[2]),beta_m=exp(par[3]),theta=exp(par[4]))
  -negbin.loglik_ckeck(params)
}
