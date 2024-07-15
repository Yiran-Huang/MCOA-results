library(glmnet)
lassopre<-function(X,Y){
  X=X[,-1]
  alpha=1
  lasso_model <- cv.glmnet(X,Y, alpha = alpha)
  best_lambda <- lasso_model$lambda.min
  lasso_fit <- glmnet(X,Y, alpha = alpha, lambda = best_lambda)
  beta=matrix(as.matrix(coef(lasso_fit)),ncol=1)
  return(beta)
}


simulation_function<-function(xcoa,xmcoa,z,fun){
  ALL=c()
  for(i in 1:10000){
    ALL=rbind(ALL,sample(0:10,11))
  }
  
  
  
  
  XCOA_loc=xcoa+1
  Xmaximin_loc=xmcoa+1
  len_sample_COA=length(XCOA_loc[,1])
  len_sample_maximin=length(Xmaximin_loc[,1])
  DCOA_loc=z
  DCOA_loc=DCOA_loc+1
  len_sample_DCOA=length(DCOA_loc[,1])
  
  XCOA=component_2_location(XCOA_loc)-1
  Xmaximin=component_2_location(Xmaximin_loc)-1
  DCOA=component_2_location(DCOA_loc)-1
  ALL_loc=component_2_location(ALL+1)
  XCOA_CP=location_2_CPdesign(XCOA_loc)
  Xmaximin_CP=location_2_CPdesign(Xmaximin_loc)
  DCOA_CP=location_2_CPdesign(DCOA_loc)
  ALL_CP=location_2_CPdesign(ALL_loc)
  XCOA_PWO=location_2_PWOdesign(XCOA_loc)
  Xmaximin_PWO=location_2_PWOdesign(Xmaximin_loc)
  DCOA_PWO=location_2_PWOdesign(DCOA_loc)
  ALL_PWO=location_2_PWOdesign(ALL_loc)
  XCOA_FO=location_2_firstorder(XCOA_loc)
  Xmaximin_FO=location_2_firstorder(Xmaximin_loc)
  DCOA_FO=location_2_firstorder(DCOA_loc)
  ALL_FO=location_2_firstorder(ALL_loc)
  XCOA_QO=location_2_quadratic(XCOA_loc)
  Xmaximin_QO=location_2_quadratic(Xmaximin_loc)
  DCOA_QO=location_2_quadratic(DCOA_loc)
  ALL_QO=location_2_quadratic(ALL_loc)
  XCOA_SO=location_2_secondorder(XCOA_loc)
  Xmaximin_SO=location_2_secondorder(Xmaximin_loc)
  DCOA_SO=location_2_secondorder(DCOA_loc)
  ALL_SO=location_2_secondorder(ALL_loc)
  XCOA_tPWO=location_2_tPWOdesign(XCOA_loc)
  Xmaximin_tPWO=location_2_tPWOdesign(Xmaximin_loc)
  DCOA_tPWO=location_2_tPWOdesign(DCOA_loc)
  ALL_tPWO=location_2_tPWOdesign(ALL_loc)
  XCOA_triPWO=location_2_tripletPWOdesign(XCOA_loc)
  Xmaximin_triPWO=location_2_tripletPWOdesign(Xmaximin_loc)
  DCOA_triPWO=location_2_tripletPWOdesign(DCOA_loc)
  ALL_triPWO=location_2_tripletPWOdesign(ALL_loc)
  Y_ALL=fun(ALL_loc)
  Y_XCOA=fun(XCOA_loc)
  Y_Xmaximin=fun(Xmaximin_loc)
  Y_DCOA=fun(DCOA_loc)
  
  #-------------CP-------------
  beta_XCOACP_hat=lassopre(XCOA_CP, Y_XCOA)
  YCOACP=ALL_CP%*%beta_XCOACP_hat
  MSE_XCOA_CP=mean((Y_ALL-YCOACP)^2)
  beta_XmaximinCP_hat=lassopre(Xmaximin_CP, Y_Xmaximin)
  YmaximinCP=ALL_CP%*%beta_XmaximinCP_hat
  MSE_Xmaximin_CP=mean((Y_ALL-YmaximinCP)^2)
  beta_DCOACP_hat=lassopre(DCOA_CP, Y_DCOA)
  YDCP=ALL_CP%*%beta_DCOACP_hat
  MSE_DCOA_CP=mean((Y_ALL-YDCP)^2)
  
  #-------------PWO-------------
  beta_XCOAPWO_hat=lassopre(XCOA_PWO, Y_XCOA)
  YCOAPWO=ALL_PWO%*%beta_XCOAPWO_hat
  MSE_XCOA_PWO=mean((Y_ALL-YCOAPWO)^2)
  beta_XmaximinPWO_hat=lassopre(Xmaximin_PWO, Y_Xmaximin)
  YmaximinPWO=ALL_PWO%*%beta_XmaximinPWO_hat
  MSE_Xmaximin_PWO=mean((Y_ALL-YmaximinPWO)^2)
  beta_DCOAPWO_hat=lassopre(DCOA_PWO, Y_DCOA)
  YDPWO=ALL_PWO%*%beta_DCOAPWO_hat
  MSE_DCOA_PWO=mean((Y_ALL-YDPWO)^2)
  
  #-------------FO-------------
  beta_XCOAFO_hat=solve(t(XCOA_FO)%*%XCOA_FO)%*%t(XCOA_FO)%*%Y_XCOA
  YCOAFO=ALL_FO%*%beta_XCOAFO_hat
  MSE_XCOA_FO=mean((Y_ALL-YCOAFO)^2)
  beta_XmaximinFO_hat=solve(t(Xmaximin_FO)%*%Xmaximin_FO)%*%t(Xmaximin_FO)%*%Y_Xmaximin
  YmaximinFO=ALL_FO%*%beta_XmaximinFO_hat
  MSE_Xmaximin_FO=mean((Y_ALL-YmaximinFO)^2)
  beta_DCOAFO_hat=solve(t(DCOA_FO)%*%DCOA_FO)%*%t(DCOA_FO)%*%Y_DCOA
  YDFO=ALL_FO%*%beta_DCOAFO_hat
  MSE_DCOA_FO=mean((Y_ALL-YDFO)^2)
  
  #-------------QO-------------
  beta_XCOAQO_hat=solve(t(XCOA_QO)%*%XCOA_QO)%*%t(XCOA_QO)%*%Y_XCOA
  YCOAQO=ALL_QO%*%beta_XCOAQO_hat
  MSE_XCOA_QO=mean((Y_ALL-YCOAQO)^2)
  beta_XmaximinQO_hat=solve(t(Xmaximin_QO)%*%Xmaximin_QO)%*%t(Xmaximin_QO)%*%Y_Xmaximin
  YmaximinQO=ALL_QO%*%beta_XmaximinQO_hat
  MSE_Xmaximin_QO=mean((Y_ALL-YmaximinQO)^2)
  beta_DCOAQO_hat=solve(t(DCOA_QO)%*%DCOA_QO)%*%t(DCOA_QO)%*%Y_DCOA
  YDQO=ALL_QO%*%beta_DCOAQO_hat
  MSE_DCOA_QO=mean((Y_ALL-YDQO)^2)
  
  #-------------SO-------------
  beta_XCOASO_hat=try(solve(t(XCOA_SO)%*%XCOA_SO)%*%t(XCOA_SO)%*%Y_XCOA,TRUE)
  if(mode(beta_XCOASO_hat)=="numeric"){
    YCOASO=ALL_SO%*%beta_XCOASO_hat
    MSE_XCOA_SO=mean((Y_ALL-YCOASO)^2)}else{MSE_XCOA_SO=Inf;YCOASO=NULL}
  beta_XmaximinSO_hat=solve(t(Xmaximin_SO)%*%Xmaximin_SO)%*%t(Xmaximin_SO)%*%Y_Xmaximin
  YmaximinSO=ALL_SO%*%beta_XmaximinSO_hat
  MSE_Xmaximin_SO=mean((Y_ALL-YmaximinSO)^2)
  beta_DCOASO_hat=solve(t(DCOA_SO)%*%DCOA_SO)%*%t(DCOA_SO)%*%Y_DCOA
  YDSO=ALL_SO%*%beta_DCOASO_hat
  MSE_DCOA_SO=mean((Y_ALL-YDSO)^2)
  
  #-------------tPWO-------------
  beta_XCOAtPWO_hat=try(solve(t(XCOA_tPWO)%*%XCOA_tPWO)%*%t(XCOA_tPWO)%*%Y_XCOA,TRUE)
  if(mode(beta_XCOAtPWO_hat)=="numeric"){
    YCOAtPWO=ALL_tPWO%*%beta_XCOAtPWO_hat
    MSE_XCOA_tPWO=mean((Y_ALL-YCOAtPWO)^2)}else{MSE_XCOA_tPWO=Inf;YCOAtPWO=NULL}
  beta_XmaximintPWO_hat=solve(t(Xmaximin_tPWO)%*%Xmaximin_tPWO)%*%t(Xmaximin_tPWO)%*%Y_Xmaximin
  YmaximintPWO=ALL_tPWO%*%beta_XmaximintPWO_hat
  MSE_Xmaximin_tPWO=mean((Y_ALL-YmaximintPWO)^2)
  beta_DCOAtPWO_hat=solve(t(DCOA_tPWO)%*%DCOA_tPWO)%*%t(DCOA_tPWO)%*%Y_DCOA
  YDtPWO=ALL_tPWO%*%beta_DCOAtPWO_hat
  MSE_DCOA_tPWO=mean((Y_ALL-YDtPWO)^2)
  
  #-------------tri-------------
  beta_XCOAtriPWO_hat=lassopre(XCOA_triPWO, Y_XCOA)
  YCOAtriPWO=ALL_triPWO%*%beta_XCOAtriPWO_hat
  MSE_XCOA_triPWO=mean((Y_ALL-YCOAtriPWO)^2)
  beta_XmaximintriPWO_hat=lassopre(Xmaximin_triPWO, Y_Xmaximin)
  YmaximintriPWO=ALL_triPWO%*%beta_XmaximintriPWO_hat
  MSE_Xmaximin_triPWO=mean((Y_ALL-YmaximintriPWO)^2)
  beta_DCOAtriPWO_hat=lassopre(DCOA_triPWO, Y_DCOA)
  YDtriPWO=ALL_triPWO%*%beta_DCOAtriPWO_hat
  MSE_DCOA_triPWO=mean((Y_ALL-YDtriPWO)^2)
  
  output=c(MSE_Xmaximin_CP,MSE_XCOA_CP,MSE_DCOA_CP,MSE_Xmaximin_PWO,MSE_XCOA_PWO,
           MSE_DCOA_PWO,MSE_Xmaximin_FO,MSE_XCOA_FO,MSE_DCOA_FO,MSE_Xmaximin_QO,MSE_XCOA_QO,
           MSE_DCOA_QO,MSE_Xmaximin_SO,MSE_XCOA_SO,MSE_DCOA_SO,MSE_Xmaximin_tPWO,MSE_XCOA_tPWO,
           MSE_DCOA_tPWO,MSE_Xmaximin_triPWO,MSE_XCOA_triPWO,MSE_DCOA_triPWO)/c(var(Y_ALL))
  output=matrix(output,nrow=1)
  colnames(output)<-c("MCOA CP","COA CP","D CP",
                   "MCOA PWO","COA PWO","D PWO",
                   "MCOA FO","COA FO","D FO",
                   "MCOA QO","COA QO","D QO",
                   "MCOA SO","COA SO","D SO",
                   "MCOA tPWO","COA tPWO","D tPWO",
                   "MCOA triPWO","COA triPWO","D triPWO")
  return(output)
}


