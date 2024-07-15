#pdist_1_2_many(x,y,p=2): x is a vector and y is a matrix. Find the distance of x to each row of y.
#genetic_on_maximin_mink: genetic algorithm, parameter p is for Lp distance.

pdist_1_2_many<-function(x,y,p=2){
  x=as.numeric(x)
  dif_abs=abs(t(t(y)-x))
  return(rowSums(dif_abs^p)^(1/p))
}

COA_trans_distance<-function(x,loc){
  batch_operator<-function(permutation){
    return(min(dist(component_2_location(x[,permutation]+1)-1)))}
  apply(loc,1,batch_operator)
}

COA_trans_distance_mink<-function(x,loc,p){
  batch_operator<-function(permutation){
    return(min(dist(component_2_location(x[,permutation]+1)-1,method='minkowski',p=p)))}
  apply(loc,1,batch_operator)
}

genetic_on_maximin_mink<-function(x,n_start=NULL,loop_time=100,cross_p=0.5,mutate_p=0.5,adjust_factor=NULL,initial=NULL,p=2){
  flag=check_COA(x)
  if(!flag){print("**************Inputting array is not a COA**************")}
  m=length(x[1,])
  n=length(x[,1])
  lambda=n/m/(m-1)
  if(is.null(n_start)){n_start=m*10}
  if(is.null(adjust_factor)){adjust_factor=lambda*m/2}
  generation<-c()
  x_distance=min(dist(x,method='minkowski',p=p))
  
  
  #generate intial individuals
  if(is.null(initial)){
    for(i in 1:n_start){generation=rbind(generation,c(1,2,sample(3:m,(m-2))))}
  }else{
    generation=initial
    n_start=length(initial[,1])
  }
  
  #fitting 
  original_fitting=COA_trans_distance_mink(x,generation,p)
  fitting=exp((original_fitting-x_distance)/x_distance*adjust_factor)
  
  
  
  #genetic algorithm
  for(loop in 1:loop_time){
    #choosing
    surviving<-sample(1:n_start,size=n_start,prob=fitting,replace = TRUE)
    generation=generation[surviving,]
    
    #crossing
    for(cross in 1:(n_start/2)){
      if(runif(1)<cross_p){
        a=generation[2*cross-1,]
        b=generation[2*cross,]
        loc=sort(sample(3:(m-1),size=2))
        d=b[loc[1]:loc[2]]
        b[loc[1]:loc[2]]=a[loc[1]:loc[2]]
        a[loc[1]:loc[2]]=d
        checklistloc=setdiff(1:m,loc[1]:loc[2])
        for(check in 1:1000){
          breakflag=TRUE
          for(checkloc in checklistloc){
            reploc=which(a[loc[1]:loc[2]]==a[checkloc])
            if(length(reploc)>0){reploc=reploc+loc[1]-1;breakflag=FALSE;break}
          }
          if(breakflag){break}else{a[checkloc]=b[reploc]}}
        for(check in 1:1000){
          breakflag=TRUE
          for(checkloc in checklistloc){
            reploc=which(b[loc[1]:loc[2]]==b[checkloc])
            if(length(reploc)>0){reploc=reploc+loc[1]-1;breakflag=FALSE;break}
          }
          if(breakflag){break}else{b[checkloc]=a[reploc]}}
        generation[2*cross-1,]=a
        generation[2*cross,]=b
      }
    }
    
    #mutate
    for (mtt in 1:n_start) {
      if(runif(1)<mutate_p){
        a=generation[mtt,]
        loc=sample(3:m,2)
        a[c(loc[1],loc[2])]=a[c(loc[2],loc[1])]
        generation[mtt,]=a
      }
    }
    
    #fitting 
    original_fitting=COA_trans_distance_mink(x,generation,p)
    fitting=exp((original_fitting-x_distance)/x_distance*adjust_factor)
    
  }
  output=list(generation,original_fitting,fitting)
  names(output)<-c("Last generation","Distance of each generation","transformed fitting value of generation")
  return(output)
}

permutation_to_PWO_efficiency<-function(initial_design,permutation,m=NULL,n=NULL,p=NULL,D_full=NULL){
  if(is.null(m)){m=ncol(initial_design)}
  if(is.null(n)){n=nrow(initial_design)}
  if(is.null(p)){p=choose(m,2)+1}
  q=p-1
  if(is.null(D_full)){D_full=((m+1)^(m-1)/(3^q))^(1/p)}
  F_PWO=location_2_PWOdesign(component_2_location(initial_design[,permutation]+1))
  D_eff=(det(t(F_PWO)%*%F_PWO/n))^(1/p)
  efficiency=D_eff/D_full
  if(is.na(efficiency)){efficiency=0}
  return(efficiency)
}



genetic_on_PWO_efficiency<-function(x,n_start=NULL,loop_time=100,cross_p=0.5,mutate_p=0.5,adjust_factor=NULL,initial=NULL){
  flag=check_COA(x)
  if(!flag){print("**************Inputting array is not a COA**************")}
  m=length(x[1,])
  n=length(x[,1])
  lambda=n/m/(m-1)
  if(is.null(n_start)){n_start=m*10}
  if(is.null(adjust_factor)){adjust_factor=lambda*m/2}
  generation<-c()
  p=choose(m,2)+1
  q=p-1
  D_full=((m+1)^(m-1)/(3^q))^(1/p)
  batch_PWO_efficiency<-function(permutation){
    permutation_to_PWO_efficiency(initial_design=x,permutation=permutation,m=m,n=n,p=p,D_full=D_full)
  }
  x_efficiency=batch_PWO_efficiency(1:m)
  
  
  #generate intial individuals
  if(is.null(initial)){
    for(i in 1:n_start){generation=rbind(generation,c(1,2,sample(3:m,(m-2))))}
  }else{
    generation=initial
    n_start=length(initial[,1])
  }
  
  #fitting 
  original_fitting=apply(generation,1,batch_PWO_efficiency)
  if(x_efficiency==0){x_efficiency=mean(original_fitting)}
  fitting=exp((original_fitting-x_efficiency)/x_efficiency*adjust_factor)
  
  
  
  #genetic algorithm
  for(loop in 1:loop_time){
    #choosing
    surviving<-sample(1:n_start,size=n_start,prob=fitting,replace = TRUE)
    generation=generation[surviving,]
    
    #crossing
    for(cross in 1:(n_start/2)){
      if(runif(1)<cross_p){
        a=generation[2*cross-1,]
        b=generation[2*cross,]
        loc=sort(sample(3:(m-1),size=2))
        d=b[loc[1]:loc[2]]
        b[loc[1]:loc[2]]=a[loc[1]:loc[2]]
        a[loc[1]:loc[2]]=d
        checklistloc=setdiff(1:m,loc[1]:loc[2])
        for(check in 1:1000){
          breakflag=TRUE
          for(checkloc in checklistloc){
            reploc=which(a[loc[1]:loc[2]]==a[checkloc])
            if(length(reploc)>0){reploc=reploc+loc[1]-1;breakflag=FALSE;break}
          }
          if(breakflag){break}else{a[checkloc]=b[reploc]}}
        for(check in 1:1000){
          breakflag=TRUE
          for(checkloc in checklistloc){
            reploc=which(b[loc[1]:loc[2]]==b[checkloc])
            if(length(reploc)>0){reploc=reploc+loc[1]-1;breakflag=FALSE;break}
          }
          if(breakflag){break}else{b[checkloc]=a[reploc]}}
        generation[2*cross-1,]=a
        generation[2*cross,]=b
      }
    }
    
    #mutate
    for (mtt in 1:n_start) {
      if(runif(1)<mutate_p){
        a=generation[mtt,]
        loc=sample(3:m,2)
        a[c(loc[1],loc[2])]=a[c(loc[2],loc[1])]
        generation[mtt,]=a
      }
    }
    
    #fitting 
    original_fitting=apply(generation,1,batch_PWO_efficiency)
    fitting=exp((original_fitting-x_efficiency)/x_efficiency*adjust_factor)
    
  }
  output=list(generation,original_fitting,fitting)
  names(output)<-c("Last generation","D_efficiency of each generation","transformed fitting value of generation")
  return(output)
}

permutation_to_SO_efficiency<-function(initial_design,permutation,m=NULL,n=NULL){
  if(is.null(m)){m=ncol(initial_design)}
  if(is.null(n)){n=nrow(initial_design)}
  F_SO=location_2_secondorder(component_2_location(initial_design[,permutation]+1))
  eigenall=eigen(t(F_SO)%*%F_SO/nrow(F_SO))[[1]]
  eigenall[eigenall<0]=0
  efficiency=exp(sum(log(eigenall))/nrow(F_SO))
  return(efficiency)
}

permutation_to_tPWO_efficiency<-function(initial_design,permutation,m=NULL,n=NULL){
  if(is.null(m)){m=ncol(initial_design)}
  if(is.null(n)){n=nrow(initial_design)}
  F_tPWO=location_2_tPWOdesign(component_2_location(initial_design[,permutation]+1))
  eigenall=eigen(t(F_tPWO)%*%F_tPWO/nrow(F_tPWO))[[1]]
  eigenall[eigenall<0]=0
  efficiency=exp(sum(log(eigenall))/nrow(F_tPWO))
  return(efficiency)
}

permutation_to_tri_efficiency<-function(initial_design,permutation,m=NULL,n=NULL){
  if(is.null(m)){m=ncol(initial_design)}
  if(is.null(n)){n=nrow(initial_design)}
  F_tr=location_2_tripletPWOdesign(component_2_location(initial_design[,permutation]+1))
  eigenall=eigen(t(F_tr)%*%F_tr/nrow(F_tr))[[1]]
  eigenall[eigenall<0]=0
  efficiency=exp(sum(log(eigenall))/nrow(F_tr))
  return(efficiency)
}
