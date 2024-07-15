simulation_MCOA_COA<-function(xx,lambda){
  m=length(xx[1,])
  dis_MCOA=0
  D_MCOA=0
  S_MCOA=0
  tPWO_MCOA=0
  tri_MCOA=0

  if(m==17){
    if(lambda>=2){looptime=2}else{looptime=5}
  }else if(m==19){
    if(lambda>=2){looptime=1}else{looptime=3}
  }else if(m==23){
    if(lambda>=2){looptime=1}else{looptime=3}
  }else{looptime=5}

  for(i in 1:looptime){
    print(paste("m=",m," lambda=",lambda,"; loop: ",i,"/",looptime,sep=""))
    x=xx
    m=length(xx[1,])
    if(lambda>1){
      for(i in 1:(lambda-1)){x=rbind(x,xx[,sample(1:m,m)])}
    }
    batch_PWO_efficiency<-function(permutation){permutation_to_PWO_efficiency(initial_design=x,permutation=permutation)}
    batch_SO_efficiency<-function(permutation){permutation_to_SO_efficiency(initial_design=x,permutation=permutation)}
    batch_tPWO_efficiency<-function(permutation){permutation_to_tPWO_efficiency(initial_design=x,permutation=permutation)}
    batch_tri_efficiency<-function(permutation){permutation_to_tri_efficiency(initial_design=x,permutation=permutation)}
    
    n=length(x[,1])
    result<-genetic_on_maximin_mink(x,p=1)
    dmax=max(result[[2]])
    Demax=0;DSOmax=0;DtPWOmax=0;Dtrimax=0
    if(dmax>=dis_MCOA){
      loc=which(result[[2]]==dmax)
      permu=matrix(result[[1]][loc,],ncol=m)
      De=apply(permu,1,batch_PWO_efficiency)
      Demax=max(De)
      permu=matrix(permu[which(De==Demax),],ncol=m)
      if(nrow(permu)>1){
        DSOe=apply(permu,1,batch_SO_efficiency)
        DSOmax=max(DSOe)
        permu=matrix(permu[which(DSOe==DSOmax),],ncol=m)
      }
      if(nrow(permu)>1){
        DtPWOe=apply(permu,1,batch_tPWO_efficiency)
        DtPWOmax=max(DtPWOe)
        permu=matrix(permu[which(DtPWOe==DtPWOmax),],ncol=m)
      }
      if(nrow(permu)>1){
        Dtrie=apply(permu,1,batch_tri_efficiency)
        Dtrimax=max(Dtrie)
        permu=matrix(permu[which(Dtrie==Dtrimax),],ncol=m)
      }
      if(nrow(permu)>1){
        permu=matrix(permu[sample(1:nrow(permu),1),],ncol=m)
      }
      permu=c(permu)
      
      flag=((dmax>(dis_MCOA+10^-5)) | (Demax>(D_MCOA+10^-5)) | (Demax>(D_MCOA-10^-5) & DSOmax>(S_MCOA+10^-5)) |
        (Demax>(D_MCOA-10^-5) & DSOmax>(S_MCOA-10^-5) & (DtPWOmax>tPWO_MCOA+10^-5))| 
          (Demax>(D_MCOA-10^-5) & DSOmax>(S_MCOA-10^-5) & (DtPWOmax>tPWO_MCOA-10^-5) & Dtrimax>(tri_MCOA+10^-5)))
      
      if(flag){
        mcoa_location=component_2_location(x[,permu]+1)
        mcoa_design=x[,permu]
        dis_MCOA=dmax
        D_MCOA=batch_PWO_efficiency(permu)
        S_MCOA=batch_SO_efficiency(permu)
        tPWO_MCOA=batch_tPWO_efficiency(permu)
        tri_MCOA=batch_tri_efficiency(permu)
      }
    }
  }
  output=list(mcoa_design,mcoa_location)
  names(output)<-c("mcoa_design","mcoa_location")
  return(output)
}
