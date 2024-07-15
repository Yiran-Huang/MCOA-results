rm(list=ls())
script_path <- as.character(sys.frame(1)$ofile)
file_name= file.path(dirname(script_path))
source(paste(file_name,'/COA_construct.R',sep=""))
source(paste(file_name,'/algorithm_on_maximin.R',sep=""))
source(paste(file_name,'/functions.R',sep=""))
source(paste(file_name,'/function_MCOA_vs_COA.R',sep=""))

set.seed(1810029)
xx<-COA_construct2(11,1)
lambda=1
output11_1<-simulation_MCOA_COA(xx,lambda)
