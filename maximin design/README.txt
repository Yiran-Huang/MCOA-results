SOURCE 'MCOA.R' yields the result of MCOA(110,11).
If you wish to generate another MCOA, you can modify the parameters. For example, if you wish to generate MCOA(2*9*8,9), you should set
----------------------------------------
xx<-COA_construct2(3,2)
lambda=2
----------------------------------------


To reduce the running time, we have decreased the number of loop iterations by a factor of 10 compared to the original program. You can modify the parameter 'looptime' within the 'function_MCOA_vs_COA.R' file.