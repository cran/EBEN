EBelasticNet.Gaussian <-
function(BASIS,Target,lambda,alpha,Epis = "no",verbose = 0 ){
	N 				= nrow(BASIS);
	K 				= ncol(BASIS);
	if (verbose>0) cat("EBEN Gaussian Model, Epis: ",Epis,"\n");
	if(Epis == "yes"){
		N_effect 		= (K+1)*K/2;
#-----------------------------------------
		Beta 			= rep(0,N_effect *5);
#-----------------------------------------	


#dyn.load("fEBLinearNeFull.dll")

		output<-.C("elasticNetLinearNeEpisEff",
			BASIS 	= as.double(BASIS),
			Target 	= as.double(Target),
			lamda 	= as.double(lambda),
			alpha 	= as.double(alpha),
			Beta 		= as.double(Beta),
			WaldScore 	= as.double(0),
			Intercept 	= as.double(0),
			N 		= as.integer(N),
			K 		= as.integer(K),
			verbose = as.integer(verbose),
			residual = as.double(0),
			PACKAGE="EBEN");
	
#dyn.unload("elasticNetLinearNeFull2.dll")	

	}else {
		N_effect 		= K;
		Beta 			= rep(0,N_effect *4);

#		dyn.load("fEBLinearNeMainEff.so")

		output<-.C("elasticNetLinearNeMainEff",
			BASIS 	= as.double(BASIS),
			Target 	= as.double(Target),
			lamda 	= as.double(lambda),
			alpha 	= as.double(alpha),
			Beta 		= as.double(Beta),
			WaldScore 	= as.double(0),
			Intercept 	= as.double(0),
			N 		= as.integer(N),
			K 		= as.integer(K),
			verbose = as.integer(verbose),
			residual = as.double(0),
			PACKAGE="EBEN");
#		dyn.unload("fEBLinearNeMainEff.so")

	}		
#-------------------------------------------------------------------
if(Epis == "yes"){	
	result 			= matrix(output$Beta,N_effect,5);
	ToKeep 			= which(result[,5]!=0);	
}else
{
	result 			= matrix(output$Beta,N_effect,4);
	ToKeep 			= which(result[,3]!=0);
}
	if(length(ToKeep)==0) { Blup = matrix(0,1,4)
	}else
	{
		nEff 	= length(ToKeep);
		Blup 		= matrix(result[ToKeep,],nEff,4);
	}
	
	fEBresult 			<- list(Blup,output$WaldScore,output$Intercept,output$residual,lambda,alpha);
	rm(list= "output")	
	names(fEBresult)		<-c("weight","WaldScore","Intercept","residVar","lambda","alpha")
	return(fEBresult)
	
}
