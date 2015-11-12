EBelasticNet.BinomialCV <- function(BASIS,Target,nFolds,Epis="no",foldId = 0)
{
  
  cat("EB-Elastic Net Linear Model, Epis: ",Epis, ";", nFolds, "fold cross-validation\n");
  N 					= nrow(BASIS);
  K 					= ncol(BASIS);
  #set.seed(proc.time())
  if(length(foldId)!=N)
  {
    if(N%%nFolds!=0){
      foldId 			= sample(c(rep(1:nFolds,floor(N/nFolds)),1:(N%%nFolds)),N);
    }else{
      foldId 			= sample(rep(1:nFolds,floor(N/nFolds)),N);
    }
  }
  lambda_Max			= log(1.1);
  for(i_b in 1:K){
    basis 			= BASIS[,i_b];
    basis 			= basis/sum(basis*basis);
    corBy 			= basis%*%(Target-mean(Target));
    if(corBy>lambda_Max) lambda_Max = corBy;
  }
  
  if(Epis == "yes"){
    for(i_b in 1:(K-1)){
      for(i_bj in (i_b + 1):K){
        basis 	= BASIS[,i_b]*BASIS[,i_bj];
        basis 	= basis/sum(basis*basis);
        corBy 	= basis%*%(Target-mean(Target));
        if(corBy>lambda_Max) lambda_Max = corBy;
      }
    }		
  }
  lambda_Max 		= lambda_Max*10;
  Lambda 				= c(1, 0.5, 0.25, 0.1, 0.05, 0.01)*lambda_Max;
  Lambda = Lambda[Lambda > 0.1]
  N_step 				= length(Lambda);
  
  step 				= 1;
  Alpha 				= seq(from = 1, to = 0.05, by = -0.2)
  nAlpha 				= length(Alpha);
  
  Likelihood 			= mat.or.vec((N_step*nAlpha),4);
  logL 				= mat.or.vec(nFolds,1);
  for(i_alpha in 1:nAlpha){
    alpha 			= Alpha[i_alpha];
    cat("Testing alpha", i_alpha, "/",nAlpha,":\t\talpha: ",alpha,"\n")
    for (i_s in 1:N_step){
      
      lambda 		= Lambda[i_s];
      
      cat("\tTesting step", step, "\t\tlambda: ",lambda,"\n")
      for(i in 1:nFolds){
        index  			= which(foldId!=i);
        Basis.Train 	= BASIS[index,];
        Target.Train 	= Target[index];
        index  			= which(foldId == i);
        Basis.Test  	= BASIS[index,];
        Target.Test 	= Target[index];
        SimF2fEB 		<-EBelasticNet.Binomial(Basis.Train,Target.Train,lambda,alpha,Epis);
        M				= length(SimF2fEB$weight)/6;
        Betas 			<- matrix(SimF2fEB$weight,nrow= M,ncol =6, byrow= FALSE);
        Mu  			= Betas[,3];
        Mu0 			= SimF2fEB$Intercept[1];
        
        rm(list="SimF2fEB");
        ntest 			= nrow(Basis.Test);
        #M 		= nrow(Betas);
        if(M==1 && Betas[1,1]== 0)
        {
          logL[i]  	= 0;
        }else
        {
          basisTest 	= matrix(rep(0,ntest*M),ntest,M);
          for(i_basis in 1:M){
            loc1 = Betas[i_basis,1];
            loc2 = Betas[i_basis,2];
            if(loc1==loc2){ 	basisTest[,i_basis] =  Basis.Test[,loc1];}
            else{			basisTest[,i_basis] =  Basis.Test[,loc1]* Basis.Test[,loc2];}
          }
          
          temp 		= exp(Mu0 + basisTest%*%Mu);
          logL[i] 	= mean(Target.Test*log(temp/(1+temp)) + (1-Target.Test)*log(1/(1+temp)));
        }
        
      }
      Likelihood[step,] = c(alpha, lambda,mean(logL),sd(logL));
      step 			= step + 1;
    }
    
  }
  index 				= which.max(Likelihood[,3]);
  Res.lambda			= Likelihood[index,2];
  Res.alpha 			= Likelihood[index,1];
  result 				<- list(Likelihood,Res.alpha,Res.lambda);
  names(result)		<-c("CrossValidation","Alpha_optimal","Lambda_optimal");
  return(result);
}