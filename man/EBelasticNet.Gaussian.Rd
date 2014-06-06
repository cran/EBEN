\name{EBelasticNet.Gaussian}
\alias{EBelasticNet.Gaussian}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{The EB Elastic Net Algorithm for Gaussian Model}
\description{General linear regression, normal-Gamma (NG) hierarchical prior for regression coefficients}
\usage{EBelasticNet.Gaussian(BASIS, Target, lambda, alpha,Epis = "no",verbose = 0)}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BASIS}{sample matrix; rows correspond to samples, columns correspond to features}
  \item{Target}{Response each individual}
  \item{lambda}{Hyperparameter controls degree of shrinkage; can be obtained via Cross Validation; lambda>0}
  \item{alpha}{Hyperparameter controls degree of shrinkage; can be obtained via Cross Validation; 0<alpha<1}
  \item{Epis}{"yes" or "no" for including two-way interactions}
  \item{verbose}{0 or 1; 1: display message; 0 no message}
}
\details{If Epis="yes", the program adds two-way interaction of K*(K-1)/2 more columns to BASIS}
\value{
	\item{weight}{the none-zero regression coefficients: col1,col2 are the indices of the bases(main if equal); col3: coefficent value; col4: variance}
	\item{WaldScore}{Wald Score}
	\item{Intercept}{Intercept}
	\item{lambda}{the hyperparameter; same as input lambda}
	\item{alpha}{the hyperparameter; same as input alpha}	
}
\references{Huang A, Xu S, Cai X (2013): Empirical Bayesian elastic net for multiple quantitative trait locus mapping. submitted.}
\author{Anhui Huang; Dept of Electrical and Computer Engineering, Univ of Miami, Coral Gables, FL}


%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
library(EBEN)
data(BASIS)
data(y)
#Blup = EBelasticNet.Gaussian(BASIS, y,lambda = 0.0072,alpha = 0.95, Epis = "no",verbose = 0)
##The above codes were designed for the example present in: 
#Huang, A., Xu, S., and Cai, X. (2014). Empirical Bayesian elastic net
# for multiple quantitative trait locus mapping. Heredity. To appear.
#To meet with the '5s' policy for CRAN submission, the following toy data was used.
N  	= 50;
K 	= 100; 
BASIS  = BASIS[1:N,1:K]
y 		= y[1:N]
Blup = EBelasticNet.Gaussian(BASIS, y,lambda = 0.0072,alpha = 0.95, Epis = "no",verbose = 0)



# t-test
intercept 		= Blup$Intercept;
M				= length(Blup$weight)/4;
Betas 			<- matrix(Blup$weight,nrow= M,ncol =4, byrow= FALSE);
N 				= nrow(BASIS)
Betas 			= rbind(c(0,0,intercept,0),Betas);
t 				= abs(Betas[,3])/(sqrt(Betas[,4])+ 1e-20);
pvalue 			= 1- pt(t,df=(N-1));
betas 			= cbind(Betas,t,pvalue);
ind1 			= which(betas[,6]>0.05);
ind2 			= which(betas[,6]<=0.05);
betas 			= rbind(betas[ind2,,drop= FALSE],betas[ind1,,drop= FALSE]);

}
