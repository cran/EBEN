\name{EBlassoNEG.Binomial}
\alias{EBlassoNEG.Binomial}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{The EBlasso Algorithm for Binomial Model with Normal-Exponential-Gamma (NEG) Prior Distribution}
\description{Generalized linear regression, normal-exponential-gamma (NEG) hierarchical prior for regression coefficients}
\usage{
EBlassoNEG.Binomial(BASIS, Target, a_gamma, b_gamma, Epis,verbose,group)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BASIS}{sample matrix; rows correspond to samples, columns correspond to features}
  \item{Target}{Class label of each individual, TAKES VALUES OF 0 OR 1}
  \item{a_gamma}{Hyperparameters control degree of shrinkage; can be obtained via Cross Validation; a_gamma>=-1}
  \item{b_gamma}{Hyperparameters control degree of shrinkage; can be obtained via Cross Validation; b_gamma>0}
  \item{Epis}{TRUE or FALSE for including two-way interactions}
  \item{verbose}{0 or 1; 1: display message; 0 no message}
  \item{group}{0 or 1; 0: No group effect; 1 two-way interaction grouped. Only valid when Epis = TRUE }
}
\details{If Epis=TRUE, the program adds two-way interaction K*(K-1)/2 more columns to BASIS}
\value{
	\item{weight}{the none-zero regression coefficients: \cr
	col1,col2 are the indices of the bases(main if equal); \cr
	col3: coefficent value; \cr
	col4: posterior variance;\cr
	col5: t-value;\cr
	col6: p-value}
	\item{logLikelihood}{log likelihood with the final regression coefficients}
	\item{WaldScore}{Wald Score}
	\item{Intercept}{Intercept}
	\item{a_gamma}{the hyperparameter; same as input}
	\item{b_gamma}{the hyperparameter; same as input}
}
\references{Huang, A., Xu, S., and Cai, X.(2012). Empirical Bayesian LASSO-logistic regression for multiple binary trait locus mapping. BMC Genetics. Submitted}
\author{Anhui Huang; Dept of Electrical and Computer Engineering, Univ of Miami, Coral Gables, FL}



%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
library(EBEN)
data(BASISbinomial)
data(yBinomial)
#reduce sample size to speed up the running time
n = 50;
k = 100;
BASIS = BASISbinomial[1:n,1:k];
y  = yBinomial[1:n];
output = EBlassoNEG.Binomial(BASIS,y,0.1,0.1,Epis = FALSE)
}

