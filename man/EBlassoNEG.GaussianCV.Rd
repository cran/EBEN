\name{EBlassoNEG.GaussianCV}
\alias{EBlassoNEG.GaussianCV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Cross Validation (CV) Function to Determine Hyperparameters of
		the EBlasso Algorithm for Gaussian Model with Normal-Exponential-Gamma (NEG) Prior Distribution}
\description{Hyperparameters control degree of shrinkage, and are obtained
				via Cross Validation. This program performs three steps of CV.\cr
			1st: a = b = 0.001, 0.01, 0.1, 1;\cr
			2nd: fix b= b1; a=[-0.5, -0.4, -0.3, -0.2, -0.1, -0.01, 0.01, 0.05, 0.1, 0.5, 1];\cr
			3rd: fix a = a2; b= 0.01 to 10 with a step size of one for b > 1 and a step size of one on the logarithmic scale for b < 1\cr
			In the 2nd step, a can take value from -1 and values in [-1, -0.5] can be added to the set in line 13 of this function (The smaller a is, the less shrinkage.)
			}
\usage{
EBlassoNEG.GaussianCV(BASIS, Target, nFolds, foldId, Epis,verbose, group)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BASIS}{sample matrix; rows correspond to samples, columns correspond to features}
  \item{Target}{Class label of each individual, TAKES VALUES OF 0 OR 1}
  \item{nFolds}{number of n-fold cv}
  \item{foldId}{random assign samples to different folds}
  \item{Epis}{TRUE or FALSE for including two-way interactions}
  \item{verbose}{from 0 to 5; larger verbose displays more messages}
  \item{group}{TRUE or FALSE; FALSE: No group effect; TRUE two-way interaction grouped. Only valid when Epis = TRUE }
}
\details{If Epis= TRUE, the program adds two-way interaction K*(K-1)/2 more columns to BASIS \cr
		Note: Given the fact that degree of shrinkage is a monotonic function of (a,b), \cr
		The function implemented a 3-step search as described in Huang, A. 2014, for full \cr
		grid search, user needs to modify the function accordingly.}
\value{
	\item{CrossValidation}{col1: hyperparameters; col2: loglikelihood mean; standard ERROR of nfold mean log likelihood}
	\item{a_optimal}{the optimal hyperparameter as computed}
	\item{b_optimal}{the optimal hyperparameter as computed}
}

\references{Huang A, Xu S, Cai X: Empirical Bayesian LASSO-logistic regression for multiple binary trait locus mapping. BMC genetics 2013, 14(1):5. \cr
Huang, A., S. Xu, et al. Whole-genome quantitative trait locus mapping reveals major role of epistasis on yield of rice. PLoS ONE 2014, 9(1): e87330.}
\author{Anhui Huang; Dept of Electrical and Computer Engineering, Univ of Miami, Coral Gables, FL}



\examples{
library(EBEN)
data(BASIS)
data(y)
#reduce sample size to speed up the running time
n = 50;
k = 100;
BASIS = BASIS[1:n,1:k];
y  = y[1:n];
\dontrun{
CV = EBlassoNEG.GaussianCV(BASIS, y, nFolds = 3,Epis = FALSE)
}
}

