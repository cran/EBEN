\name{EBlassoNEG.Gaussian}
\alias{EBlassoNEG.Gaussian}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The EBlasso Algorithm for Gaussian Model with Normal-Exponential-Gamma (NEG) Prior Distribution
}
\description{
General linear regression, normal-exponential-gamma (NEG) hierarchical prior for regression coefficients
}
\usage{
EBlassoNEG.Gaussian(BASIS, Target, a_gamma, b_gamma, Epis, verbose, group)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BASIS}{sample matrix; rows correspond to samples, columns correspond to features}
  \item{Target}{Response each individual}
  \item{a_gamma}{Hyperparameters control degree of shrinkage; can be obtained via Cross Validation}
  \item{b_gamma}{Hyperparameters control degree of shrinkage; can be obtained via Cross Validation}
    \item{Epis}{TRUE or FALSE for including two-way interactions}
  \item{verbose}{from 0 to 5; larger verbose displays more messages}
  \item{group}{0 or 1; 0: No group effect; 1 two-way interaction grouped. Only valid when Epis = TRUE }
}
\details{If Epis=TURE, the program adds two-way interaction K*(K-1)/2 more columns to BASIS\cr
		for memory efficient, the function pass n_effect to C. n_effect > n_true effects, which is \cr
		a rough guess on how many variables will be selected by the function \cr
		by providing a relative 'small' n_effect, the function will not allocate \cr
		a large trunck of memory during computation.
		}
\value{
	\item{weight}{the none-zero regression coefficients: \cr
	col1,col2 are the indices of the bases(main if equal); \cr
	col3: coefficent value; \cr
	col4: posterior variance;\cr
	col5: t-value;\cr
	col6: p-value}
	\item{WaldScore}{Wald Score}
	\item{Intercept}{Intercept}
	\item{residVar}{residual variance}
	\item{a_gamma}{the hyperparameter; same as input}
	\item{b_gamma}{the hyperparameter; same as input}
}
\references{Cai, X., Huang, A., and Xu, S. (2011). Fast empirical Bayesian LASSO for multiple quantitative trait locus mapping. BMC Bioinformatics 12, 211.}
\author{Anhui Huang; Dept of Electrical and Computer Engineering, Univ of Miami, Coral Gables, FL}


\examples{
library(EBEN)
data(BASIS)
data(y)
n = 50;
k = 100;
BASIS = BASIS[1:n,1:k];
y  = y[1:n];
output = EBlassoNEG.Gaussian(BASIS, y, a_gamma = 0.1, b_gamma = 0.1)

}
