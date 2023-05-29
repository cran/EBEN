\name{EBEN-internal}
\title{Internal EBEN functions}
\alias{ijIndex}
\alias{CVonePair}
\alias{lambdaMax}
\description{Internal EBEN functions}
\usage{
	ijIndex(trueLoc,K)
	CVonePair(X,y,nFolds,foldId,hyperpara,Epis,prior,family,verbose,group)
	lambdaMax(X,y,Epis)
}
\author{Anhui Huang and Dianting Liu \cr Dept of Electrical and Computer Engineering, Univ of Miami, Coral Gables, FL}

\details{These are not intended for use by users.
	\code{ijIndex}Function for looking at the pair of interaction terms.
	\code{CVonePair}Function performs nFolds CV for the given one pair of hyperparameter.
	\code{lambdaMax}Function calculate the maximum lambda for EBlasso-NE and EBEN in CV.
}
\keyword{internal}
