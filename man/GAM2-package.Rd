\name{GAM2-package}
\alias{GAM2-package}
\alias{GAM2}
\docType{package}
\title{
  Fast GAM Implementation
}
\description{
   The package contains General Additive Modeling (GAM) functionality by using the smoothing splines method.

}
\details{
SmoothingSpline function contains the formula of the form y ~ f1(X1)+f2(X2)+...+fp(Xp), where f can be linear, categoric or splines. Its input includes train_data X: Attribute matrix (n*p), and y: the response vector; lambda: a vector of smoothing parameters; max_iter: maximum iteration to train the model; tol: the convergence tolerence.
}
\author{
Haodong Liang & Ziqi Wei

Maintainer: Ziqi Wei <ziqiwei@umich.edu>
}
\references{
  An Introduction to Statistical Learning
}
\keyword{ package }
\examples{
  \dontrun{
     #' y <- mpg
     #' X <- cbind(cyl, hp, drat)
     #' Coef <- SmoothingSpline(X, y, lambda=c(100,100,100))
  }
}
