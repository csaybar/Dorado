% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/best_method.R
\name{best_method}
\alias{best_method}
\title{Best Method Interpolation}
\usage{
best_method(gauge, cov, idpR = seq(0.8, 3.5, 0.1), formula)
}
\arguments{
\item{gauge}{Is an object of SpatialPointsDataFrame-class.}

\item{cov}{Is an object of RasterLayer.}

\item{idpR}{Is vector numeric of the power coeficient to evaluate.}

\item{formula}{that defines the dependent variable as a linear model
of independent variables; suppose the dependent variable has
name 'z', for Kriging with Extrenal Drift (KED) use the formula
'z~x+y+....', you do not need define}
}
\value{
a List that contains: \code{Interpol} is the KED result in Raster,
 \code{params} being \code{residual} that are spatial residual obtains in the cross-validation (residual),
 \code{MSE} is the Residual Mean squared error and finally \code{var} is  the variogram model.
}
\description{
Best Method Interpolation
}
\examples{
data("Dorado")
gauge <- mean_doble_Station(gauge = Dorado$gauge,cov = Dorado$TRMM)
sat <- Dorado$TRMM
x <- best_method(gauge = gauge,cov = sat,formula = PP_anual~Precipitacion_Anual)
}
\author{
Cesar Aybar <aybar1994@gmail.com>
}
