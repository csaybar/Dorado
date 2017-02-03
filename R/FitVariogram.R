#' Fit Varigoram using autofitVariogram modified
#' @author Paul Hiemstra
#' @description This function is a very slight variation about the function autofitVariogram of packages
#' automap. Only added the param boundaries so the user can define it.
#' @details   Geostatistical routines are used from package \code{gstat}.
#'  A few simple choices are made when estimating the inital guess for \code{fit.variogram}.
#'  The initial sill is estimated as the \code{mean} of the \code{max} and the \code{median}
#'  of the semi-variance. The inital range is defined as 0.10 times the diagonal of the bounding
#'  box of the data. The initial nugget is defined as the \code{min} of the the semi-variance.
#'  There are five different types of models that are often used:
#'  \describe{
#'    \item{Sph}{A shperical model.}
#'    \item{Exp}{An exponential model.}
#'    \item{Gau}{A gaussian model.}
#'    \item{Mat}{A model of the Matern familiy}
#'    \item{Ste}{Matern, M. Stein's parameterization}
#'    }
#'A list of all permitted variogram models is available by typing vgm() into the R console.
#'      \code{autofitVariogram} iterates over the variogram models listed in \code{model} and picks the model
#'      that has the smallest residual sum of squares with the sample variogram. For the Matern model, all the
#'      kappa values in \code{kappa} are tested.
#'
#'      Note that when using the power model, and not specifying starting values yourself, the sill is set to 1,
#'      the range to 1 and the nugget to 0. This is because the normal initial values for those paramters don't
#'      work well with the power model. I consider this a temporary solution, any suggestions are appreciated.
#'
#' It is possible to pass anisotropy parameters to \code{autofitVariogram}. However, \code{autofitVariogram} does not fit anisotropic variogram models. The function sees the anisotropic sample variogram as one big sample variogram. So it fits an average isotropic variogram model from the anisotropic sample variogram. A warning is issued when a users passes \code{alpha} to \code{autofitVariogram}.
#' @return  An object of type \code{autofitVariogram} is returned. This object contains the experimental variogram,
#'  the fitted variogram model and the sums of squares (\code{sserr}) between the sample variogram and the
#'  fitted variogram model.
#' @param formula that defines the dependent variable as a linear model
#' of independent variables; suppose the dependent variable has
#' name 'z', for ordinary and simple kriging use the formula
#' 'z~1'; for simple kriging also define 'beta' (see below); for
#' universal kriging, suppose 'z' is linearly dependent on 'x'
#' and 'y', use the formula 'z~x+y'.
#' @param input_data An object of \link[sp]{SpatialPointsDataFrame-class}.
#' @param model The list of variogrammodels that will be tested.
#' @param kappa Smoothing parameter of the Matern model. Provide a list if you want to check
#' more than one value.
#' @param boundaries Are the Boundaries for the bins in km.
#' @param fix.values Can be used to fix a variogram parameter to a certain value. It
#' consists of a list with a length of three. The items describe the
#' fixed value for the nugget, range and sill respectively. They need to be given in that order.
#' Setting the value to NA means that the value is not fixed.
#' @param verbose logical, if TRUE the function will give extra feedback on the fitting process
#' @param GLS.model If a variogram model is passed on through this parameter a Generalized Least Squares
#' sample variogram is calculated.
#' @param start_vals Can be used to give the starting values for the variogram fitting. The items describe the
#' fixed value for the nugget, range and sill respectively. They need to be given in that order.
#' Setting the value to NA means that the value will be automatically chosen.
#' @param ... parameters that are passed on to \link[gstat]{variogram} variogram when calculating the sample variogram
#' @param miscFitOptions A list with named arguments that provide additional control over the fitting process.
#' For example: \code{list(merge.small.bins = TRUE)}. If the list is empty, autofitVariogram
#' uses default values. The following parameters can be set:
#' \describe{
#'  \item{\code{merge.small.bins}:}{logical, when TRUE, the function checks if there are bins with less than 5 points.
#'  If so, the first two bins are merged and the check is repeated. This is done until all bins have more
#'  than \code{min.np.bin} points.}
#'  \item{\code{min.np.bin}:}{integer, the minimum number of points allowed in a bin before
#'  we start merging bins. See also \code{merge.small.bins}.}
#'  }
#' @seealso \code{\link[gstat]{fit.variogram}}, \code{\link{autoKrige}}, \code{\link{posPredictionInterval}}
#' @importFrom sp spDists bbox is.projected
#' @importFrom gstat variogram gstat fit.variogram vgm
#' @import gstat
#' @export


FitVariogram = function(formula, input_data, model = c("Sph", "Exp", "Gau", "Ste"),
                            kappa = c(0.05, seq(0.2, 2, 0.1), 5, 10), fix.values = c(NA,NA,NA),
                            verbose = FALSE, GLS.model = NA, start_vals = c(NA,NA,NA),
                            miscFitOptions = list(),boundaries,...)
  # This function automatically fits a variogram to input_data
{
  # Check for anisotropy parameters
  if('alpha' %in% names(list(...)))
    warning('Anisotropic variogram model fitting not supported, see the documentation of autofitVariogram for more details.')
  # Take the misc fit options and overwrite the defaults by the user specified ones
  miscFitOptionsDefaults = list(merge.small.bins = TRUE, min.np.bin = 5)
  miscFitOptions = modifyList(miscFitOptionsDefaults, miscFitOptions)

  # Create boundaries
  longlat = !is.projected(input_data)
  if(is.na(longlat)) longlat = FALSE
  diagonal = spDists(t(bbox(input_data)), longlat = longlat)[1,2]                # 0.35 times the length of the central axis through the area
  if(!missing(boundaries))  boundaries =  boundaries         # Boundaries for the bins in km
  else{
    boundaries = c(2, 4, 6, 9, 12, 15, 25, 35, 50, 65, 80, 100) *
    diagonal * 0.35/100
  }
  # If you specifiy a variogram model in GLS.model the Generelised least squares sample variogram is constructed
  if(!is(GLS.model, "variogramModel")) {
    experimental_variogram = variogram(formula, input_data,boundaries = boundaries, ...)
  } else {
    if(verbose) cat("Calculating GLS sample variogram\n")
    g = gstat(NULL, "bla", formula, input_data, model = GLS.model, set = list(gls=1))
    experimental_variogram = variogram(g, boundaries = boundaries, ...)
  }

  # request by Jon Skoien
  if(miscFitOptions[["merge.small.bins"]]) {
    if(verbose) cat("Checking if any bins have less than 5 points, merging bins when necessary...\n\n")
    while(TRUE) {
      if(length(experimental_variogram$np[experimental_variogram$np < miscFitOptions[["min.np.bin"]]]) == 0 | length(boundaries) == 1) break
      boundaries = boundaries[2:length(boundaries)]
      if(!is(GLS.model, "variogramModel")) {
        experimental_variogram = variogram(formula, input_data,boundaries = boundaries, ...)
      } else {
        experimental_variogram = variogram(g, boundaries = boundaries, ...)
      }
    }
  }

  # set initial values
  if(is.na(start_vals[1])) {  # Nugget
    initial_nugget = min(experimental_variogram$gamma)
  } else {
    initial_nugget = start_vals[1]
  }
  if(is.na(start_vals[2])) { # Range
    initial_range = 0.1 * diagonal   # 0.10 times the length of the central axis through the area
  } else {
    initial_range = start_vals[2]
  }
  if(is.na(start_vals[3])) { # Sill
    initial_sill = mean(c(max(experimental_variogram$gamma), median(experimental_variogram$gamma)))
  } else {
    initial_sill = start_vals[3]
  }

  # Determine what should be automatically fitted and what should be fixed
  # Nugget
  if(!is.na(fix.values[1]))
  {
    fit_nugget = FALSE
    initial_nugget = fix.values[1]
  } else
    fit_nugget = TRUE

  # Range
  if(!is.na(fix.values[2]))
  {
    fit_range = FALSE
    initial_range = fix.values[2]
  } else
    fit_range = TRUE

  # Partial sill
  if(!is.na(fix.values[3]))
  {
    fit_sill = FALSE
    initial_sill = fix.values[3]
  } else
    fit_sill = TRUE

  getModel = function(psill, model, range, kappa, nugget, fit_range, fit_sill, fit_nugget, verbose)
  {
    if(verbose) debug.level = 1 else debug.level = 0
    if(model == "Pow") {
      warning("Using the power model is at your own risk, read the docs of autofitVariogram for more details.")
      if(is.na(start_vals[1])) nugget = 0
      if(is.na(start_vals[2])) range = 1    # If a power mode, range == 1 is a better start value
      if(is.na(start_vals[3])) sill = 1
    }
    obj = try(fit.variogram(experimental_variogram,
                            model = vgm(psill=psill, model=model, range=range,
                                        nugget=nugget,kappa = kappa),
                            fit.ranges = c(fit_range), fit.sills = c(fit_nugget, fit_sill),
                            debug.level = 0),
              TRUE)
    if("try-error" %in% class(obj)) {
      #print(traceback())
      warning("An error has occured during variogram fitting. Used:\n",
              "\tnugget:\t", nugget,
              "\n\tmodel:\t", model,
              "\n\tpsill:\t", psill,
              "\n\trange:\t", range,
              "\n\tkappa:\t",ifelse(kappa == 0, NA, kappa),
              "\n  as initial guess. This particular variogram fit is not taken into account. \nGstat error:\n", obj)
      return(NULL)
    } else return(obj)
  }


  # Automatically testing different models, the one with the smallest sums-of-squares is chosen
  test_models = model
  SSerr_list = c()
  vgm_list = list()
  counter = 1

  for(m in test_models) {
    if(m != "Mat" && m != "Ste") {        # If not Matern and not Stein
      model_fit = getModel(initial_sill - initial_nugget, m, initial_range, kappa = 0, initial_nugget, fit_range, fit_sill, fit_nugget, verbose = verbose)
      if(!is.null(model_fit)) {	# skip models that failed
        vgm_list[[counter]] = model_fit
        SSerr_list = c(SSerr_list, attr(model_fit, "SSErr"))}
      counter = counter + 1
    } else {                 # Else loop also over kappa values
      for(k in kappa) {
        model_fit = getModel(initial_sill - initial_nugget, m, initial_range, k, initial_nugget, fit_range, fit_sill, fit_nugget, verbose = verbose)
        if(!is.null(model_fit)) {
          vgm_list[[counter]] = model_fit
          SSerr_list = c(SSerr_list, attr(model_fit, "SSErr"))}
        counter = counter + 1
      }
    }
  }

  # Check for negative values in sill or range coming from fit.variogram
  # and NULL values in vgm_list, and remove those with a warning
  strange_entries = sapply(vgm_list, function(v) any(c(v$psill, v$range) < 0) | is.null(v))
  if(any(strange_entries)) {
    if(verbose) {
      print(vgm_list[strange_entries])
      cat("^^^ ABOVE MODELS WERE REMOVED ^^^\n\n")
    }
    warning("Some models where removed for being either NULL or having a negative sill/range/nugget, \n\tset verbose == TRUE for more information")
    SSerr_list = SSerr_list[!strange_entries]
    vgm_list = vgm_list[!strange_entries]
  }

  if(verbose) {
    cat("Selected:\n")
    print(vgm_list[[which.min(SSerr_list)]])
    cat("\nTested models, best first:\n")
    tested = data.frame("Tested models" = sapply(vgm_list, function(x) as.character(x[2,1])),
                        kappa = sapply(vgm_list, function(x) as.character(x[2,4])),
                        "SSerror" = SSerr_list)
    tested = tested[order(tested$SSerror),]
    print(tested)
  }

  result = list(exp_var = experimental_variogram, var_model = vgm_list[[which.min(SSerr_list)]], sserr = min(SSerr_list))
  class(result) = c("autofitVariogram","list")

  return(result)
}
