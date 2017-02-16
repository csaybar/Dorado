#' Best Method Interpolation
#' @author Cesar Aybar <aybar1994@gmail.com>
#' @param gauge Is an object of SpatialPointsDataFrame-class.
#' @param cov Is an object of RasterLayer.
#' @param formula that defines the dependent variable as a linear model
#' of independent variables; suppose the dependent variable has
#' name 'z', for Kriging with Extrenal Drift (KED) use the formula
#' 'z~x+y+....', you do not need define
#' @param idpR Is vector numeric of the power coeficient to evaluate.
#' @importFrom automap autofitVariogram
#' @importFrom raster extract projection writeRaster
#' @importFrom sp coordinates
#' @importFrom gstat krige.cv idw krige
#' @importFrom stats as.formula
#' @return a List that contains: \code{Interpol} is the KED result in Raster,
#'  \code{params} being \code{residual} that are spatial residual obtains in the cross-validation (residual),
#'  \code{MSE} is the Residual Mean squared error and finally \code{var} is  the variogram model.
#'  @export

best_method <- function(gauge, cov, idpR = seq(0.8, 3.5, 0.1), formula) {

  # Proccesing Data -------------------------------------------- Extract
  # Data
  ext <- raster::extract(cov, gauge, cellnumber = F, sp = T)
  station <- gauge
  linear <- na.omit(ext@data) %>% tbl_df %>% mutate_all(as.character) %>%
    mutate_all(as.numeric)
  llm <- lm(formula, linear)
  station$residuals <- llm$residuals


  # Define Grid
  # -------------------------------------------------------------
  point <- rasterToPoints(cov) %>% data.frame
  coordinates(point) <- ~x + y
  projection(point) <- projection(cov)

  # RIDW
  # --------------------------------------------------------------------
  idpRange <- idpR
  mse_RIDW <- rep(NA, length(idpRange))
  message("Analyzing RIDW method wait ....")
  for (i in 1:length(idpRange)) {
    mse_RIDW[i] <- mean(krige.cv(residuals ~ 1, station, nfold = nrow(station),
                                 nmax = Inf, set = list(idp = idpRange[i]), verbose = F)$residual^2)
  }
  poss_RIDW <- which(mse_RIDW %in% min(mse_RIDW))
  best_mse_RIDW <- mse_RIDW[poss_RIDW]
  idp_RIDW <- idpRange[poss_RIDW]

  # IDW
  # --------------------------------------------------------------------

  namesF <- unlist(strsplit(as.character(formula), " "))
  max_k <- floor(length(namesF)/2) + 1
  name_cov <- namesF[!namesF %in% c("~", "+", "-", "*", "/")][1]
  name_cov_lm <- namesF[!namesF %in% c("~", "+", "-", "*", "/")][2:max_k]
  form <- paste0(name_cov, "~1")
  mse_IDW <- rep(NA, length(idpRange))
  message("Analyzing IDW method wait ....")
  for (i in 1:length(idpRange)) {
    mse_IDW[i] <- mean(krige.cv(as.formula(form), station, nfold = nrow(station),
                                nmax = Inf, set = list(idp = idpRange[i]), verbose = F)$residual^2)
  }

  poss_IDW <- which(mse_IDW %in% min(mse_IDW))
  best_mse_IDW <- mse_IDW[poss_IDW]
  idp_IDW <- idpRange[poss_IDW]

  # KO
  # --------------------------------------------------------------------

  vm_fit <- variogram_fit(gauge, cov, as.formula(form))
  message("Analyzing KO method wait ....")
  Zs.cv_KO <- krige.cv(as.formula(form), ext, vm_fit$ftvariogram$var_model,
                       nfold = nrow(ext), nmax = Inf)
  best_mse_KO <- mean(Zs.cv_KO$residual^2)

  # KED
  # --------------------------------------------------------------------

  vm.fit <- variogram_fit(gauge, cov, formula)
  ext <- vm.fit$ext
  vm.fit <- vm.fit$ftvariogram
  message("Analyzing KED method wait ....")
  Zs.cv_KED <- krige.cv(formula, ext, vm.fit$ftvariogram$var_model, nfold = nrow(ext),
                        nmax = Inf)
  best_mse_KED <- mean(Zs.cv_KED$residual^2)

  # LM
  # --------------------------------------------------------------------
  message("Analyzing LM method wait ....")
  cov_lm <- cov[[name_cov_lm]]
  OBSp <- sum(stack(mapply(function(i) cov_lm[[i]] * llm$coefficients[i +
                                                                        1], 1:nlayers(cov_lm)))) + llm$coefficients[1]
  best_mse_LM <- mean((extract(OBSp, gauge) - gauge@data[, 1])^2)

  # BIAS
  # --------------------------------------------------------------------
  message("Analyzing BIAS method wait ....")
  biasC <- sum(gauge@data[, 1])/sum(extract(cov_lm[[1]], gauge))
  bias_map <- cov_lm[[1]] * biasC
  best_mse_bias <- mean((extract(bias_map, gauge) - gauge@data[, 1])^2)


  listT <- list(ked_mse = best_mse_KED, idw = list(idw_mse = best_mse_IDW,
                                                   idp = idp_IDW),
                ridw = list(ridw_mse = best_mse_RIDW, idp = idp_RIDW),
                ko = best_mse_KO, lm = best_mse_LM, bias = best_mse_bias)

  mse_T <- c(listT$ked_mse, listT$idw$idw_mse, listT$ridw$ridw_mse, listT$ko,
             listT$lm, listT$bias)
  minT <- min(mse_T)
  position_method <- which(mse_T == minT)
  methods <- c("ked", "idw", "ridw", "ko", "lm", "bias")
  if (position_method == 1)
    print("El mejor metodo es KED")
  if (position_method == 2)
    print("El mejor metodo es IDW")
  if (position_method == 3)
    print("El mejor metodo es RIDW")
  if (position_method == 4)
    print("El mejor metodo es KO")
  if (position_method == 5)
    print("El mejor metodo es LM")
  if (position_method == 6)
    print("El mejor metodo es BIAS")
  listT$bestMethod <- methods[position_method]
  return(listT)
}
