computePolarizationMeasure <- function (
  profiles=NULL,
  measure,
  lambda=NULL,
  alpha=NULL,
  beta=NULL,
  error=10^(-3),
  zeta=getzeta(profiles),
  weights = {
    if(!is.null(profiles))
      profiles$freq
    else
      rep(1, nrow(zeta))
  },
  linext=lingen(zeta),
  nit = floor({n <- nrow(zeta); n^5*log(n)+n^4*log(error^(-1))}),
  maxint = 2^31-1
)

{

  start_time <- Sys.time()


  valid_measures <- c("abulnaga", "apouey", "blairlacy1", "blairlacy2", "kobusmilos", "kvalseth",
                     "leik", "leti", "reardon1", "reardon3", "reardon4", "reardon2", "berrymielke")

  aliases <- c("reardon2", "berrymielke")

  if (!(measure %in% valid_measures))
    stop(paste0("measure must be one of: ", paste(valid_measures[order(valid_measures)], collapse=", ")))


  if (measure=="abulnaga") {
    if (is.null(alpha) | is.null(beta)) stop("You must specify both alpha and beta for measure abulnaga")
    else if (alpha>1 | beta>1) stop("alpha and beta must be <= 1 for measure abulnaga")
  }

  else if (measure=="apouey") {
    if (is.null(alpha)) stop("You must specify alpha for measure apouey")
    else if (!(0 < alpha && alpha < 1)) stop("alpha must be 0 < alpha < 1 for measure apouey")
  }

  else if (measure=="kobusmilos") {
    if (is.null(alpha) | is.null(beta)) stop("You must specify both alpha and beta for measure kobusmilos")
    else if (alpha<0 | beta <0) stop("alpha and beta must be >= 0 for measure kobusmilos")
  }


  else if (measure =="berrymielke") measure <- "leti"
  else if (measure == "reardon2") measure <- "blairlacy1"


  if (measure == "leti") measure_real <- "blairlacy1"
  else measure_real <- measure

  real_measures <- c("abulnaga", "apouey", "blairlacy1", "blairlacy2", "kobusmilos", "kvalseth",
                     "leik", "reardon1", "reardon3", "reardon4")

  real_measure_index <- which(real_measures == measure_real) - 1

  if (!is.null(lambda)) {
    reduceLE <- function(x,y) {
      ord <- rownames(x)
      return(x * y[ord, ord])
    }
    lst <- LE(lambda)
    lstZeta <- LE2incidence(lst, varlen=apply(profiles$profiles, 2, function(x) length(unique(x))))
    zeta <- Reduce(reduceLE, lstZeta)
    class(zeta) <- "incidence"
  }

  if (is.null(alpha)) alpha <- 0
  if (is.null(beta)) beta <- 0

  n <- nrow(zeta)
  ntot <- sum(weights)

  nitot <- nit
  nit <- rep(maxint, nitot %/% maxint)
  remainder <- nitot %% maxint

  if (remainder > 0) nit <- c(nit, remainder)

  pb <- txtProgressBar(style = 3, min = 0, max = nitot)
  cont <- 0

  ## initialize arrays to calculate quantiles
  p <- c(0.1, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.9)
  quantiles_n <- 1:25

  quantiles_f <- rep(0.0,25)
  quantiles_f[25] <- 1
  quantiles_f[seq(3,24,2)] <- p
  quantiles_f[seq(2,24,2)] <- (quantiles_f[seq(1,23,2)] + quantiles_f[seq(3,25,2)]) / 2


  quantiles_d <- 1 + 2*(length(p) + 1)*quantiles_f
  quantiles_q <- rep(0.0, 25)

  values_for_bd <- sapply(1:n, function(i) i*(n-i))
  ##

  l <- list(
    zeta = zeta,
    measure = real_measure_index,
    alpha = alpha,
    beta = beta,
    linext = linext,
    n = n,
    nit = 0,
    relfreqs = weights/ntot,
    means = numeric(),
    variances = numeric(),
    rolling_mean = numeric(),
    polarization_mean = 0.0,
    polarization_variance = 0.0,
    quantiles_q = quantiles_q,
    quantiles_n = quantiles_n,
    quantiles_f = quantiles_f,
    quantiles_d = quantiles_d,
    polar_min = 0,
    polar_max = 0,
    firstit = 1,
    probvec_bd = values_for_bd
  )
  class(l) <- "pre_parsec_polarization"

  for(j in nit) {
    l$nit <- j
    l <- runC(l)
    cont <- cont + j
    setTxtProgressBar(pb, cont)

  }

  close(pb)


  measures_with_scaling <- head(valid_measures, -length(aliases))

  scaling_index <- which(measures_with_scaling == measure)

  scaling_factors <- c(1, -(2^alpha/(n-1)), 4/(n-1), -sqrt(4/(n-1)), 2, -1, 2/(n-1),
                       if (ntot %% 2 == 0) 4/(n-1) else 4*ntot^2/((n-1)*(ntot^2-1)),
                       -1/(n-1), 2/(n-1), -1/(n-1))

  sum_factors <- c(0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1)

  l$polarization_mean <- l$polarization_mean/nitot*scaling_factors[scaling_index] + sum_factors[scaling_index]

  quantiles <- l$quantiles_q[seq(3,24,2)]*scaling_factors[scaling_index] + sum_factors[scaling_index]
  minpol <- l$polar_min*scaling_factors[scaling_index] + sum_factors[scaling_index]
  maxpol <- l$polar_max*scaling_factors[scaling_index] + sum_factors[scaling_index]

  regularize_variance <- function(list, nit, nitot, scaling_factors, scaling_index) {
    rolling_variances <- c(list$variances[1])
    iterations <- length(nit)
    if (iterations > 1) {
      for (i in 1:(iterations-1)){
        mean_diff <- list$means[i+1] - list$rolling_means[i]
        rolling_variances[i+1] <- rolling_variances[i] + list$variances[i+1]  +
          mean_diff^2*nit[i]*nit[i+1]/sum(nit[1:i])
      }
    }
    tail(rolling_variances, 1)*scaling_factors[scaling_index]^2/(nitot-1)
  }

  l$polarization_variance <- regularize_variance(l, nit, nitot, scaling_factors, scaling_index)


  l$nit <- nitot
  names(l$linext) <- rownames(zeta)[l$linext][l$linext]

  end_time <- Sys.time()

  et <- as.double(difftime(end_time, start_time, units = "secs"))


  #########################
  # CREAZIONE DELL'OUTPUT #
  #########################

  res <- list(
    profiles = profiles,
    measure = measure,
    number_of_profiles = l$n,
    number_of_variables = ncol(profiles$profiles),
    incidence = l$zeta,
    cover = incidence2cover(l$zeta),
    number_of_iterations = l$nit,
    prof_w = weights,
    mean = l$polarization_mean,
    sd = sqrt(l$polarization_variance),
    quantiles = quantiles,
    min = minpol,
    max = maxpol,
    time = et
  )

  class(res) <- "parsec_polarization"

  return(res)
}
