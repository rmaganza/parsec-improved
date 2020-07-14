runC.pre_parsec_polarization <-
  function(l) {
    v <- .C("bd_polarization",
            measure = as.integer(l$measure),
            alpha = as.double(l$alpha),
            beta = as.double(l$beta),
            linext = as.integer(l$linext - 1),
            n_scal = as.integer(l$n),
            nit_scal = as.integer(l$nit),
            zeta = as.integer(t(l$zeta)),
            relfreqs = as.double(l$relfreqs),
            polarization_mean = 0.0,
            polarization_variance = 0.0,
            quantiles_q = as.double(l$quantiles_q),
            quantiles_n = as.integer(l$quantiles_n),
            quantiles_f = as.double(l$quantiles_f),
            quantiles_d = as.double(l$quantiles_d),
            polar_min = as.double(l$polar_min),
            polar_max = as.double(l$polar_max),
            computequantiles = as.integer(l$computequantiles),
            firstit = as.integer(l$firstit),
            probvec_bd = as.double(l$probvec)
    )
    # aggiorno solo quello che serve
    l$linext <- v$linext + 1

    l$means <- c(l$means, v$polarization_mean)
    l$polarization_mean <- l$polarization_mean + v$polarization_mean * l$nit
    l$rolling_mean <- c(l$rolling_mean, l$polarization_mean)
    l$variances <- c(l$variances, v$polarization_variance)

    l$quantiles_q <- v$quantiles_q
    l$quantiles_n <- v$quantiles_n
    l$quantiles_f <- v$quantiles_f
    l$quantiles_d <- v$quantiles_d
    l$polar_min <- v$polar_min
    l$polar_max <- v$polar_max
    l$firstit <- 0
    return(l)
  }
