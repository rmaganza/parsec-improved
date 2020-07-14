summary.parsec_polarization <-
  function(object, ...) {
    if (!is.null(object$profiles)) {
      res <- object$profiles$profiles
    } else {
      res <- data.frame(tmp = 1:object$number_of_profiles)
    }

    res$weights <- object$prof_w

    if (object$number_of_profiles < 11)
      print(res)
    else {
      print(res[1:5,])
      cat("\n...\n\n")
      print(res[object$number_of_profiles - 4:0,])
    }

    try (
      {
        cat("\nComputed polarization measure", object$measure, "in", object$time, "seconds")
      })

    try(
      {
        cat("\n\nMean: ", object$mean)
        cat("\nSD: ", object$sd)
        cat("\nMin: ", object$min)
        cat("\nMax: ", object$max)
      }
    )
    invisible(res)
  }
