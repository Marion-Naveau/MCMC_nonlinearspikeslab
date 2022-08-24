#' Title
#'
#' @param p
#' @param N
#' @param ni
#' @param times
#' @param alpha
#' @param beta
#' @param Gamma
#' @param sigma
#'
#' @return
#' @export
#'
#' @examples
SimulateDataOneVariableParam <- function(p = 5, N = 10, ni = 30, times = seq(0, 200, length.out = N), alpha = 15, beta = rbinom(n = p, size = 1, prob = 0.2), Gamma = 0.1, sigma = 0.1, param = "final size") {
  V <- cbind(rep(1, ni), rnorm(n = ni * p) %>% matrix(ncol = p))
  phi <- V %*% as.matrix(c(alpha, beta)) + rnorm(ni, sd = Gamma)

  if (param == "final size") {
    res <- SimulateDataOneVariableParamFinalSize(V = V, phi = phi, ni = ni, times = times, beta = beta, sigma = sigma, alpha_ = alpha)
  } else if (param == "delay") {
    res <- SimulateDataOneVariableParamDelay(V = V, phi = phi, ni = ni, times = times, beta = beta, sigma = sigma, alpha = alpha)
  } else {
    stop("simulation function not implemented yet")
  }

  res %>%
    mutate(alpha = alpha) %>%
    return()
}

#' Title
#'
#' @param N
#' @param ni
#' @param times
#' @param alpha
#' @param beta
#' @param Gamma
#' @param sigma
#' @param param
#' @param Lmax
#' @param Tmax
#' @param psi1
#' @param psi2
#' @param psi3
#'
#' @return
#' @export
#'
#' @examples
SimulateDataOneVariableParam_reparam <- function(J = 10, n = 30, times = seq(0, 200, length.out = J), beta = c(0.3, 0.3, 0, 0, 0, 0), Gamma = 0.1, sigma = 0.1, param = "delay", Lmax = 15, psi1 = 1, psi2 = 0.5, psi3 = 0.5) {
  p <- length(beta)
  Tmax <- max(times)
  V <- cbind(rep(1, n), rnorm(n = n * p) %>% matrix(ncol = p))
  phi <- V %*% as.matrix(c(1, beta)) + rnorm(n, sd = Gamma)

  if (param == "delay") {
    res <- SimulateDataOneVariableParamDelay_reparam(V = V, phi = phi, ni = n, times = times, beta = beta, sigma = sigma, Lmax = Lmax, Tmax = Tmax, psi1 = psi1, psi2 = psi2, psi3 = psi3)
  } else {
    stop("simulation function not implemented yet")
  }

  res %>%
    return()
}



SimulateDataOneVariableParamFinalSize <- function(V, phi, ni, times, beta, sigma, phi2 = 0, alpha_) {
  phi3 <- median(times) / 4
  lapply(1:ni, function(i) {
    tibble::tibble(times = times, i = i, phi_i = phi[i], gphi = LoglogisticIncreasing(t = times, phi1 = phi[i], phi2 = phi2, phi3 = phi3))
  }) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(Y = gphi + rnorm(n = length(gphi), sd = sigma)) %>%
    dplyr::left_join(tibble::rowid_to_column(V %>% tibble::as_tibble(), var = "i")) %>%
    dplyr::mutate(phi2 = 0, phi3 = phi3) %>%
    dplyr::bind_cols(., tibble::tibble(alpha = rep(alpha_, nrow(.))), lapply(1:nrow(.), function(i) {
      beta %>%
        t() %>%
        tibble::as_tibble() %>%
        (function(df) {
          df %>%
            names() %>%
            gsub("V", "beta", .) %>%
            setNames(df, .)
        })
    }) %>% dplyr::bind_rows()) %>%
    dplyr::select(-V1)
}

SimulateDataOneVariableParamDelay <- function(V, phi, ni, times, beta, sigma, phi1 = 15, alpha) {
  phi3 <- median(times) / 4
  lapply(1:ni, function(i) {
    tibble::tibble(times = times, i = i, phi_i = phi[i], gphi = LoglogisticIncreasing(t = times, phi1 = phi1, phi2 = phi[i], phi3 = phi3))
  }) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(Y = gphi + rnorm(n = length(gphi), sd = sigma)) %>%
    dplyr::left_join(tibble::rowid_to_column(V %>% tibble::as_tibble(), var = "i")) %>%
    dplyr::mutate(phi2 = 0, phi3 = phi3) %>%
    dplyr::bind_cols(., tibble::tibble(alpha = rep(alpha, nrow(.))), lapply(1:nrow(.), function(i) {
      beta %>%
        t() %>%
        tibble::as_tibble() %>%
        (function(df) {
          df %>%
            names() %>%
            gsub("V", "beta", .) %>%
            setNames(df, .)
        })
    }) %>% dplyr::bind_rows()) %>%
    dplyr::select(-V1)
}

SimulateDataOneVariableParamDelay_reparam <- function(V, phi, ni, times, beta, sigma, Lmax = 15, Tmax = 10, psi1 = 1, psi2 = 0.5, psi3 = 0.5) {
  lapply(1:ni, function(i) {
    tibble::tibble(times = times, i = i, phi_i = phi[i], gphi = LoglogisticIncreasing(t = times, phi1 = Lmax * psi1, phi2 = Tmax * psi2 * phi[i], phi3 = Tmax * psi3))
  }) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(Y = gphi + rnorm(n = length(gphi), sd = sigma)) %>%
    dplyr::left_join(tibble::rowid_to_column(V %>% tibble::as_tibble(), var = "i")) %>%
    dplyr::bind_cols(., lapply(1:nrow(.), function(i) {
      beta %>%
        t() %>%
        tibble::as_tibble() %>%
        (function(df) {
          df %>%
            names() %>%
            gsub("V", "beta", .) %>%
            setNames(df, .)
        })
    }) %>% dplyr::bind_rows()) %>%
    dplyr::select(-V1) %>%
    mutate(Lmax = Lmax, Tmax = Tmax, psi1 = psi1, psi2 = psi2, psi3 = psi3)
}


#' Title
#'
#' @param data
#'
#' @return
#' @export
#'
#' @examples
PlotSimulatedData <- function(data, maxspecies = Inf) {
  data %>%
    filter(i <= maxspecies) %>%
    ggplot2::ggplot(ggplot2::aes(x = times, y = Y, colour = factor(i))) +
    ggplot2::theme_bw() +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    scale_colour_discrete(name = "Individual") +
    xlab("Time") +
    theme(legend.position = "none")
}
