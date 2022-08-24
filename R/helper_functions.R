#' @importFrom magrittr %>%
LoglogisticIncreasing <- function(t, phi1, phi2, phi3) {
  phi1 / (1 + exp((phi2 - t) / phi3))
}
