#' Title
#'
#' @param mcmc.out
#'
#' @return
#' @export
#'
#' @importFrom tibble as_tibble
#' @importFrom dplyr bind_rows select starts_with mutate
#'
#' @examples
posterior_on_number_of_covariates_selected <- function(mcmc.out, warmup_prop = 0.5) {
  mcmc.out$samples %>%
    add_chain_id_iter() %>%
    filter(it > warmup_prop * max(it)) %>%
    select(starts_with("delta")) %>%
    mutate(nvars = rowSums(.)) %>%
    (function(df) df$nvars) %>%
    tibble(nvars = .) %>%
    ggplot(aes(x = nvars, y = ..density..)) +
    theme_bw() +
    geom_histogram() +
    ggtitle("Posterior density on the number of covariates selected") +
    ylab("") +
    xlab("Number of covariates selected")
}


#' Title
#'
#' @param mcmc.out
#'
#' @return
#' @export
#'
#' @examples
traceplot <- function(mcmc.out, parnames = NULL) {
  to_plot <- mcmc.out$samples %>%
    add_chain_id_iter() %>%
    select(!starts_with("yhat")) %>%
    select(!contains("prior"))

  if (!is.null(parnames)) to_plot <- select(to_plot, starts_with(c(parnames, "it", "chain_id")))

  to_plot %>%
    gather(param, value, -it, -chain_id) %>%
    ggplot(aes(x = it, y = value, colour = factor(chain_id))) +
    theme_bw() +
    facet_wrap(~param, scales = "free_y") +
    geom_line()
}


#' Title
#'
#' @param l
#'
#' @return
#' @export
#'
#' @examples
posterior_predictive_estimate <- function(data, mcmc_out, warmup_prop = 0.5, maxspecies = Inf) {
  # predicted <- data %>%
  #   bind_cols(mcmc_out$samples %>%
  #     add_chain_id_iter() %>%
  #     filter(it > warmup_prop * max(it)) %>%
  #     select(starts_with("yhat")) %>%
  #     colMeans() %>%
  #     tibble(y_meanpred = .))
  #
  # predicted %>%
  #   ggplot(aes(x = times)) +
  #   theme_bw() +
  #   facet_wrap(~i) +
  #   geom_point(aes(y = Y)) +
  #   geom_line(aes(y = y_meanpred), linetype = "dotted")

  predicted <- data %>%
    bind_cols(mcmc_out$samples %>%
      add_chain_id_iter() %>%
      filter(it > warmup_prop * max(it)) %>%
      select(starts_with("yhat")) %>%
      gather(variable, value) %>%
      group_by(variable) %>%
      # mutate(data_index = as.numeric(factor(variable))) %>%
      # group_by(data_index) %>%
      summarise(
        mean = mean(value),
        med = median(value),
        firstq = quantile(value, probs = 0.25),
        thirdq = quantile(value, probs = 0.75),
        infCI = quantile(value, probs = 0.025),
        supCI = quantile(value, probs = 0.975)
      ) %>%
      mutate(data_index = variable %>% sapply(function(str) stringr::str_extract_all(str, pattern = "[0-9]+") %>% .[[1]])) %>%
      mutate(data_index = as.numeric(data_index)) %>%
      arrange(data_index))

  predicted %>%
    filter(i <= maxspecies) %>%
    ggplot(aes(x = times)) +
    theme_bw() +
    facet_wrap(~i) +
    geom_point(aes(y = Y), colour = "red") +
    geom_line(aes(y = mean), linetype = "dotted") +
    geom_ribbon(aes(ymin = infCI, ymax = supCI), linetype = "dotted")
}

#' Title
#'
#' @param mcmc_out
#' @param parnames
#' @param warmup_prop
#'
#' @return
#' @export
#'
#' @examples
PriorPosteriorPlotScalarParams <- function(mcmc_out, parnames, warmup_prop = 0.5) {
  ExtractSamplesPriorPosteriorScalarParams(mcmc_out, parnames, warmup_prop = warmup_prop) %>%
    mutate(type = factor(type, levels = c("prior", "posterior"))) %>%
    ggplot(aes(x = value, y = ..density.., colour = type, fill = type)) +
    theme_bw() +
    facet_wrap(~variable, scales = "free") +
    geom_density(alpha = 0.5) +
    ylab("") +
    xlab("Parameter value") +
    scale_colour_discrete(name = "") +
    scale_fill_discrete(name = "")
}

#' Title
#'
#' @param mcmc_out
#' @param parnames
#' @param warmup_prop
#'
#' @return
#' @export
#'
#' @examples
PriorPosteriorPlotLargeVectorParams <- function(mcmc_out, parnames, warmup_prop = 0.5) {
  to_plot <- ExtractSamplesPriorPosteriorLargeVectorParams(mcmc_out, parnames, warmup_prop = warmup_prop)

  ggplot(data = to_plot, aes(x = value, y = ..density.., colour = type, fill = type)) +
    theme_bw() +
    facet_wrap(~param, scales = "free") +
    geom_density(data = to_plot %>% filter(type == "prior"), alpha = 0.5) +
    geom_histogram(data = to_plot %>% filter(type == "estimate")) +
    ylab("") +
    xlab("Parameter value") +
    scale_colour_discrete(name = "", labels = c("Estimates", "Prior Distribution")) +
    scale_fill_discrete(name = "", labels = c("Estimates", "Prior Distribution"))
}


#' Title
#'
#' @param mcmc_out
#' @param parnames
#' @param warmup_prop
#'
#' @return
#' @export
#'
#' @examples
PriorPosteriorPlotSmallVectorParams <- function(mcmc_out, parnames, warmup_prop = 0.5) {
  to_plot <- ExtractSamplesPriorPosteriorSmallVectorParams(mcmc_out, parnames, warmup_prop = warmup_prop)

  to_plot %>%
    mutate(type = factor(type, levels = c("prior", "posterior"))) %>%
    ggplot(aes(x = value, y = ..density.., colour = type, fill = type)) +
    theme_bw() +
    facet_wrap(~variable, scales = "free") +
    geom_density(alpha = 0.5) +
    ylab("") +
    xlab("Parameter value") +
    scale_colour_discrete(name = "") +
    scale_fill_discrete(name = "")
}

#' Title
#'
#' @param mcmc_out
#' @param parnames
#' @param final_size
#'
#' @return
#' @export
#'
#' @examples
CorrelationPlot = function(mcmc_out, parnames, final_size = 2000, warmup_prop = 0.){
  nit = mcmc_out$samples %>% first %>% nrow()

  startit = max(1, round(nit*warmup_prop))

  iters_to_keep = round(seq(startit, nit, length.out = final_size))

  mcmc_out$samples %>%
    lapply(function(x) x[iters_to_keep,]) %>%
    mapply(function(x, c) as_tibble(x) %>% mutate(chain = c), ., seq_along(.), SIMPLIFY = F) %>%
    bind_rows() %>%
    select(parnames, chain) %>%
    GGally::ggpairs(mapping = aes(color = factor(chain)), columns = parnames) +
    theme_bw()

}
