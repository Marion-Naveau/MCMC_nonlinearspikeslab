add_chain_id_iter <- function(codamcmc) {
  codamcmc %>%
    lapply(as_tibble) %>%
    (function(l) {
      lapply(X = 1:length(l), FUN = function(idx) {
        l[[idx]] %>%
          mutate(chain_id = idx) %>%
          rowid_to_column("it")
      })
    }) %>%
    bind_rows()
}

ExtractSamplesPriorPosteriorScalarParams <- function(mcmc_out, parnames, warmup_prop = 0.5) {
  nsamples <- mcmc_out$samples %>%
    first() %>%
    nrow()

  retained_smpls <- mcmc_out$samples %>%
    lapply(as_tibble) %>%
    lapply(rowid_to_column) %>%
    bind_rows() %>%
    filter(rowid >= warmup_prop * nsamples)
  retained_smpls %>%
    select(all_of(parnames)) %>%
    mutate(type = "posterior") %>%
    gather(variable, value, -type) %>%
    bind_rows(retained_smpls %>%
                select(all_of(paste(parnames, "_prior", sep = ""))) %>%
                setNames(parnames) %>%
                mutate(type = "prior") %>%
                gather(variable, value, -type))
}

ExtractSamplesPriorPosteriorLargeVectorParams <- function(mcmc_out, parnames, warmup_prop = 0.5) {
  nsamples <- mcmc_out$samples %>%
    first() %>%
    nrow()

  retained_smpls <- mcmc_out$samples %>%
    lapply(as_tibble) %>%
    lapply(rowid_to_column) %>%
    bind_rows() %>%
    filter(rowid >= warmup_prop * nsamples)

  parnames %>%
    lapply(function(parname) {
      retained_smpls %>%
        select(starts_with(parname)) %>%
        select(-contains("prior")) %>%
        colMeans() %>%
        tibble(param = parname, value = ., type = "estimate") %>%
        bind_rows(retained_smpls %>%
                    select(all_of(paste(parname, "_prior", sep = ""))) %>%
                    setNames(c("value")) %>%
                    mutate(type = "prior", param = parname))
    }) %>%
    bind_rows()
}


ExtractSamplesPriorPosteriorSmallVectorParams <- function(mcmc_out, parnames, warmup_prop = 0.5) {
  nsamples <- mcmc_out$samples %>%
    first() %>%
    nrow()

  retained_smpls <- mcmc_out$samples %>%
    lapply(as_tibble) %>%
    lapply(rowid_to_column) %>%
    bind_rows() %>%
    filter(rowid >= warmup_prop * nsamples)

  parnames %>%
    lapply(function(parname) {
      retained_smpls %>%
        select(starts_with(parname)) %>%
        select(-contains("prior")) %>%
        mutate(type = "posterior") %>%
        (function(df){
          vector_length = ncol(df)-1
        bind_rows(df, retained_smpls %>%
                    select(all_of(paste(parname, "_prior", sep = ""))) %>%
                    (function(col){
                      lapply(1:vector_length, function(x) col) %>%
                        bind_cols()
                    }) %>%
                    mutate(type = "prior") %>%
                    setNames(names(df))
                  ) %>%
                    gather(variable, value, -type)
        })
    }) %>%
    bind_rows()
}
