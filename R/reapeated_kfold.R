#' Split data to K for N times
#' @description Create K fold splitted stan inputs N times
#' @param K number of folds
#' @param N_rep number of repetition
#' @param N_obs number of observation in the full dataset
#' @param data a list that is supposed to feed into Stan
#' @param seed seed
#' @param cores number of cores to be used. On Windows, will fall back to single core.
#' @return a nested list with 3 elements:
#'
#' - keptin: list of keptin(1) or holdout(0) in each folds
#' - holdout: negated version of keptin
#' - inputs: processed dataset with a keptin variable injected in
#'
#' @export
repeated_kfold <- function(K = 10, N_rep = 4, N_obs, data, seed = 123, cores = 10){
  if (K == 1) return(
    list(keptin = list(rep(1L, N_obs)), holdout = list(rep(0L, N_obs)),
         inputs = list(modifyList(data, list(keptin = rep(1L, N_obs)))))
  )
  holdout <- vector("list", N_rep * K)
  keptin <- vector("list", N_rep * K)
  for(n in seq_len(N_rep)){
    set.seed(seed + n - 1)
    hh <- loo::kfold_split_random(K = K, N = N_obs)
    foldkept <- matrix(1, nrow = N_obs, ncol = K)
    for(i in 1L:N_obs) foldkept[i, hh[i]] <- 0
    foldhold  <- 1L - foldkept
    keptin[((n-1)*K+1):(n*K)] <- split(foldkept,rep(1:ncol(foldkept),each=nrow(foldkept)))
    holdout[((n-1)*K+1):(n*K)]  <- split(foldhold,rep(1:ncol(foldhold),each=nrow(foldhold)))
  }
  inputs <- pbmcapply::pbmclapply(seq_len(K * N_rep),
                                  mc.cores=cores,
                                  function(k) modifyList(data, list(keptin=keptin[[k]])))
  list(keptin=keptin, holdout=holdout, inputs = inputs)
}

#' Stan fit with K fold CV
#' @description Perfoming Stan fit for K-Fold CV data created by repeated_kfold
#' @param file a stan file
#' @param sampler a precompiled stanmodel or CmdStanModel object
#' @param list_of_datas list of data, if outputted by repeated_kfold: outputs$inputs
#' @param include_paths path to included file. See rstan/cmdstanr help documentation
#' @param sample_dir directory to put sample csv files
#' @param backend either rstan or cmdstanr. if sampler is available, will be chosen accordingly.
#' @param chains number of MCMC chains
#' @param cores number of cores to use
#' @param seed seed
#' @param pars list of parameters to be included
#' @param merge whether to merge the samples to one
#' @param ... other parameters passed to engine call (rstan::sampling or CmdStanModel$sample)
#' @return
#' If merge == TRUE: a stanfit object (if backend = rstan) or sample array (if backend = cmdstanr)
#' If merge == FALSE: a list of stanfit objects (if backend = rstan) or list of CmdStanMCMC objects (if backend = cmdstanr)
#' @export
stan_kfold <- function(file, sampler, list_of_datas, include_paths=NULL, sample_dir = NULL, backend = c("cmdstanr", "rstan"), chains, cores = getOption('mc.cores', default = 1), seed, pars = NA, merge = TRUE,...){
  backend <- match.arg(backend)
  badRhat <- 1.1 # don't know why we need this?
  n_fold <- length(list_of_datas)
  if (missing(sampler)){
    model <- if (backend == "cmdstanr")
      cmdstanr::cmdstan_model(stan_file=file, include_paths = include_paths)
    else rstan::stan_model(file=file)
  } else {
    model <- sampler
    backend2 <- switch(class(sampler)[[1]], stanmodel = "rstan", CmdStanModel = "cmdstanr")
    if (backend != backend2){
      cat("backend set to ", backend2)
      backend <- backend2
    }
  }
  if (length(list_of_datas) == 1L)
    return(normal_stan(model, list_of_datas[[1]], include_path, sample_dir, backend,chains, cores, seed, pars,...))
  # First parallelize all chains:
  future::plan(future::multicore, workers=cores, gc=TRUE)
  wd <- getwd()
  progressr::handlers(progressr::handler_progress)
  progressr::with_progress({
    p <- progressr::progressor(steps = n_fold*chains)
    sflist <- furrr::future_map(seq_len(n_fold*chains), function(i){
      setwd(wd)
      k <- ceiling(i / chains)
      if (backend == "rstan"){
        if (length(sample_dir)){
          sample_file <- file.path(sample_dir, paste0("data_", k, "_chain_", chains-(k*chains-i), ".csv"))
          sf <- rstan::sampling(model, data = list_of_datas[[k]], chains=1, seed = seed,
                                chain_id = i, sample_file = sample_file, pars = pars,...)
        } else

          sf <- rstan::sampling(model, data = list_of_datas[[k]], chains=1, seed = seed,
                                chain_id = i, pars = pars, ...)
      }

      else {
        m <- model$clone()
        sf <- m$sample(data = list_of_datas[[k]],
                       chains = 1, seed = seed, chain_ids = i,
                       output_dir = sample_dir,...)
      }
      p()
      sf
    }, .options = furrr::furrr_options(seed=TRUE))
  }, enable = TRUE)

  # Then merge the K * chains to create K stanfits:
  if (!merge) return(sflist)
  stanfit <- vector('list', n_fold)
  if (backend=="rstan"){
    for (k in seq_len(n_fold)){
      inchains <- (chains*k - (chains - 1)):(chains*k)
      #  Merge `chains` of each fold
      stanfit[[k]] <- rstan::sflist2stanfit(sflist[inchains])
    }
    return(stanfit)
  }

  # if backend='cmdstanr'
  draws <- if (!all(is.na(pars))){
    pbmcapply::pbmclapply(sflist,
                          \(x) posterior::as_draws(x$draws(variables=pars)), mc.cores = cores)
  } else {
    pbmcapply::pbmclapply(sflist,
                          \(x) posterior::as_draws(x$draws()), mc.cores = cores)
  }
  for (k in seq_len(n_fold)) {
    inchains <- (chains*k - (chains - 1)):(chains*k)
    stanfit[[k]] <- do.call(posterior::bind_draws,
                            append(draws[inchains], list(along = 'chain')))
    # dimnames(stanfit[[k]])[[2]] <- seq_len(dim(stanfit[[k]])[2])
    # names(dimnames(stanfit[[k]])) <- c('iteration', 'chain', 'variable')
  }

  stanfit
}

# Normal stan
normal_stan <- function(model, data, include_paths=NULL, sample_dir = NULL, backend, chains, cores, seed, pars,...){
  sample_file <- file.path(sample_dir, 'data_chain')
  if (backend == 'rstan') {
    sf <- rstan::sampling(model, data=data, sample_file = sample_file, chains=chains, cores=cores, seed=seed, pars=pars,... )
  } else {
    m <- model$clone()
    sf <- m$sample(data = data,
                   chains = 1, seed = seed,
                   output_dir = sample_dir,...)
  }

  sf
}

#' Extract loglik from K-Fold CV
#' Copied from AhVehtari
#' @description This function extract loglik from repeated K Fold CV, modified from code of Stan dev
#' @param list_of_stanfits list of stanfits, output from stan_kfold
#' @param list_of_holdout a list of holdoult, output from repeated_kfold$holdout
#extract log-likelihoods of held-out data
extract_log_lik_K <- function(list_of_stanfits, list_of_holdout, loglik_par = 'log_lik',...){
  K <- length(list_of_stanfits)
  holdout_matrix <- simplify2array(list_of_holdout)
  Nrep <- sum(holdout_matrix[1,])
  K <- ncol(holdout_matrix)/Nrep
  extract_log_lik <-
    if (inherits(list_of_stanfits[[1]], 'stanfit')) loo::extract_log_lik
  else \(x) {
    ll <- grep(paste0('^',loglik_par,'\\[?'), dimnames(x)[[3]])
    x[,,ll]
  }
  list_of_log_liks <- lapply(seq_len(Nrep*K), function(k){
    extract_log_lik(list_of_stanfits[[k]],...)
  })
  # `log_lik_heldout` will include the loglike of all the held out data of all the folds.
  # We define `log_lik_heldout` as a (samples x N_obs) matrix
  # (similar to each log_lik matrix)
  N_obs <- dim(list_of_log_liks[[1]])[2]
  # log_lik_heldout <- list_of_log_liks[[1]] * NA
  log_lik_heldout <-
    do.call(rbind,
            lapply(seq_len(Nrep),
                   function(n){
                     log_lik_heldout_n <- list_of_log_liks[[1]] * NA
                     for (k in ((n-1)*K+1):(n*K)){
                       log_lik <- list_of_log_liks[[k]]
                       samples <- dim(log_lik)[1]
                       N_obs <- dim(log_lik)[2]
                       # This is a matrix with the same size as log_lik_heldout
                       # with 1 if the data was held out in the fold k
                       heldout <- matrix(rep(list_of_holdout[[k]], each = samples), nrow = samples)
                       # Sanity check that the previous log_lik is not being overwritten:
                       if(any(!is.na(log_lik_heldout_n[heldout==1]))){
                         warning("Heldout log_lik has been overwritten!!!!")
                       }
                       # We save here the log_lik of the fold k in the matrix:
                       log_lik_heldout_n[heldout==1] <- log_lik[heldout==1]
                     }
                     log_lik_heldout_n
                   })
    )

  attr(log_lik_heldout, "K") <- K
  class(log_lik_heldout) <- 'kfoldll'
  return(log_lik_heldout)
}

#compute ELPD
#' Computation of ELPD
#' @description method for computation of elpd in loo
#' @param log_lik_heldout come from extract_log_lik_K
#' @return an object of class c("kfold", "loo"). See \link[loo]{kfold}
#' @importFrom loo kfold
#' @method kfold kfoldll
#' @export
kfold.kfoldll <- function(log_lik_heldout)  {
  # library(matrixStats)
  logColMeansExp <- function(x) {
    # should be more stable than log(colMeans(exp(x)))
    S <- nrow(x)
    matrixStats::colLogSumExps(x, na.rm=TRUE) - log(S)
  }
  # See equation (20) of @VehtariEtAl2016
  pointwise <-  matrix(logColMeansExp(log_lik_heldout), ncol= 1)
  colnames(pointwise) <- "elpd"
  # See equation (21) of @VehtariEtAl2016
  elpd_kfold <- sum(pointwise)
  se_elpd_kfold <-  sqrt(ncol(log_lik_heldout) * var(pointwise))
  out <- list(
    pointwise = structure(pointwise, dimnames = list(NULL, "elpd_kfold")),
    estimates = structure(cbind(elpd_kfold, se_elpd_kfold), dimnames = list('elpd_kfold', c("Estimate", "SE")))
  )
  attr(out, "K") <- attr(log_lik_heldout, "K")
  class(out) <- c("kfold", "loo")
  return(out)
}

#' Extract hold-out individuals' parameters from K Fold
#' @description Function to extract hold-out individual parameters aka parameters that have the same dim as number of observation
#' @param list_of_stanfits a list of fits coming out from stan_kfold
#' @param list_of_holdouts a list of holdout coming from repeate_kfold$holdout
#' @param pars parameters
#' @param ... additional params passed to the extract method
#' @param holdout Default to TRUE. whether to extract the holdout part or keptin part.
#' @return a list of array of extracted values
#' @export
extract_K_fold <- function(list_of_stanfits, list_of_holdouts, pars = NULL, ...,  holdout=TRUE){
  K = length(list_of_holdouts)
  Nrep = sum(simplify2array(list_of_holdouts)[1,])
  holdout <- as.numeric(holdout)
  stopifnot(length(list_of_stanfits) == K)
  D = if (holdout) 1 else (K/Nrep - 1)
  FUN =
    if (inherits(list_of_stanfits[[1]], 'stanfit')) rstan::extract
  else \(x, pars, ...) {
    ll <- grep(paste0('^',pars,'\\[?'), dimnames(x)[[3]])
    x[,,ll]
  }
  par_extract_list <- lapply(list_of_stanfits, FUN = FUN, pars=pars, ...)
  extract_pars <- names(par_extract_list[[1]])
  extract_holdout <- lapply(extract_pars, function(p) {
    extract_holdout_par <- array(NA, dim = c(dim(par_extract_list[[1]][[p]])[1] * Nrep * D, dim(par_extract_list[[1]][[p]])[-1]))
    for (n in seq_len(Nrep)){
      for (k in seq_len(K2 <- K/Nrep)){
        holdout_k <- list_of_holdouts[[(n-1)*K2+k]]
        for (d in seq_len(D)){
          k_par <- par_extract_list[[(n-1)*K2+k]][[p]]
          par_dims <- dim(k_par)
          if (length(par_dims) >= 2 && par_dims[min(length(par_dims),2)] == length(holdout_k)) {
            commas <- paste(rep(',', length(par_dims)-2), collapse = " ")
            call_string <- glue::glue('extract_holdout_par[(dim(par_extract_list[[1]][[p]])[1]*(n-1)*d+1):(dim(par_extract_list[[1]][[p]])[1]*n*d),
                                    holdout_k=={holdout} {commas}] <- k_par[,holdout_k=={holdout} {commas}]')
            eval(parse(text=call_string))
          } else {
            stop('Par(s) is/are not individual parameter(s). Please use the rstan::extract or $draws() instead!')
          }
        }
      }
    }

    extract_holdout_par
  })
  setNames(extract_holdout, extract_pars)
}

#' Extract variables from a stanfit of KFold stanfit
#' @description Function to extract variables from draws array
#' @param x a draw array
#' @param variables a vector of variables
#' @param merge_chains merge all chains to one
#' @return an array of 3 dimensions: iteration, chain, and variables if merge = FALSE or a matrix of 2 dimensions: draw and variable if merge = TRUE
#' @return a draw matrix
#' @export
extract_variables <- function(x, variables, ..., merge_chains = FALSE){
  match.vars <-
    unlist(lapply(variables,
                  function(v)
                    grep(paste0('^', v, '\\[?'), dimnames(x)[[3]], value=TRUE)))
  if (merge_chains) return({
    m <- sapply(match.vars,
                function(v) posterior::extract_variable(x, v, ...))
    names(dimnames(m)) <- c('draw', 'variable')
    m
  })

  x[,,match.vars]
}

