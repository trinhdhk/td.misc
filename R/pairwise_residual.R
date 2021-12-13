#' Parwise correlation plot for LCA
#' @description Function to calculate pairwise residual correlation
#' @param obs a matrix of observed binary outcomes
#' @param pred a matrix of predicted probability for corresponding outcome
#' @param B number of bootstrap
#' @seealso (ttps://www.jstor.org/stable/2533043?origin=crossref)
#' @export
pairwise_corr_plot <-
  function(obs, pred, B = 0){
    obs <- as.matrix(obs); pred <- as.matrix(pred)
    obs <- obs[complete.cases(obs),]
    pred <- pred[complete.cases(obs),]
    if (length(colnames(obs)) && length(colnames(pred))) pred <- pred[names(obs)]

    com.matrix <- function(M)
      do.call(expand.grid,split(M,rep(1:nrow(M),ncol(M))))

    pair_corr <- function(p, corr_list){
      pair <- pair.matrix[p,]
      i <- pair[[1]]; j <- pair[[2]]
      mu_i <- corr_list$mu_i[[i]]
      mu_j <- corr_list$mu_i[[j]]
      mu_ij <- corr_list$mu_ij[[p]]
      (mu_ij  - mu_i*mu_j)/sqrt(mu_i*(1-mu_i)*mu_j*(1-mu_j))
    }
    n.test <- ncol(obs)
    pair.matrix <- com.matrix(rbind(seq_len(n.test), seq_len(n.test)))
    pair.matrix <- subset(pair.matrix, pair.matrix[,1]!=pair.matrix[,2])
    # browser()
    pred.sample <- do.call(
      rbind,
      apply(pred, 1,
            function(x)
              sapply(x, function(y) rbinom(100,1,y)),
            simplify = FALSE))

    corr_obs <- corr_pred <- list()
    corr_obs$mu_i <- apply(obs, 2, mean)
    corr_obs$mu_ij <- apply(pair.matrix, 1,
                             function(pair)
                               mean(obs[,pair[1]]*obs[,pair[2]]))

    corr_obs$corr_ij <-
      sapply(seq_len(nrow(pair.matrix)), pair_corr, corr_list = corr_obs)

    corr_pred$mu_i <- apply(pred.sample, 2, mean)
    corr_pred$mu_ij <- apply(pair.matrix, 1,
                              function(pair)
                                mean(pred.sample[,pair[1]]*pred.sample[,pair[2]]))

    corr_pred$corr_ij <-
      sapply(seq_len(nrow(pair.matrix)), pair_corr, corr_list = corr_pred)

    corr_resid <- corr_pred$corr_ij - corr_obs$corr_ij
    names(corr_resid) <- names(corr_pred$corr_ij) <- names(corr_obs$corr_ij) <-
                                 apply(pair.matrix,1,
                                       function(x){
                                         if (length(colnames(obs))){
                                           name.obs <- colnames(obs)
                                           return(paste(name.obs[x[[1]]], name.obs[x[[2]]], sep='-'))
                                         }
                                         paste(x[[1]],x[[2]],sep='-')
                                       }
                                 )


    return(
      list(
        corr_pred = corr_pred$corr_ij,
        corr_obs = corr_obs$corr_ij,
        corr_resid = corr_resid
      )
    )
  }
