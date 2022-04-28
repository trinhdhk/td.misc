#' Function to create a side by side mcmc intervals
#' @import ggplot2
#' @description This function is similar to \code{bayesplot::mcmc_intervals} but
#' provides support for pooling several intervals from several fits.
#' @param fits a list of stan fits, must be wrapped in \code{list()}
#' @param ... other params passed to \code{bayesplot::mcmc_intervals_data}
#' @param multi_point_est default to FALSE, whether to set show multiple point est, which are mean and median
#' @export
mcmc_intervals_multi <-
  rlang::new_function(
    args = append(
      list(fits = list(), multi_point_est = FALSE),
      formals(bayesplot::mcmc_intervals)[-1]
      ),
    body = quote({
      # require(ggplot2)
      w <- .15 * length(fits)

      if (!length(names(fits))){
        names(fits) <- seq_along(fits)
      }

      plot_call <- sys.call()
      plot_call[[1]] <- quote(bayesplot::mcmc_intervals_data)
      if (!'fits' %in% names(plot_call)) plot_call[[2]] <- NULL
      else plot_call[['fits']] <- NULL
      plot_call[names(plot_call)=='multi_point_est'] <- NULL
      if (!'point_est' %in% names(plot_call)) multi_point_est <- FALSE
      point_est <- if (multi_point_est) unique(point_est) else NULL
      est  <- sapply(seq_along(fits),
                     \(i) {
                       this_call <- plot_call
                       this_call[['x']] <- fits[[i]]
                       if (multi_point_est) {
                         this_call[['point_est']] <- point_est[[1]]
                         dt <- eval(this_call) |>
                           dplyr::mutate(Model = names(fits)[[i]])
                         for (j in seq_along(point_est)[-1]){
                           this_call[['point_est']] <- point_est[[j]]
                           dt[[paste0('m',j)]] <- eval(this_call)$m
                         }
                       } else{
                         dt <- eval(this_call) |>
                           dplyr::mutate(Model = names(fits)[[i]])
                       }
                       dt
                     }, simplify=FALSE, USE.NAMES = TRUE)
      est_dat <- dplyr::bind_rows(est) |> dplyr::mutate(Model = factor(Model, levels = names(fits)))
      point_shape = c(21, 3, 4)
      plt <- ggplot(est_dat, aes(y=parameter, color=Model, fill=Model)) +
        geom_linerange(aes(xmin=ll, xmax=hh), position=position_dodge(width=w)) +
        geom_linerange(aes(xmin=l, xmax=h), size=1, position=position_dodge(width=w)) +
        geom_point(aes(x=m), size=point_size, shape=point_shape[1], position=position_dodge(width=w)) +
        xlab('Estimate') + ylab('Parameter')

      if (multi_point_est){
        for (j in seq_along(point_est)[-1]){
          plt <- plt +
            geom_point(aes(x=.data[[paste0('m',j)]]), size=point_size, shape=point_shape[j], color=grey(.1), position=position_dodge(width=w))
        }
      }

      plt
    })
  )

#' Change y tick labels of bayesplot
#' @description Misc. function to change y tick marks of bayesplot
#' @param mcmc_plot an plot returned by bayepslot mcmc_
#' @param ... tick labels to be changed to
#' @param labs same as ... but accept a vector of character
#' @param top_down default to TRUE which respects bayesplot's ordering, otherwise ggplot's ordering.
#' @export
change_ylabs <- function(mcmc_plot, ..., labs = character(), top_down = TRUE){
  labs <- c(..., labs)
  params <- unique(mcmc_plot$data$parameter)
  stopifnot(length(labs) == length(params))
  mcmc_plot + ggplot2::scale_y_discrete(limits = if(top_down) rev(params) else params,
                                        labels = if(top_down) rev(labs) else labs)
}

