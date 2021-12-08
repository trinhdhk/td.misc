#' Calibration curve for binary outcome
#' @description Helper Function to create calibration curve for binary outcomes
#' @param pred predictor values/probabilities
#' @param obs observed values
#' @param pred_rep for K-Fold CV, plot the predicted values for all repetitions
#' @param method either loess or splines
#' @param se whether to draw se ribbon
#' @param hist whether to draw histograms for each class
#' @param hist.normalize whether to normalize the height of histograms to make them reflects the true unbalance of the two class.
#' @param ... control parameters passed to the smoothing method, commonly span and knots
#' @param theme theme function. default to theme_bw
#' @return a ggplot
#' @export
binary_calibration <- function(pred, pred_rep = NULL, obs, title = NULL, method=c("loess","splines"), se = TRUE, hist = TRUE, hist.normalize = TRUE, yscale=.1,..., theme = ggplot2::theme_bw){
  require(ggplot2)
  maincolor <- "#CD113B"
  subcolor1 <- "#111111"
  subcolor2 <- "#999999"

  method <- match.arg(method)

  # browser()
  p <- ggplot(mapping=aes(x = pred))
  add_pred <- function(pr, line = FALSE){
    pr <- pr
    if (method == "loess"){
      geom_smooth(aes(x = pr, y = as.numeric(obs)), method = 'loess', color = subcolor1, fill = subcolor2, alpha=if (line) 1 else .5/length(pred_rep), size=line, se = !line,...)
    } else {
      stat_smooth(aes(x = pr, y = as.numeric(obs)), method="glm", formula=as.logical(y)~splines::ns(qlogis(x), knots), method.args=list(family='binomial'), color = subcolor1, fill = subcolor2, size=line, alpha=if (line) 1 else .6/length(pred_rep), se = !line, ...)
    }
  }

  for (pr in pred_rep) p <- p + add_pred(pr)
  for (pr in pred_rep) p <- p + add_pred(pr, line = .3)


  if (method == "loess"){
    p <- p + geom_smooth(aes(y = as.numeric(obs)), method = 'loess', color=maincolor, fill = subcolor2, se = se, size=.9, alpha=.3, ...)
  }
  else {
    p <- p + stat_smooth(aes(y = as.numeric(obs)), fullrange = TRUE, method="glm", formula=as.logical(y)~splines::ns(qlogis(x), knots), method.args=list(family='binomial'), color = maincolor, fill = subcolor2, size=.9, alpha = .3, se = se,...)
  }
  p <- p +
    geom_line(aes(y = pred), color = "#000000") +
    ylab("Observed") + xlab("Predicted")+
    scale_y_continuous(breaks = seq(0, 1, length.out = 5), limits = c(0,1), oob = scales::squish)+
    scale_x_continuous(breaks = seq(0, 1, length.out = 5))+
    ggtitle(title)

  if (hist){
    p0 <-  thinhist_subplot.binary(x = pred, y = as.numeric(obs), normalize = hist.normalize, yscale = yscale, plot = FALSE)
    p <- p + p0[[1]] + p0[[2]]
  }

  # Find average calibration
  obs_p <- sum(obs)/length(obs)
  pred_p <- mean(pred)

  # Find minimum calibration
  # browser()
  pred.logit <- qlogis(pred)
  calib.fit <- glm(as.logical(obs)~pred.logit, family='binomial')

  p +
    ggplot2::annotate(
      "text", x = .99, y = .05, hjust='right', vjust='bottom',
      label = paste0("Observed prevalence: ", format(obs_p, digits = 2), '\n',
                     'Predicted prevalence: ', format(pred_p, digits=2), '\n',
                     "Calibration intercept: ", format(coef(calib.fit)[['(Intercept)']], digits = 2),'\n',
                     'Calibration slope: ', format(coef(calib.fit)[['pred.logit']], digits=2)),
      parse = F, size = 2.8, colour = classifierplots:::fontgrey_str
    ) +
    theme()
}

# Thin histogram
thinhist_subplot <- function(x, normalize.factor = NULL, digits = 2, yscale = 1, zero_y = 0, plot = TRUE){
  round_factor = 10^digits
  sub = floor(min(x*round_factor))/round_factor
  sup = ceiling(max(x*round_factor))/round_factor
  breaks = seq(sub, sup, 10^-digits)
  x.breaks = hist(x, breaks, plot = FALSE)
  if (!length(normalize.factor)) normalize.factor <- max(x.breaks$count)
  counts = (x.breaks$count / normalize.factor) * yscale + zero_y
  require(ggplot2)
  p <- geom_segment(aes(y = zero_y, x = x.breaks$mids, yend = counts, xend = x.breaks$mids))
  if (!plot) return(p)
  ggplot() + p  +
    ylab("") + xlab("")
}

thinhist_subplot.binary <- function(x, y, normalize = TRUE, digits = 2, yscale = 1, plot = TRUE){
  round_factor = 10^digits
  sub = floor(min(x*round_factor))/round_factor
  sup = ceiling(max(x*round_factor))/round_factor
  breaks = seq(sub, sup, 10^-digits)
  x.breaks = hist(x, breaks, plot = FALSE)
  normalize.factor <- if (normalize) quantile(x.breaks$count, .95, na.rm=TRUE) else NULL
  p <-
    list(
      thinhist_subplot(x[y==0], normalize.factor=normalize.factor, digits = digits, yscale = yscale, zero_y=0, plot=FALSE),
      thinhist_subplot(x[y==1], normalize.factor=normalize.factor, digits = digits, yscale = -yscale, zero_y=1, plot=FALSE)
    )
  if (!plot) return(p)
  ggplot() + p[[1]] + p[[2]] +
    ylab("") + xlab("")
}
