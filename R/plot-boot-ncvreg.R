#' Title
#'
#' @param eb_boot 
#' @param n 
#' @param quiet 
#'
#' @return
#' @export
#'
#' @examples
plot.boot.ncvreg <- function(eb_boot, n = 30, quiet = TRUE, ci_method = "quantile", original_data = NULL) {
  
  plot_res <- ci.boot.ncvreg(eb_boot, quiet = quiet, ci_method = ci_method, original_data = original_data) %>%
    dplyr::arrange(desc(abs(estimate))) %>%
    head(n)
  
  plot_res$variable <- factor(plot_res$variable, levels = rev(plot_res$variable))
  
  plot_res %>%
    ggplot() +
    geom_errorbar(aes(xmin = lower, xmax = upper, y = variable)) +
    geom_point(aes(x = estimate, y = variable)) +
    theme_bw() +
    labs(y = "Variable", x = "Estimate")
  
}
