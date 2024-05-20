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
plot.boot.ncvreg <- function(eb_boot, n = 30, alpha = 0.2, original_order = FALSE, absolute_order = FALSE, quiet = TRUE) {
  
  plot_res <- ci.boot.ncvreg(eb_boot, quiet = quiet, alpha = alpha) 
  
  if (!original_order) {
    plot_res <- plot_res %>%
      dplyr::arrange(desc(abs(estimate))) 
  }
  
  plot_res <- plot_res %>%
    head(n)
  
  if (!original_order) {
    if (!absolute_order) {
      plot_res <- plot_res %>%
        dplyr::arrange(desc(estimate)) 
    }
  }
  
  plot_res$variable <- factor(plot_res$variable, levels = rev(plot_res$variable))
  
  plot_res %>%
    ggplot() +
    geom_errorbar(aes(xmin = lower, xmax = upper, y = variable)) +
    geom_point(aes(x = estimate, y = variable)) +
    theme_bw() +
    labs(y = "Variable", x = "Estimate")
  
}
