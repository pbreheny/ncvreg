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
plot.boot_ncvreg <- function(eb_boot, method = "hybrid", n = 30, alpha = 0.2, order = "original", absolute_order = FALSE, quiet = TRUE) {
  
  # Ensure the order variable is one of the accepted options or a vector of variable names
  if (!is.character(order) && !is.vector(order)) {
    stop("The 'order' parameter must be 'original', 'descending', 'ascending', or a vector of variable names.")
  }
  
  plot_res <- ci.boot_ncvreg(eb_boot, quiet = quiet, alpha = alpha, methods = method)
  
  # Handle the ordering of the results
  if (length(order) == 1) {
    if (order == "descending") {
      if (absolute_order) {
        plot_res <- plot_res %>%
          dplyr::arrange(desc(abs(estimate)))
      } else {
        plot_res <- plot_res %>%
          dplyr::arrange(desc(estimate))
      }
    } else if (order == "ascending") {
      if (absolute_order) {
        plot_res <- plot_res %>%
          dplyr::arrange(abs(estimate))
      } else {
        plot_res <- plot_res %>%
          dplyr::arrange(estimate)
      }
    } 
  } else {
    # Ensure the provided vector is a proper subset of the variables in plot_res
    if (!all(order %in% plot_res$variable)) {
      stop("All elements in the 'order' vector must be present in the variable names of the plot_res data frame.")
    }
    plot_res$variable <- factor(plot_res$variable, levels = order)
    plot_res <- plot_res %>%
      dplyr::arrange(factor(variable, levels = order))
  }
  
  plot_res <- plot_res %>%
    head(n) 

  if (order == "original") {
    plot_res$variable <- factor(plot_res$variable, levels = rev(plot_res$variable))
  }

  # Plot the results
  plot_res %>%
    ggplot() +
    geom_errorbar(aes(xmin = lower, xmax = upper, y = variable)) +
    geom_point(aes(x = estimate, y = variable)) +
    theme_bw() +
    labs(y = "Variable", x = "Estimate")
}
