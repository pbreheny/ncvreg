#' Plot Bootstrap Results for NCV Regression
#'
#' This function plots the bootstrap confidence intervals and estimates for
#' an NCV regression model. It provides options for ordering the results by
#' different criteria and allows for displaying a subset of variables.
#'
#' @param boot An object containing the bootstrap results from an NCV regression.
#' @param n An integer specifying the number of variables to plot. Default is 30.
#' @param alpha A numeric value for the confidence level, with default 0.2.
#' @param order A character string specifying the order of the variables. Options
#'        are "original", "descending", "ascending", or a vector of variable names.
#' @param absolute_order A logical indicating whether to order by the absolute
#'        value of the estimates. Default is FALSE.
#' @param quiet A logical indicating whether to suppress messages during execution.
#'        Default is TRUE.
#' @return A ggplot object displaying the bootstrap confidence intervals and estimates.
#' @export
#'
#' @examples
#' \dontrun{
#'   boot_results <- boot_ncvreg(...) # Assuming boot_ncvreg() produces the required object
#'   plot.boot_ncvreg(boot_results, n = 20, order = "descending", alpha = 0.05)
#' }
plot.boot_ncvreg <- function(boot, n = 30, alpha = 0.2, order = "original", absolute_order = FALSE, quiet = TRUE) {
  
  # Ensure the order variable is one of the accepted options or a vector of variable names
  if (!is.character(order) && !is.vector(order)) {
    stop("The 'order' parameter must be 'original', 'descending', 'ascending', or a vector of variable names.")
  }
  
  plot_res <- confint(boot, quiet = quiet, alpha = alpha)
  
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
