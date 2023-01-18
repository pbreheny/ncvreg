ggplot_ncvreg <- function(x, Y, xlab, col, ...) {
  colnames(Y) <- x
  DF <- as.data.frame.table(Y)
  DF$Var2 <- as.numeric(levels(DF$Var2))[DF$Var2]
  ggplot2::ggplot(DF, ggplot2::aes(Var2, Freq, group=Var1, color=Var1)) +
    ggplot2::geom_line() +
    ggplot2::theme_minimal() + 
    ggplot2::scale_x_reverse() +
    ggplot2::xlab(xlab) + 
    ggplot2::ylab(expression(hat(beta))) + 
    ggplot2::theme(axis.title.y=ggplot2::element_text(angle=0, vjust=0.5), legend.title=ggplot2::element_blank())
}
