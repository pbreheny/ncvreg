\name{Lung}
\alias{Lung}
\title{VA lung cancer data set}
\description{Data from a randomised trial of two treatment regimens for
  lung cancer. This is a standard survival analysis data set from the
  classic textbook by Kalbfleisch and Prentice.
  \itemize{
    \item X: A design matrix with 137 observations (rows) and 9
    predictor variables (columns).  The columns are:
    \itemize{
      \item trt: Treatment indicator (1=control group, 2=treatment
      group)
      \item karno: Karnofsky performance score (0=bad, 100=good)
      \item diagtime: Time from diagnosis to randomization (months)
      \item age: Age (years)
      \item prior: Prior therapy (0=no, 1=yes)
      \item squamous: Indicator for whether the cancer type is squamous
      cell carcinoma (0=no, 1=yes)
      \item small: Indicator for whether the cancer type is small cell
      lung cancer (0=no, 1=yes) 
      \item adeno: Indicator for whether the cancer type is
      adenocarcinoma (0=no, 1=yes) 
      \item large: Indicator for whether the cancer type is large cell
      carcinoma (0=no, 1=yes) 
    }
    \item y: A two column matrix (\code{Surv} object) containing the
    follow-up time (in days) and an indicator variable for whether the
    patient died while on the study or not.
  }
}
\usage{data(Lung)}
\format{A list containing the design matrix \code{X} and response matrix
  \code{y}}
\source{\url{http://CRAN.R-project.org/package=survival}}
\references{
  \itemize{
    \item Kalbfleisch D and Prentice RL (1980), \emph{The Statistical
  Analysis of Failure Time Data}. Wiley, New York.
  }
}
\seealso{\code{\link{ncvsurv}}}
\keyword{datasets}
