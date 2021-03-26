\name{Prostate}
\alias{Prostate}
\alias{prostate}
\title{Factors associated with prostate specific antigen}
\description{Data from a study by by Stamey et al. (1989) to examine the
  association between prostate specific antigen (PSA) and several
  clinical measures that are potentially associated with PSA in men who
  were about to receive a radical prostatectomy.  The variables are as
  follows: 

  \itemize{
    
    \item X: A design matrix with 97 instances (rows) and 8
    predictor variables (columns).  The columns are:
    \itemize{
      \item lcavol: Log cancer volume
      \item lweight: Log prostate weight
      \item age: The man's age
      \item lbph: Log of the amount of benign hyperplasia
      \item svi: Seminal vesicle invasion; 1=Yes, 0=No
      \item lcp: Log of capsular penetration
      \item gleason: Gleason score
      \item pgg45: Percent of Gleason scores 4 or 5
    }
    \item y: Log PSA}}
\usage{data(Prostate)}
\source{\url{https://web.stanford.edu/~hastie/ElemStatLearn/}}
\references{
  \itemize{
    \item Hastie T, Tibshirani R, and Friedman J. (2001). \emph{The
      Elements of Statistical Learning}.  Springer.
    \item Stamey T, et al. (1989).  Prostate specific antigen in the
    diagnosis and treatment of adenocarcinoma of the
    prostate. II. Radical prostatectomy treated patients. \emph{Journal 
      of Urology}, \strong{16}: 1076-1083.
  }
}
\keyword{datasets}
