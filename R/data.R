#' Risk factors associated with heart disease
#' 
#' Data from a subset of the Coronary Risk-Factor Study baseline survey,
#' carried out in rural South Africa.
#' 
#' @format A list of two objects: `y` and `X`
#' \describe{
#'   \item{y}{Coronary heart disease at baseline; 1=Yes 0=No}
#'   \item{X}{A matrix with 462 observations (rows) and 9 predictor variables
#'   (columns). The remainder of this list describes the columns of `X`}
#'   \item{sbp}{Systolic blood pressure}
#'   \item{tobacco}{Cumulative tobacco consumption, in kg}
#'   \item{ldl}{Low-density lipoprotein cholesterol}
#'   \item{adiposity}{Adipose tissue concentration}
#'   \item{famhist}{Family history of heart disease (1=Present, 0=Absent)}
#'   \item{typea}{Score on test designed to measure type-A behavior}
#'   \item{obesity}{Obesity}
#'   \item{alcohol}{Current consumption of alcohol}
#'   \item{age}{Age of subject}
#' }
#'
#' @aliases heart
#' @references
#' \itemize{
#'   \item Hastie T, Tibshirani R, and Friedman J. (2001). *The Elements of
#'   Statistical Learning*.  Springer.
#'   \item Rousseauw J, et al. (1983). Coronary risk factor screening in three
#'   rural communities. *South African Medical Journal*, **64**: 430-436.
#' }
#' @source \url{https://web.stanford.edu/~hastie/ElemStatLearn/}
"Heart"

#' VA lung cancer data set
#' 
#' Data from a randomised trial of two treatment regimens for lung cancer. This
#' is a standard survival analysis data set from the classic textbook by
#' Kalbfleisch and Prentice.
#' 
#' @format A list of two objects: `y` and `X`
#' \describe{
#'   \item{y}{A two column matrix (`Surv` object) containing the follow-up
#'   time (in days) and an indicator variable for whether the patient died
#'   while on the study or not.}
#'   \item{X}{A matrix with 137 observations (rows) and 9 predictor variables
#'   (columns). The remainder of this list describes the columns of `X`}
#'   \item{trt}{Treatment indicator (1=control group, 2=treatment group)}
#'   \item{karno}{Karnofsky performance score (0=bad, 100=good)}
#'   \item{diagtime}{Time from diagnosis to randomization (months)}
#'   \item{age}{Age (years, at baseline)}
#'   \item{prior}{Prior therapy (0=no, 1=yes)}
#'   \item{squamous}{Indicator for whether the cancer type is squamous cell
#'   carcinoma (0=no, 1=yes)}
#'   \item{small}{Indicator for whether the cancer type is small cell lung
#'   cancer (0=no, 1=yes)}
#'   \item{adeno}{Indicator for whether the cancer type is adenocarcinoma
#'   (0=no, 1=yes)}
#'   \item{large}{Indicator for whether the cancer type is large cell carcinoma
#'   (0=no, 1=yes)}
#' }
#' 
#' @seealso `ncvsurv()`
#' @references
#' \itemize{
#'   \item Kalbfleisch D and Prentice RL (1980), *The Statistical Analysis of
#'   Failure Time Data*. Wiley, New York.
#' }
#' @source \url{https://cran.r-project.org/package=survival}
"Lung"

#' Factors associated with prostate specific antigen
#' 
#' Data from a study by by Stamey et al. (1989) to examine the association
#' between prostate specific antigen (PSA) and several clinical measures that
#' are potentially associated with PSA in men who were about to receive a
#' radical prostatectomy.
#' 
#' @format A list of two objects: `y` and `X`
#' \describe{
#'   \item{y}{Log PSA}
#'   \item{X}{A matrix with 97 instances (rows) and 8 predictor variables
#'   (columns). The remainder of this list describes the columns of `X`}
#'   \item{lcavol}{Log cancer volume}
#'   \item{lweight}{Log prostate weight}
#'   \item{age}{The man's age (years)}
#'   \item{lbph}{Log of the amount of benign hyperplasia}
#'   \item{svi}{Seminal vesicle invasion (1=Yes, 0=No)}
#'   \item{lcp}{Log of capsular penetration}
#'   \item{gleason}{Gleason score}
#'   \item{pgg45}{Percent of Gleason scores 4 or 5}
#' }
#' 
#' @aliases prostate
#' @references
#' \itemize{
#'   \item Hastie T, Tibshirani R, and Friedman J. (2001). *The Elements of
#'   Statistical Learning*.  Springer.
#'   \item Stamey T, et al. (1989). Prostate specific antigen in the diagnosis
#'   and treatment of adenocarcinoma of the prostate. II. Radical prostatectomy
#'   treated patients. *Journal of Urology*, **16**: 1076-1083.
#' }
#' @source \url{https://web.stanford.edu/~hastie/ElemStatLearn/}
"Prostate"
