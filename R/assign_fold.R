#' Assign folds for cross-validation
#' 
#' If y has only two unique values, fold assignments are chosen so that the
#' balance between outcomes is the same in each fold. This is useful for
#' logistic regression and time-to-event data (to balance the fraction of
#' observations that are censored).
#' 
#' @param y      Either (i) the vector of outcomes or (ii) a vector such as `1:n`
#' @param folds  Number of folds
#' @param seed   (optional) set a seed for reproducibility
#'   
#' @returns A vector of integers indicating fold assignments
#' 
#' @seealso `[cv.ncvreg()]`
#' 
#' @examples
#' assign_fold(rnorm(11), 2)
#' assign_fold(1:41, 7)
#' assign_fold(1:41, 7) |> table()
#' data(Heart)
#' assign_fold(Heart$y, 7) |> head()
#' assign_fold(Heart$y, 7) |> table(Heart$y)
#' @export assign_fold

assign_fold <- function(y, folds, seed) {
  if (!missing(seed)) {
    original_seed <- .GlobalEnv$.Random.seed
    on.exit(.GlobalEnv$.Random.seed <- original_seed)
    set.seed(seed)
  }
  n <- length(y)
  if (length(unique(y)) == 2) {
    u <- unique(y)
    ind1 <- which(y==u[1])
    ind2 <- which(y==u[2])
    n1 <- length(ind1)
    n2 <- length(ind2)
    fold1 <- 1:n1 %% folds
    fold2 <- (n1 + 1:n2) %% folds
    fold1[fold1 == 0] <- folds
    fold2[fold2 == 0] <- folds
    fold <- double(n)
    fold[y==u[1]] <- sample(fold1)
    fold[y==u[2]] <- sample(fold2)
  } else {
    fold <- sample(1:length(y) %% folds)
    fold[fold==0] <- folds
  }
  fold
}
