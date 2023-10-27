
#' calculates parameter-optimized fit (including LL) using a dirichlet multinomial model given a matrix of data (k), a model matrix with the features (mm), and a concentration parameter.
#'
#' @param k A matrix with the data to fit a dirichlet multinomial model.
#' @param mm A model matrix with the features.
#' @param conc The concentration parameter that accounts for overdispersion. To estimate this value use the estimate_overdispersion() function.
#'
#' @return A list with beta = a matrix with the best set of parameters found; and ll = the log-likelihood.
#' @export
dirmnom_fit <- function(k, mm, conc) {
    fracs <- k / rowSums(k)
    beta0 <- rbind(apply(fracs, 2, function(y) glm.fit(mm, y, family = binomial(link = "logit"), weights = rowSums(k))$coefficients))
    ores <- optim(
        beta0[, -ncol(beta0)] + beta0[, ncol(beta0)],
        ll_ <- function(beta) {
            mm %*% matrix(as.vector(beta), nrow = ncol(mm)) -> eta
            cbind(exp(eta), 1) -> pihat
            pihat / rowSums(pihat) -> pihat
            sum(sapply(1:nrow(k), function(i) ddirmnom(k[i, ], pihat[i, ] * conc, log = TRUE)))

            -ll_
        }
    )
    return(list(beta = matrix(ores$par, nrow = ncol(mm)), ll = -ores$value))
}