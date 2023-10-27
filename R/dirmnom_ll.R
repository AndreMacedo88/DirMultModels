

#' calculates a parameter-non-optimized log-likelihood using a dirichlet multinomial model given a matrix of data (k), a model matrix with the features (mm), and a concentration parameter.
#' This is useful for data sets that fail to converge the optimization process. e.g. those with not so many samples. 
#'
#' @param k A matrix with the data to fit a dirichlet multinomial model.
#' @param mm A model matrix with the features.
#' @param conc The concentration parameter that accounts for overdispersion. To estimate this value use the estimate_overdispersion() function.
#'
#' @return The log-likelihood for the fit.
#' @export
rough_dirmnom_ll <- function(k, mm, conc) {
    # we need to calculate the probabilities of usage of each column (UTR site) based on pairwise multinomial fits to the data
    fracs <- k / rowSums(k)
    beta0 <- rbind(apply(fracs, 2, function(y) glm.fit(mm, y, family = binomial(link = "logit"), weights = rowSums(k))$coefficients))
    # we apply a softmax function to get the probabilities
    eta <- exp(mm %*% beta0)
    pihat <- eta / rowSums(eta)
    # finally, we get the log likelihood of fitting a dirichlet multinomial to the data by summing up the likelihoods of all samples
    # the sapply will apply the ddirmnom function to each sample (row), using the counts (k) and estimated probabilities (pihat) multiplied by a concentration factor (conc) that accounts for overdispersion
    return(sum(sapply(1:nrow(k), function(i) ddirmnom(k[i, ], pihat[i, ] * conc, log = TRUE))))
}