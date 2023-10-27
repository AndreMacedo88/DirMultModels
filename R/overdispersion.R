#!/usr/bin/env Rscript

# For the theory behind this estimator for overdispersion, refer to the thesis by Afroz, 2018, page 35 for example, with changes for multinomial as explained in page 57.

#' Estimates overdispersion in multinomial models
#' 
#' This is an adaptation of code found spread out in the package mclogit.
#' The methods for estimating overdispersion are described in Afroz 2018. Summarized in page 54 and 55.
#'
#' @param y The vector of multinomial observations.
#' @param fitted_values The fitted mean values.
#' @param coefs The Coefficients or the number of coefficients.
#' @param method The method of estimating overdispersion. One of "Afroz", "Pearson", "Fletcher", or "Deviance". "Afroz" is recommended.
#' @param number_samples The number of samples that were fitted. Usually the number of rows in the matrix fitted (\code{y}). Only necessary if providing a flattened \code{y}.
#' 
#' @return The estimation for overdispersion.
#' @export
estimate_overdispersion <- function(y, fitted_values, coefs, method = "Afroz", number_samples = NULL, weight = NULL) {

    # prepare list of weights
    if (!is.null(weight)) {
        if (length(dim(weight)) > 1) {
            if (length(dim(y)) > 1) {
                w <- as.vector(t(rep(weight, each = ncol(y))))
            } else {
                y_ <- matrix(y, nrow = number_samples, byrow = TRUE)
                w <- as.vector(t(rep(weight, each = ncol(y_))))
            }
        } else {
            w <- weight
        }
    } else {
        w <- t(rowSums(y))
        w <- as.vector(t(rep(w, each = ncol(y))))
    }

    # prepare observations and number of samples
    if (length(dim(y) > 1)) {
        n <- nrow(y)
        y <- as.vector(t(sweep(y, 1, rowSums(y), FUN = "/")))
    } else {
        n <- number_samples
    }

    N <- length(w)

    if (length(coefs) == 1) {
        p <- coefs
    } else {
        p <- length(coefs)
    }

    if (length(dim(fitted_values) > 1)) {
        fitted_values <- as.vector(t(fitted_values))
    }

    res.df <- N - n - p
    if (method == "Deviance") {
        Dresid <- 2 * w * y * (log(y) - log(fitted_values))
        Dresid[w == 0 | y == 0] <- 0
        D <- sum(Dresid)
        phi <- D / res.df
    } else {
        X2 <- sum(w * (y - fitted_values)^2 / fitted_values)
        phi.pearson <- X2 / (N - n - p)
        if (method %in% c("Afroz", "Fletcher")) {
            s.bar <- sum((y - fitted_values) / fitted_values) / (N - n)
        }
        phi <- switch(method,
            Pearson = phi.pearson,
            Afroz = phi.pearson / (1 + s.bar),
            Fletcher = phi.pearson - (N - n) * s.bar / (N - n - p)
        )
    }
    return(phi)
}