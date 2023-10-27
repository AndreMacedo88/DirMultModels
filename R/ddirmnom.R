
#' calculates the PDF of a dirichlet-multinomial compound distribution for the data x
#' @param x The
#' @param alpha The
#' @param log Weather
#'
#' @return The probability density value.
#' @export
ddirmnom <- function(x, alpha, log = FALSE) {
    n <- sum(x)
    sumalpha <- sum(alpha)
    if (!log) {
        gamma(sumalpha) * gamma(n + 1) / gamma(n + sumalpha) *
            prod(gamma(x + alpha) / gamma(alpha) / gamma(x + 1))
    } else {
        lgamma(sumalpha) + lgamma(n + 1) - lgamma(n + sumalpha) +
            sum(lgamma(x + alpha) - lgamma(alpha) - lgamma(x + 1))
    }
}