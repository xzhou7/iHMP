
gpd_permutation_test <- function(x0, y, minM_epdf = 10,
                                 Nexc = 250, Nexc_shrink = 10, Nexc_alpha = 0.05,
                                 yfun, ystart = 200, ygrow = 100, ymax = 5000) {

    # Implementation of the GPD-based p-value estimation algorithm from
    # Knijnenburg, Wessels, Reinders, and Shmulevich (2009) Fewer
    # permutations, more accurate P-values. Bioinformatics 25(12): i161–i168.
    # DOI: 10.1093/bioinformatics/btp211

    if (missing(y)) {
        # Automatic sampling of y
        stopifnot(!missing(yfun))
        y <- sapply(rep(NA, ystart), function(i)yfun())

        while (sum(y > x0) < minM_epdf && length(y) < ymax) {
            y <- c(y, sapply(rep(NA, ygrow), function(i)yfun()))
        }
    }

    # Use the ecdf when there are at least 10 (default) exceedances
    M <- sum(y > x0)
    if (M >= minM_epdf) {
        return (M / length(y))
    }

    # Not enough exceedances to be confident about the p-value.. use the GPD
    library("extRemes")
    while (T) {
        # Fit a GPD
        t <- mean(y[max(ceiling(length(y)/2), length(y) - Nexc) + c(0, 1)])
        fit <- fevd(y, threshold = t, type="GP")

        # Test goodness-of-fit
        gof <- ks.test(pextRemes(fit, y[y>t]), punif)
        if (Nexc <= 50 || gof$p.value >= Nexc_alpha) {
            if (Nexc <= 50) {
                warning("GPD fit to tail samples fits poorly even with (<= 50) null samples. Consider adding more null samples.")
            }
            # Didn't reject.. fit is good enough
            break
        } else {
            Nexc <- Nexc - Nexc_shrink
        }
    }

    # Get the p-value estimate from the fit
    p <- mean(y>t) * pextRemes(fit, x0, lower.tail=F)

    return (p)
}

