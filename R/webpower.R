print.webpower <- function(x, ...) {
    cat(x$method, "\n\n", sep = "")
    note <- x$note
    url <- x$url
    x[c("method", "note", "url", "alternative", "family", "parameter")] <- NULL
    x <- do.call(cbind, x)
    rownames(x) <- rep("   ", nrow(x))
    print(x)
    if (!is.null(note)) 
        cat("\nNOTE: ", note, sep = "")
    cat("\nURL: ", url,  sep = "")
    cat("\n")
    invisible(x)
}

plot.webpower <- function(x, xvar = NULL, yvar = NULL, xlab = NULL, ylab = NULL, 
    ...) {
	wp <- x
    wp[c("method", "note", "url", "alternative", "family", "parameter")] <- NULL
    wp.matrix <- do.call(cbind, wp)
    check <- is.null(xvar) + is.null(yvar)
    if (check == 2) {
        if (!is.null(xlab)) 
            xlab <- xlab else xlab <- "Sample size"
        if (!is.null(ylab)) 
            ylab <- ylab else ylab <- "Power"
        plot(wp$n, wp$power, xlab = xlab, ylab = ylab, ...)
        lines(wp$n, wp$power, ...)
    } else if (check == 1) {
        stop("Both x and y need to be specified!")
    } else {
        if (!is.null(xlab)) 
            xlab <- xlab else xlab <- xvar
        if (!is.null(ylab)) 
            ylab <- ylab else ylab <- yvar
        plot(wp.matrix[, xvar], wp.matrix[, yvar], xlab = xlab, ylab = ylab, 
            ...)
        lines(wp.matrix[, xvar], wp.matrix[, yvar], ...)
    }
}


wp.anova <- function(k = NULL, n = NULL, f = NULL, alpha = 0.05, power = NULL, 
    type = c("overall", "two.sided", "greater", "less")) {
    type <- type[1]
    
    if (sum(sapply(list(k, n, f, power, alpha), is.null)) != 1) 
        stop("exactly one of k, n, f, power, and alpha must be NULL")
    if (!is.null(f) && min(f) < 0) 
        stop("f must be positive")
    if (!is.null(k) && min(k) < 2) 
        stop("number of groups must be at least 2")
    if (!is.null(n) && min(n) < 2) 
        stop("number of observations in each group must be at least 2")
    if (!is.null(alpha) && !is.numeric(alpha) || any(0 > alpha | alpha > 
        1)) 
        stop(sQuote("alpha"), " must be numeric in [0, 1]")
    if (!is.null(power) && !is.numeric(power) || any(0 > power | power > 
        1)) 
        stop(sQuote("power"), " must be numeric in [0, 1]")
    if (type == "overall") 
        p.body <- quote({
            lambda <- n * f^2
            pf(qf(alpha, k - 1, n - k, lower = FALSE), k - 1, n - k, lambda, 
                lower = FALSE)
        })
    if (type == "two.sided") 
        p.body <- quote({
            lambda <- n * f^2
            pf(qf(alpha, 1, n - k, lower = FALSE), 1, n - k, lambda, lower = FALSE)
        })
    if (type == "greater") 
        p.body <- quote({
            lambda <- sqrt(n) * f
            pt(qt(alpha, n - k, lower = FALSE), n - k, lambda, lower = FALSE)
        })
    if (type == "less") 
        p.body <- quote({
            lambda <- sqrt(n) * f
            pt(qt(alpha, n - k), n - k, lambda)
        })
    if (is.null(power)) 
        power <- eval(p.body) else if (is.null(k)) 
        k <- uniroot(function(k) eval(p.body) - power, c(2 + 1e-10, 100))$root else if (is.null(n)) 
        n <- uniroot(function(n) eval(p.body) - power, c(2 + k + 1e-10, 
            1e+05))$root else if (is.null(f)) 
        f <- uniroot(function(f) eval(p.body) - power, c(1e-07, 1e+07))$root else if (is.null(alpha)) 
        alpha <- uniroot(function(alpha) eval(p.body) - power, c(1e-10, 
            1 - 1e-10))$root else stop("internal error")
    NOTE <- "n is the total sample size"
    if (type == "overall") {
        NOTE <- paste(NOTE, " (overall)", sep = "")
    } else {
        NOTE <- paste(NOTE, " (contrast, ", type, ")", sep = "")
    }
    METHOD <- "Power for One-way ANOVA"
    URL <- "http://psychstat.org/anova"
    structure(list(k = k, n = n, f = f, alpha = alpha, power = power, note = NOTE, 
        method = METHOD, url = URL), class = "webpower")
}


wp.prop <- function(h = NULL, n1 = NULL, n2 = NULL, alpha = 0.05, power = NULL, 
    type = c("1p", "2p", "2p2n"), alternative = c("two.sided", "less", 
        "greater")) {
    type <- type[1]
    if (type == "1p") 
        return(wp.one.prop(h = h, n = n1, alpha = alpha, power = power, 
            alternative = alternative))
    if (type == "2p") 
        return(wp.two.prop(h = h, n = n1, alpha = alpha, power = power, 
            alternative = alternative))
    if (type == "2p2n") 
        return(wp.two.prop.two.n(h = h, n1 = n1, n2 = n2, alpha = alpha, 
            power = power, alternative = alternative))
}

wp.one.prop <- function(h = NULL, n = NULL, alpha = 0.05, power = NULL, 
    alternative = c("two.sided", "less", "greater")) {
    if (sum(sapply(list(h, n, power, alpha), is.null)) != 1) 
        stop("exactly one of h, n, power, and alpha must be NULL")
    if (!is.null(n) && min(n) < 1) 
        stop("number of observations in each group must be at least 1")
    if (!is.null(alpha) && !is.numeric(alpha) || any(0 > alpha | alpha > 
        1)) 
        stop(sQuote("alpha"), " must be numeric in [0, 1]")
    if (!is.null(power) && !is.numeric(power) || any(0 > power | power > 
        1)) 
        stop(sQuote("power"), " must be numeric in [0, 1]")
    alternative <- match.arg(alternative)
    tside <- switch(alternative, less = 1, two.sided = 2, greater = 3)
    if (tside == 2 && !is.null(h)) 
        h <- abs(h)
    if (tside == 2) {
        p.body <- quote({
            pnorm(qnorm(alpha/2, lower = FALSE) - h * sqrt(n), lower = FALSE) + 
                pnorm(qnorm(alpha/2, lower = TRUE) - h * sqrt(n), lower = TRUE)
        })
    }
    if (tside == 3) {
        p.body <- quote({
            pnorm(qnorm(alpha, lower = FALSE) - h * sqrt(n), lower = FALSE)
        })
    }
    if (tside == 1) {
        p.body <- quote({
            pnorm(qnorm(alpha, lower = TRUE) - h * sqrt(n), lower = TRUE)
        })
    }
    if (is.null(power)) 
        power <- eval(p.body) else if (is.null(h)) {
        if (tside == 2) {
            h <- uniroot(function(h) eval(p.body) - power, c(1e-10, 10))$root
        }
        if (tside == 1) {
            h <- uniroot(function(h) eval(p.body) - power, c(-10, 5))$root
        }
        if (tside == 3) {
            h <- uniroot(function(h) eval(p.body) - power, c(-5, 10))$root
        }
    } else if (is.null(n)) 
        n <- uniroot(function(n) eval(p.body) - power, c(2 + 1e-10, 1e+05))$root else if (is.null(alpha)) 
        alpha <- uniroot(function(alpha) eval(p.body) - power, c(1e-10, 
            1 - 1e-10))$root else stop("internal error")
    METHOD <- "Power for one-sample proportion test"
    NOTE <- NULL
    URL <- "http://psychstat.org/prop"
    structure(list(h = h, n = n, alpha = alpha, power = power, url = URL, 
        method = METHOD, note = NOTE), class = "webpower")
}

wp.two.prop <- function(h = NULL, n = NULL, alpha = 0.05, power = NULL, 
    alternative = c("two.sided", "less", "greater")) {
    if (sum(sapply(list(h, n, power, alpha), is.null)) != 1) 
        stop("exactly one of h, n, power, and alpha must be NULL")
    if (!is.null(n) && min(n) < 1) 
        stop("number of observations in each group must be at least 1")
    if (!is.null(alpha) && !is.numeric(alpha) || any(0 > alpha | alpha > 
        1)) 
        stop(sQuote("alpha"), " must be numeric in [0, 1]")
    if (!is.null(power) && !is.numeric(power) || any(0 > power | power > 
        1)) 
        stop(sQuote("power"), " must be numeric in [0, 1]")
    alternative <- match.arg(alternative)
    tside <- switch(alternative, less = 1, two.sided = 2, greater = 3)
    if (tside == 2 && !is.null(h)) 
        h <- abs(h)
    if (tside == 3) {
        p.body <- quote({
            pnorm(qnorm(alpha, lower = FALSE) - h * sqrt(n/2), lower = FALSE)
        })
    }
    if (tside == 2) {
        p.body <- quote({
            pnorm(qnorm(alpha/2, lower = FALSE) - h * sqrt(n/2), lower = FALSE) + 
                pnorm(qnorm(alpha/2, lower = TRUE) - h * sqrt(n/2), lower = TRUE)
        })
    }
    if (tside == 1) {
        p.body <- quote({
            pnorm(qnorm(alpha, lower = TRUE) - h * sqrt(n/2), lower = TRUE)
        })
    }
    if (is.null(power)) 
        power <- eval(p.body) else if (is.null(h)) {
        if (tside == 2) {
            h <- uniroot(function(h) eval(p.body) - power, c(1e-10, 10))$root
        }
        if (tside == 1) {
            h <- uniroot(function(h) eval(p.body) - power, c(-10, 5))$root
        }
        if (tside == 3) {
            h <- uniroot(function(h) eval(p.body) - power, c(-5, 10))$root
        }
    } else if (is.null(n)) 
        n <- uniroot(function(n) eval(p.body) - power, c(2 + 1e-10, 1e+05))$root else if (is.null(alpha)) 
        alpha <- uniroot(function(alpha) eval(p.body) - power, c(1e-10, 
            1 - 1e-10))$root else stop("internal error")
    
    NOTE <- "Sample sizes for EACH group"
    URL <- "http://psychstat.org/prop2p"
    METHOD <- "Power for two-sample proportion (equal n)"
    structure(list(h = h, n = n, alpha = alpha, power = power, url = URL, 
        method = METHOD, note = NOTE), class = "webpower")
}

wp.two.prop.two.n <- function(h = NULL, n1 = NULL, n2 = NULL, alpha = 0.05, 
    power = NULL, alternative = c("two.sided", "less", "greater")) {
    if (sum(sapply(list(h, n1, n2, power, alpha), is.null)) != 1) 
        stop("exactly one of h, n1, n2, power, and alpha must be NULL")
    if (!is.null(n1) && min(n1) < 2) 
        stop("number of observations in the first group must be at least 2")
    if (!is.null(n2) && min(n2) < 2) 
        stop("number of observations in the second group must be at least 2")
    if (!is.null(alpha) && !is.numeric(alpha) || any(0 > alpha | alpha > 
        1)) 
        stop(sQuote("alpha"), " must be numeric in [0, 1]")
    if (!is.null(power) && !is.numeric(power) || any(0 > power | power > 
        1)) 
        stop(sQuote("power"), " must be numeric in [0, 1]")
    alternative <- match.arg(alternative)
    tside <- switch(alternative, less = 1, two.sided = 2, greater = 3)
    if (tside == 2 && !is.null(h)) 
        h <- abs(h)
    if (tside == 3) {
        p.body <- quote({
            pnorm(qnorm(alpha, lower = FALSE) - h * sqrt((n1 * n2)/(n1 + 
                n2)), lower = FALSE)
        })
    }
    if (tside == 1) {
        p.body <- quote({
            pnorm(qnorm(alpha, lower = TRUE) - h * sqrt((n1 * n2)/(n1 + 
                n2)), lower = TRUE)
        })
    }
    if (tside == 2) {
        p.body <- quote({
            pnorm(qnorm(alpha/2, lower = FALSE) - h * sqrt((n1 * n2)/(n1 + 
                n2)), lower = FALSE) + pnorm(qnorm(alpha/2, lower = TRUE) - 
                h * sqrt((n1 * n2)/(n1 + n2)), lower = TRUE)
        })
    }
    if (is.null(power)) 
        power <- eval(p.body) else if (is.null(h)) {
        if (tside == 2) {
            h <- uniroot(function(h) eval(p.body) - power, c(1e-10, 10))$root
        }
        if (tside == 1) {
            h <- uniroot(function(h) eval(p.body) - power, c(-10, 5))$root
        }
        if (tside == 3) {
            h <- uniroot(function(h) eval(p.body) - power, c(-5, 10))$root
        }
    } else if (is.null(n1)) 
        n1 <- uniroot(function(n1) eval(p.body) - power, c(2 + 1e-10, 1e+05))$root else if (is.null(n2)) 
        n2 <- uniroot(function(n2) eval(p.body) - power, c(2 + 1e-10, 1e+05))$root else if (is.null(alpha)) 
        alpha <- uniroot(function(alpha) eval(p.body) - power, c(1e-10, 
            1 - 1e-10))$root else stop("internal error")
    NOTE <- "Sample size for each group"
    URL <- "http://psychstat.org/prop2p2n"
    METHOD <- "Power for two-sample proportion (unequal n)"
    structure(list(h = h, n1 = n1, n2 = n2, alpha = alpha, power = power, 
        url = URL, method = METHOD, note = NOTE), class = "webpower")
}

wp.t <- function(n1 = NULL, n2 = NULL, d = NULL, alpha = 0.05, power = NULL, 
    type = c("two.sample", "one.sample", "paired", "two.sample.2n"), alternative = c("two.sided", 
        "less", "greater"), tol = .Machine$double.eps^0.25) {
    type <- type[1]
    if (type == "two.sample.2n") {
        return(wp.t2(n1 = n1, n2 = n2, d = d, alpha = alpha, power = power, 
            alternative = alternative, tol = tol))
    } else {
        return(wp.t1(n = n1, d = d, alpha = alpha, power = power, type = type, 
            alternative = alternative, tol = tol))
    }
    
}

wp.t1 <- function(n = NULL, d = NULL, alpha = 0.05, power = NULL, type = c("two.sample", 
    "one.sample", "paired"), alternative = c("two.sided", "less", "greater"), 
    tol = .Machine$double.eps^0.25) {
    if (sum(sapply(list(n, d, power, alpha), is.null)) != 1) 
        stop("exactly one of n, d, power, and alpha must be NULL")
    if (!is.null(alpha) && !is.numeric(alpha) || any(0 > alpha | alpha > 
        1)) 
        stop(sQuote("alpha"), " must be numeric in [0, 1]")
    if (!is.null(power) && !is.numeric(power) || any(0 > power | power > 
        1)) 
        stop(sQuote("power"), " must be numeric in [0, 1]")
    type <- match.arg(type)
    alternative <- match.arg(alternative)
    tsample <- switch(type, one.sample = 1, two.sample = 2, paired = 1)
    ttside <- switch(alternative, less = 1, two.sided = 2, greater = 3)
    tside <- switch(alternative, less = 1, two.sided = 2, greater = 1)
    if (tside == 2 && !is.null(d)) 
        d <- abs(d)
    if (ttside == 1) {
        p.body <- quote({
            nu <- (n - 1) * tsample
            pt(qt(alpha/tside, nu, lower = TRUE), nu, ncp = sqrt(n/tsample) * 
                d, lower = TRUE)
        })
    }
    if (ttside == 2) {
        p.body <- quote({
            nu <- (n - 1) * tsample
            qu <- qt(alpha/tside, nu, lower = FALSE)
            pt(qu, nu, ncp = sqrt(n/tsample) * d, lower = FALSE) + pt(-qu, 
                nu, ncp = sqrt(n/tsample) * d, lower = TRUE)
        })
    }
    if (ttside == 3) {
        p.body <- quote({
            nu <- (n - 1) * tsample
            pt(qt(alpha/tside, nu, lower = FALSE), nu, ncp = sqrt(n/tsample) * 
                d, lower = FALSE)
        })
    }
    if (is.null(power)) 
        power <- eval(p.body) else if (is.null(n)) 
        n <- uniroot(function(n) eval(p.body) - power, c(2 + 1e-10, 1e+07), 
            tol = tol, extendInt = "upX")$root else if (is.null(d)) {
        if (ttside == 2) {
            d <- uniroot(function(d) eval(p.body) - power, c(1e-07, 10), 
                tol = tol, extendInt = "upX")$root
        }
        if (ttside == 1) {
            d <- uniroot(function(d) eval(p.body) - power, c(-10, 5), tol = tol, 
                extendInt = "upX")$root
        }
        if (ttside == 3) {
            d <- uniroot(function(d) eval(p.body) - power, c(-5, 10), tol = tol, 
                extendInt = "upX")$root
        }
    } else if (is.null(alpha)) 
        alpha <- uniroot(function(alpha) eval(p.body) - power, c(1e-10, 
            1 - 1e-10), tol = tol, extendInt = "upX")$root else stop("internal error", domain = NA)
    
    NOTE <- switch(type, paired = "n is number of *pairs*", two.sample = "n is number in *each* group", 
        NULL)
    METHOD <- paste(switch(type, one.sample = "One-sample", two.sample = "Two-sample", 
        paired = "Paired"), "t-test")
    URL <- "http://psychstat.org/ttest"
    structure(list(n = n, d = d, alpha = alpha, power = power, alternative = alternative, 
        note = NOTE, method = METHOD, url = URL), class = "webpower")
}

wp.t2 <- function(n1 = NULL, n2 = NULL, d = NULL, alpha = 0.05, power = NULL, 
    alternative = c("two.sided", "less", "greater"), tol = .Machine$double.eps^0.25) {
    if (sum(sapply(list(n1, n2, d, power, alpha), is.null)) != 1) 
        stop("exactly one of n1, n2, d, power, and alpha must be NULL")
    if (!is.null(alpha) && !is.numeric(alpha) || any(0 > alpha | alpha > 
        1)) 
        stop(sQuote("alpha"), " must be numeric in [0, 1]")
    if (!is.null(power) && !is.numeric(power) || any(0 > power | power > 
        1)) 
        stop(sQuote("power"), " must be numeric in [0, 1]")
    if (!is.null(n1) && min(n1) < 2) 
        stop("number of observations in the first group must be at least 2")
    if (!is.null(n2) && min(n2) < 2) 
        stop("number of observations in the second group must be at least 2")
    alternative <- match.arg(alternative)
    tsample <- 2
    ttside <- switch(alternative, less = 1, two.sided = 2, greater = 3)
    tside <- switch(alternative, less = 1, two.sided = 2, greater = 1)
    if (tside == 2 && !is.null(d)) 
        d <- abs(d)
    if (ttside == 1) {
        p.body <- quote({
            nu <- n1 + n2 - 2
            pt(qt(alpha/tside, nu, lower = TRUE), nu, ncp = d * (1/sqrt(1/n1 + 
                1/n2)), lower = TRUE)
        })
    }
    if (ttside == 2) {
        p.body <- quote({
            nu <- n1 + n2 - 2
            qu <- qt(alpha/tside, nu, lower = FALSE)
            pt(qu, nu, ncp = d * (1/sqrt(1/n1 + 1/n2)), lower = FALSE) + 
                pt(-qu, nu, ncp = d * (1/sqrt(1/n1 + 1/n2)), lower = TRUE)
        })
    }
    if (ttside == 3) {
        p.body <- quote({
            nu <- n1 + n2 - 2
            pt(qt(alpha/tside, nu, lower = FALSE), nu, ncp = d * (1/sqrt(1/n1 + 
                1/n2)), lower = FALSE)
        })
    }
    if (is.null(power)) 
        power <- eval(p.body) else if (is.null(n1)) 
        n1 <- uniroot(function(n1) eval(p.body) - power, c(2 + 1e-10, 1e+07), 
            tol = tol, extendInt = "upX")$root else if (is.null(n2)) 
        n2 <- uniroot(function(n2) eval(p.body) - power, c(2 + 1e-10, 1e+07), 
            tol = tol, extendInt = "upX")$root else if (is.null(d)) {
        if (ttside == 2) {
            d <- uniroot(function(d) eval(p.body) - power, c(1e-07, 10), 
                tol = tol, extendInt = "upX")$root
        }
        if (ttside == 1) {
            d <- uniroot(function(d) eval(p.body) - power, c(-10, 5), tol = tol, 
                extendInt = "upX")$root
        }
        if (ttside == 3) {
            d <- uniroot(function(d) eval(p.body) - power, c(-5, 10), tol = tol, 
                extendInt = "upX")$root
        }
    } else if (is.null(alpha)) 
        alpha <- uniroot(function(alpha) eval(p.body) - power, c(1e-10, 
            1 - 1e-10), tol = tol, extendInt = "upX")$root else stop("internal error", domain = NA)
    NOTE <- "n1 and n2 are number in *each* group"
    METHOD <- "Unbalanced two-sample t-test"
    URL <- "http://psychstat.org/ttest2n"
    structure(list(n1 = n1, n2 = n2, d = d, alpha = alpha, power = power, 
        alternative = alternative, note = NOTE, method = METHOD, url = URL), 
        class = "webpower")
}

## function for repeated measures anova analysis
wp.rmanova <- function(n = NULL, ng = NULL, nm = NULL, f = NULL, nscor = 1, 
    alpha = 0.05, power = NULL, type = 0) {
    if (sum(sapply(list(n, ng, nm, f, nscor, power, alpha), is.null)) != 
        1) 
        stop("exactly one of n, ng, nm, power, and alpha must be NULL")
    if (!is.null(f) && min(f) < 0) 
        stop("Effect size must be positive")
    if (!is.null(n) && min(n) < 1) 
        stop("Sample size has to be larger than 1")
    if (!is.null(ng) && min(ng) < 1) 
        stop("Number of groups have to be at least 1")
    if (!is.null(nm) && min(nm) < 2) 
        stop("number of measurements must be at least 1")
    if (!is.null(nscor) && !is.numeric(nscor) || any(0 > nscor | nscor > 
        1)) 
        stop("Nonsphericity correction must be numeric in [0, 1]")
    if (!is.null(alpha) && !is.numeric(alpha) || any(0 > alpha | alpha > 
        1)) 
        stop(sQuote("alpha"), " must be numeric in [0, 1]")
    if (!is.null(power) && !is.numeric(power) || any(0 > power | power > 
        1)) 
        stop(sQuote("power"), " must be numeric in [0, 1]")
    ## function to evaluate
    if (type == 0) {
        p.body <- quote({
            df1 <- ng - 1
            df2 <- n - ng
            lambda <- f^2 * n * nscor
            pf(qf(alpha, df1, df2, lower = FALSE), df1, df2, lambda, lower = FALSE)
        })
    }
    if (type == 1) {
        p.body <- quote({
            df1 <- (nm - 1) * nscor
            df2 <- (n - ng) * df1
            lambda <- f^2 * n * nscor
            pf(qf(alpha, df1, df2, lower = FALSE), df1, df2, lambda, lower = FALSE)
        })
    }
    if (type == 2) {
        p.body <- quote({
            df1 <- (ng - 1) * (nm - 1) * nscor
            df2 <- (n - ng) * (nm - 1) * nscor
            lambda <- f^2 * n * nscor
            pf(qf(alpha, df1, df2, lower = FALSE), df1, df2, lambda, lower = FALSE)
        })
    }
    if (is.null(power)) 
        power <- eval(p.body) else if (is.null(n)) {
        n <- uniroot(function(n) eval(p.body) - power, c(5 + ng, 1e+07))$root
    } else if (is.null(ng)) {
        ng <- uniroot(function(ng) eval(p.body) - power, c(1 + 1e-10, 
            1e+05))$root
    } else if (is.null(nm)) {
        nm <- uniroot(function(nm) eval(p.body) - power, c(1e-07, 1e+07))$root
    } else if (is.null(f)) {
        f <- uniroot(function(f) eval(p.body) - power, c(1e-07, 1e+07))$root
    } else if (is.null(alpha)) {
        alpha <- uniroot(function(alpha) eval(p.body) - power, c(1e-10, 
            1 - 1e-10))$root
    } else stop("internal error")
    
    if (type == 0) 
        NOTE <- "Power analysis for between-effect test"
    if (type == 1) 
        NOTE <- "Power analysis for within-effect test"
    if (type == 2) 
        NOTE <- "Power analysis for interaction-effect test"
    
    URL <- "http://psychstat.org/rmanova"
    METHOD <- "Repeated-measures ANOVA analysis"
    structure(list(n = n, f = f, ng = ng, nm = nm, nscor = nscor, alpha = alpha, 
        power = power, note = NOTE, method = METHOD, url = URL), class = "webpower")
}

## Effect size and ICC calculation for CRT with 2 arms based on
## empirical data The data has to be in text format where the first
## column of the data is the ID variable, the second column represents
## cluster, the third column is the outcome variable, and the fourth
## column is the condition variable (0 for control, 1 for condition).
## The first line of the data should be the variable names.
wp.effect.CRT2arm <- function(file) {
    if (is.character(file)){
		dat <- read.table(file, header = TRUE)
	}else{
		dat <- file
	}
    
    J <- length(unique(dat[, 2]))
    n <- nrow(dat)/J
    
    dat.trt <- dat[dat[, 4] == 1, ]
    dat.ctl <- dat[dat[, 4] == 0, ]
    
    mutj <- aggregate(dat.trt[, 3], by = list(dat.trt[, 2]), FUN = mean)
    mucj <- aggregate(dat.ctl[, 3], by = list(dat.ctl[, 2]), FUN = mean)
    
    mu_j <- rbind(mutj, mucj)
    names(mu_j) <- c(names(dat)[2], "mu")
    dat.mu <- merge(dat, mu_j, by = names(dat)[2])
    SSW <- sum((dat.mu[, 3] - dat.mu[5])^2)/(n * J - J)
    
    mut <- mean(mutj[, 2])
    muc <- mean(mucj[, 2])
    SSB <- (sum((mutj[, 2] - mut)^2) + sum((mucj[, 2] - muc)^2))/(J - 2)
    
    f <- (mut - muc)/sqrt(SSB + (1 - 1/n) * SSW)
    rho <- (SSB - SSW/n)/(SSB + (1 - 1/n) * SSW)
    
    return(list(f = f, ICC = rho))
}


## Effect size and ICC calculation for CRT with 3 arms based on
## empirical data The data has to be in text format, where the first
## column of the data is the ID variable, the second column represents
## cluster, the third column is the outcome variable, and the fourth
## column is the condition variable (0 for control, 1 for treatment1, 2
## for treatment2).  The first line of the data should be the variable
## names.
wp.effect.CRT3arm <- function(file) {
    if (is.character(file)){
		dat <- read.table(file, header = TRUE)
	}else{
		dat <- file
	}
	 
    J <- length(unique(dat[, 2]))
    n <- nrow(dat)/J
    
    dat.ctl <- dat[dat[, 4] == 0, ]
    dat.trt1 <- dat[dat[, 4] == 1, ]
    dat.trt2 <- dat[dat[, 4] == 2, ]
    
    mut1j <- aggregate(dat.trt1[, 3], by = list(dat.trt1[, 2]), FUN = mean)
    mut2j <- aggregate(dat.trt2[, 3], by = list(dat.trt2[, 2]), FUN = mean)
    mucj <- aggregate(dat.ctl[, 3], by = list(dat.ctl[, 2]), FUN = mean)
    
    mu_j <- rbind(mut1j, mut2j, mucj)
    names(mu_j) <- c(names(dat)[2], "mu")
    dat.mu <- merge(dat, mu_j, by = names(dat)[2])
    SSW <- sum((dat.mu[, 3] - dat.mu[5])^2)/(n * J - J)
    
    mut1 <- mean(mut1j[, 2])
    mut2 <- mean(mut2j[, 2])
    muc <- mean(mucj[, 2])
    SSB <- (sum((mut1j[, 2] - mut1)^2) + sum((mut2j[, 2] - mut2)^2) + sum((mucj[, 
        2] - muc)^2))/(J - 3)
    tau_hat <- SSB - SSW/n
    
    f1 <- (0.5 * mut1 + 0.5 * mut2 - muc)/sqrt(tau_hat + SSW)
    f2 <- (mut1 - mut2)/sqrt(tau_hat + SSW)
    f3 <- sqrt((((0.5 * (mut1 + mut2) - muc)^2)/4.5 + ((mut1 - mut2)^2)/6)/(tau_hat + 
        SSW))
    rho <- tau_hat/(tau_hat + SSW)
    
    return(list(f1 = f1, f2 = f2, f3 = f3, ICC = rho))
    # f1: effect size of treatment main effect f2: effect size of comparing
    # two treatments f3: effect size of omnibus test
}



## Effect size for MRT with 2 arms based on empirical data The data has
## to be in text format where the first column of the data is the ID
## variable, the second column represents cluster, the third column is
## the outcome variable, and the fourth column is the condition variable
## (0 for control, 1 for condition).  The first line of the data should
## be the variable names.
wp.effect.MRT2arm <- function(file) {
    if (is.character(file)){
		dat <- read.table(file, header = TRUE)
	}else{
		dat <- file
	}
    J <- length(unique(dat[, 2]))
    n <- nrow(dat)/J
    
    mutcj <- aggregate(dat[, 3], by = list(dat[, 2], dat[, 4]), FUN = mean)
    
    mutj <- mutcj[mutcj[, 2] == 1, ]
    mucj <- mutcj[mutcj[, 2] == 0, ]
    
    beta0 <- aggregate(mutcj[, 3], by = list(mutcj[, 1]), FUN = mean)
    names(beta0) <- c(names(dat)[2], "beta0")
    dat1 <- merge(dat, beta0, by = names(dat)[2])  #merge beta0
    
    mut_c <- merge(mutj, mucj, names(mutj)[1])
    mut_c$beta1 <- mut_c[, 3] - mut_c[, 5]
    mut_c1 <- mut_c[, c(1, 6)]
    names(mut_c1) <- c(names(dat)[2], "beta1")
    dat2 <- merge(dat1, mut_c1, by = names(dat)[2])  #merge beta1
    
    mud <- sum(mut_c[, 3] - mut_c[, 5])/J
    Xij <- dat[, 4] - 0.5
    sg2 <- sum((dat2[, 3] - dat2[, 5] - dat2[, 6] * Xij)^2)/(J * (n - 2))
    
    f <- mud/sqrt(sg2)
    return(list(f = f))
}



## Effect size for MRT with 3 arms based on empirical data The data has
## to be in text format, where the first column of the data is the ID
## variable, the second column represents cluster, the third column is
## the outcome variable, and the fourth column is the condition variable
## (0 for control, 1 for treatment1, 2 for treatment2).  The first line
## of the data should be the variable names.
wp.effect.MRT3arm <- function(file) {
    if (is.character(file)){
		dat <- read.table(file, header = TRUE)
	}else{
		dat <- file
	}
    J <- length(unique(dat[, 2]))
    n <- nrow(dat)/J
    
    mutcj <- aggregate(dat[, 3], by = list(dat[, 2], dat[, 4]), FUN = mean)
    
    mut1j <- mutcj[mutcj[, 2] == 1, ]
    names(mut1j)[3] <- "y1"
    mut2j <- mutcj[mutcj[, 2] == 2, ]
    names(mut2j)[3] <- "y2"
    mucj <- mutcj[mutcj[, 2] == 0, ]
    names(mucj)[3] <- "y0"
    
    beta0 <- aggregate(mutcj[, 3], by = list(mutcj[, 1]), FUN = mean)
    names(beta0) <- c(names(dat)[2], "beta0")
    dat1 <- merge(dat, beta0, by = names(dat)[2])  #merge beta0
    
    
    mut_c1 <- merge(mut1j, mucj, names(mut1j)[1])
    mut_c <- merge(mut2j, mut_c1, names(mut2j)[1])
    mut_c <- mut_c[, c(1, 3, 5, 7)]
    
    mut_c$beta1 <- (mut_c$y1 + mut_c$y2)/2 - mut_c$y0
    mut_c$beta2 <- (mut_c$y1 - mut_c$y2)
    
    
    names(mut_c)[1] <- names(dat)[2]
    dat2 <- merge(dat1, mut_c, by = names(dat)[2])  #merge beta1 and beta2
    
    mud1 <- sum(mut_c$y1 - mut_c$y0)/J
    mud2 <- sum(mut_c$y2 - mut_c$y0)/J
    
    X1ij <- rep(0, nrow(dat))  # 0.5 for trt1 and trt2, -1 for ctl
    X2ij <- rep(0, nrow(dat))  # 1 for trt1, -1 for trt2, 0 for ctl
    X1ij[dat[, 4] == 1] <- 0.5
    X1ij[dat[, 4] == 2] <- 0.5
    X1ij[dat[, 4] == 0] <- -1
    X2ij[dat[, 4] == 1] <- 1
    X2ij[dat[, 4] == 2] <- -1
    X2ij[dat[, 4] == 0] <- 0
    
    sg2 <- sum((dat2[, 3] - dat2$beta0 - dat2$beta1 * X1ij - dat2$beta2 * 
        X2ij)^2)/(J * (n - 3))
    
    f1 <- (mud1 + mud2)/2/sqrt(sg2)  #effect size of treatment main effect
    f2 <- (mud1 - mud2)/sqrt(sg2)  #effect size of comparing the two treatments
    return(list(f1 = f1, f2 = f2))
}

nuniroot <- function(f, interval, maxlength = 100) {
    if (length(interval)!=2) stop("Please provide an interval with two values such as c(0,1).")
	x <- seq(min(interval), max(interval), length = maxlength)
    f.out <- f(x)
    if (min(f.out) * max(f.out) > 0) 
        stop("The specified parameters do not yield valid results. Please try to supply a different interval, e.g., using interval=c(0,1), for your parameter.") else {
        low <- max(f.out[f.out < 0])
        high <- min(f.out[f.out > 0])
        interval <- c(x[f.out == low][1], x[f.out == high][1])
        uniroot(f, interval)
    }
}

####################################################### function for cluster randomized trials with 2 arms##
wp.crt2arm <- function(n = NULL, f = NULL, J = NULL, icc = NULL, power = NULL, 
    alpha = 0.05, alternative = c("two.sided", "one.sided"), interval = NULL) {
    alternative <- alternative[1]
    if (sum(sapply(list(J, n, f, icc, alpha, power), is.null)) != 1) 
        stop("exactly one of J, n, f, icc, and alpha,power must be NULL")
    if (!is.null(f) && min(f) < 0) 
        stop("Effect size must be positive")
    if (!is.null(n) && min(n) < 1) 
        stop("Sample size has to be larger than 1")
    if (!is.null(J) && min(J) < 2) 
        stop("Number of clusters has to be at least 2")
    if (!is.null(icc) && !is.numeric(icc) || any(0 > icc | icc > 1)) 
        stop("Intraclass correlation must be numeric in [0, 1]")
    if (!is.null(alpha) && !is.numeric(alpha) || any(0 > alpha | alpha > 
        1)) 
        stop("alpha must be numeric in [0, 1]")
    if (!is.null(power) && !is.numeric(power) || any(0 > power | power > 
        1)) 
        stop("Power must be numeric in [0, 1]")
    
    ## function to evaluate two-sided
    if (alternative == "two.sided") {
        p.body <- quote({
            df <- J - 2
            lambda <- sqrt(J * f^2/(4 * icc + 4 * (1 - icc)/n))
            1 - pt(qt(1 - alpha/2, df), df, lambda) + pt(-qt(1 - alpha/2, 
                df), df, lambda)
        })
    } else {
        ## one-sided
        p.body <- quote({
            df <- J - 2
            lambda <- sqrt(J * f^2/(4 * icc + 4 * (1 - icc)/n))
            1 - pt(qt(1 - alpha, J - 2), df, lambda)
        })
    }
    
    #### 
    if (is.null(power)) 
        power <- eval(p.body) else if (is.null(J)) 
        J <- nuniroot(function(J) eval(p.body) - power, ifelse(rep(is.null(interval), 2), c(2 + 1e-10, 1000), interval))$root else if (is.null(n)) 
        n <- nuniroot(function(n) eval(p.body) - power, ifelse(rep(is.null(interval), 2), c(1, 1e+06), interval))$root else if (is.null(f)) 
        f <- nuniroot(function(f) eval(p.body) - power, ifelse(rep(is.null(interval), 2), c(1e-07, 1e+07), interval))$root else if (is.null(icc)) 
        icc <- nuniroot(function(icc) eval(p.body) - power, ifelse(rep(is.null(interval), 2), c(0, 1), interval))$root else if (is.null(alpha)) 
        alpha <- nuniroot(function(alpha) eval(p.body) - power, ifelse(rep(is.null(interval), 2), c(1e-10, 1 - 1e-10), interval))$root else stop("internal error")
    
    NOTE <- "n is the number of subjects per cluster."
    METHOD <- "Cluster randomized trials with 2 arms"
    URL <- "http://psychstat.org/crt2arm"
    structure(list(J = J, n = n, f = f, icc = icc, power = power, alpha = alpha, note = NOTE, method = METHOD, url = URL), class = "webpower")
}

####################################################### function for cluster randomized trials with 3 arms##
wp.crt3arm <- function(n = NULL, f = NULL, J = NULL, icc = NULL, power = NULL, 
    alpha = 0.05, alternative = c("two.sided", "one.sided"), type = c("main", 
        "treatment", "omnibus"), interval = NULL) {
    side <- alternative[1]
    type <- type[1]
    
    if (sum(sapply(list(J, n, f, icc, alpha, power), is.null)) != 1) 
        stop("exactly one of J, n, f, icc, and alpha,power must be NULL")
    if (!is.null(f) && min(f) < 0) 
        stop("Effect size must be positive")
    if (!is.null(n) && min(n) < 1) 
        stop("Sample size has to be larger than 1")
    if (!is.null(J) && min(J) < 3) 
        stop("Number of clusters has to be at least 3")
    if (!is.null(icc) && !is.numeric(icc) || any(0 > icc | icc > 1)) 
        stop("Intraclass correlation must be numeric in [0, 1]")
    if (!is.null(alpha) && !is.numeric(alpha) || any(0 > alpha | alpha > 
        1)) 
        stop("alpha must be numeric in [0, 1]")
    if (!is.null(power) && !is.numeric(power) || any(0 > power | power > 
        1)) 
        stop("Power must be numeric in [0, 1]")
    
    
    ## compare the average treatment with control
    if (type == "main") {
        if (side == "two.sided") {
            ## two-sided
            p.body <- quote({
                df <- J - 3
                lambda1 <- sqrt(J) * f/sqrt(4.5 * (icc + (1 - icc)/n))
                t0 <- qt(1 - alpha/2, df)
                1 - pt(t0, df, lambda1) + pt(-t0, df, lambda1)
            })
        } else {
            ## one-sided
            p.body <- quote({
                df <- J - 3
                lambda1 <- sqrt(J) * f/sqrt(4.5 * (icc + (1 - icc)/n))
                t0 <- qt(1 - alpha, df)
                1 - pt(t0, df, lambda1)
            })
        }
    }
    ## compare the two treatments
    if (type == "treatment") {
        if (side == "two.sided") {
            ## two-sided
            p.body <- quote({
                df <- J - 3
                lambda2 <- sqrt(J) * f/sqrt(6 * (icc + (1 - icc)/n))
                t0 <- qt(1 - alpha/2, df)
                1 - pt(t0, df, lambda2) + pt(-t0, df, lambda2)
            })
        } else {
            ## one-sided
            p.body <- quote({
                df <- J - 3
                lambda2 <- sqrt(J) * f/sqrt(6 * (icc + (1 - icc)/n))
                t0 <- qt(1 - alpha, df)
                1 - pt(t0, df, lambda2)
            })
        }
    }
    
    ## omnibus test
    if (type == "omnibus") {
        p.body <- quote({
            df1 <- 2
            df2 <- J - 3
            lambda3 <- J * f^2/(icc + (1 - icc)/n)
            f0 <- qf(1 - alpha, df1, df2)
            1 - pf(f0, df1, df2, lambda3)
        })
    }
    #### 
    if (is.null(power)) 
        power <- eval(p.body) else if (is.null(J)) 
        J <- nuniroot(function(J) eval(p.body) - power, ifelse(rep(is.null(interval), 2), c(3 + 1e-10, 1000), interval))$root else if (is.null(n)) 
        n <- nuniroot(function(n) eval(p.body) - power, ifelse(rep(is.null(interval), 2), c(2 + 1e-10, 1e+06), interval))$root else if (is.null(f)) 
        f <- nuniroot(function(f) eval(p.body) - power, ifelse(rep(is.null(interval), 2), c(1e-07, 1e+07), interval))$root else if (is.null(icc)) 
        icc <- nuniroot(function(icc) eval(p.body) - power, ifelse(rep(is.null(interval), 2), c(0, 1), interval))$root else if (is.null(alpha)) 
        alpha <- nuniroot(function(alpha) eval(p.body) - power, ifelse(rep(is.null(interval), 2), c(1e-10, 1 - 1e-10), interval))$root else stop("internal error")
    
    NOTE <- "n is the number of subjects per cluster."
    METHOD <- "Cluster randomized trials with 3 arms"
    URL <- "http://psychstat.org/crt3arm"
    structure(list(J = J, n = n, f = f, icc = icc, power = power, alpha = alpha, 
        note = NOTE, method = METHOD, url = URL), class = "webpower")
    
}

######################################################### function for multisite randomized trials with 2 arms##
wp.mrt2arm <- function(n = NULL, f = NULL, J = NULL, tau00 = NULL, tau11 = NULL, 
    sg2 = NULL, power = NULL, alpha = 0.05, alternative = c("two.sided", 
        "one.sided"), type = c("main", "site", "variance"), interval = NULL) {
    type <- type[1]
    side <- alternative[1]
    if (sum(sapply(list(f, n, J, power), is.null)) != 1 & type == "main") 
        stop("exactly one of f, J, n and power must be NULL")
    if (sum(sapply(list(n, J, power), is.null)) != 1 & type != "main") 
        stop("exactly one of J, n and power must be NULL")
    if (!is.null(f) && min(f) < 0) 
        stop("Effect size must be positive")
    if (!is.null(n) && min(n) < 2) 
        stop("Sample size has to be larger than 2")
    if (!is.null(J) && min(J) < 2) 
        stop("Number of sites has to be at least 2")
    if (!is.null(tau00) && min(tau00) < 0) 
        stop("Variance of site means must be positive")
    if (!is.null(tau11) && min(tau11) < 0) 
        stop("Variance of treatment main effects across sites must be positive")
    if ((!is.null(sg2) && min(sg2) < 0) || is.null(sg2)) 
        stop("Between-person variation must be a positive number")
    if (is.null(alpha) || (!is.null(alpha) && !is.numeric(alpha) || any(0 > 
        alpha | alpha > 1))) 
        stop("alpha must be numeric in [0, 1]")
    if (!is.null(power) && !is.numeric(power) || any(0 > power | power > 
        1)) 
        stop("Power must be numeric in [0, 1]")
    
    if (is.null(tau11) && (type == "main" || type == "variance")) 
        stop("For this type of test, variance of treatment main effects across sites must be specified")
    if (is.null(tau00) && type == "site") 
        stop("For this type of test, variance of site means must be specified")
    
    ## test treatment main effect #no tau00 needed
    if (type == "main") {
        if (side == "two.sided") {
            ## two-sided
            p.body <- quote({
                df <- J - 1
                lambda1 <- sqrt(J) * f/sqrt(4/n + tau11/sg2)
                t0 <- qt(1 - alpha/2, df)
                1 - pt(t0, df, lambda1) + pt(-t0, df, lambda1)
            })
        } else {
            ## one-sided #no tau00 needed
            p.body <- quote({
                df <- J - 1
                lambda1 <- sqrt(J) * f/sqrt(4/n + tau11/sg2)
                t0 <- qt(1 - alpha, df)
                1 - pt(t0, df, lambda1)
            })
        }
    }
    ## test site variability #no tau11 and f needed
    if (type == "site") {
        p.body <- quote({
            df1 <- J - 1
            df2 <- J * (n - 2)
            f0 <- qf(1 - alpha, df1, df2)
            1 - pf(f0/(n * tau00/sg2 + 1), df1, df2)
        })
    }
    
    ## test variance of treatment effects #no tau00 and f needed
    if (type == "variance") {
        p.body <- quote({
            df1 <- J - 1
            df2 <- J * (n - 2)
            f0 <- qf(1 - alpha, df1, df2)
            1 - pf(f0/(n * tau11/sg2/4 + 1), df1, df2)
        })
    }
    #### 
    if (is.null(power)) 
        power <- eval(p.body) else if (is.null(J)) 
        J <- nuniroot(function(J) eval(p.body) - power, ifelse(rep(is.null(interval), 2), c(1 + 1e-10, 1000), interval))$root else if (is.null(n)) 
        n <- nuniroot(function(n) eval(p.body) - power, ifelse(rep(is.null(interval), 2), c(3 - 1e-10, 1e+06), interval))$root else if (is.null(f) & type == 1) 
        f <- nuniroot(function(f) eval(p.body) - power, ifelse(rep(is.null(interval), 2), c(1e-07, 1e+07), interval))$root else stop("internal error")
    
    NOTE <- "n is the number of subjects per cluster"
    METHOD <- "Multisite randomized trials with 2 arms"
    URL <- "http://psychstat.org/mrt2arm"
    structure(list(J = J, n = n, f = f, tau00 = tau00, tau11 = tau11, sg2 = sg2, 
        power = power, alpha = alpha, note = NOTE, method = METHOD, url = URL), 
        class = "webpower")
    
}

######################################################### function for multisite randomized trials with 3 arms##
wp.mrt3arm <- function(n = NULL, f1 = NULL, f2 = NULL, J = NULL, tau = NULL, 
    sg2 = NULL, power = NULL, alpha = 0.05, alternative = c("two.sided", 
        "one.sided"), type = c("main", "treatment", "omnibus"), interval = NULL) {
    type <- type[1]
    side <- alternative[1]
    
    if (type == "main") 
        if (sum(sapply(list(f1, n, J, power, alpha), is.null)) != 1) 
            stop("exactly one of f1, J, n, power and alpha must be NULL")
    if (type == "treatment") 
        if (sum(sapply(list(f2, n, J, power, alpha), is.null)) != 1) 
            stop("exactly one of f2, J, n, power and alpha must be NULL")
    if (type == "omnibus") {
        if (sum(sapply(list(n, J, power, alpha), is.null)) != 1) 
            stop("exactly one of J, n, power and alpha must be NULL")
        if (is.null(f1) || is.null(f2)) 
            stop("f1 & f2 must be given")
    }
    
    
    if (!is.null(f1) && min(f1) < 0) 
        stop("Effect size must be positive")
    if (!is.null(f2) && min(f2) < 0) 
        stop("Effect size must be positive")
    if (!is.null(n) && min(n) < 3) 
        stop("Sample size has to be larger than 3")
    if (!is.null(J) && min(J) < 2) 
        stop("Number of sites has to be at least 2")
    if (!is.null(tau) && min(tau) < 0) 
        stop("Variance of treatment main effects across sites must be positive")
    if ((!is.null(sg2) && min(sg2) < 0) || is.null(sg2)) 
        stop("Between-person variation must be a positive number")
    if (!is.null(alpha) && !is.numeric(alpha) || any(0 > alpha | alpha > 
        1)) 
        stop("alpha must be numeric in [0, 1]")
    if (!is.null(power) && !is.numeric(power) || any(0 > power | power > 
        1)) 
        stop("Power must be numeric in [0, 1]")
    
    
    ## Test treatment main effect
    if (type == "main") {
        if (side == "two.sided") {
            ## two-sided
            p.body <- quote({
                df <- J - 1
                lambda1 <- sqrt(J) * f1/sqrt(4.5/n + 1.5 * tau/sg2)
                t0 <- qt(1 - alpha/2, df)
                1 - pt(t0, df, lambda1) + pt(-t0, df, lambda1)
            })
        } else {
            ## one-sided
            p.body <- quote({
                df <- J - 1
                lambda1 <- sqrt(J) * f1/sqrt(4.5/n + 1.5 * tau/sg2)
                t0 <- qt(1 - alpha, df)
                1 - pt(t0, df, lambda1)
            })
        }
    }
    ## Compare two treatments
    if (type == "treatment") {
        if (side == "two.sided") {
            ## two-sided
            p.body <- quote({
                df <- J - 1
                lambda2 <- sqrt(J) * f2/sqrt(6/n + 2 * tau/sg2)
                t0 <- qt(1 - alpha/2, df)
                1 - pt(t0, df, lambda2) + pt(-t0, df, lambda2)
            })
        } else {
            ## one-sided
            p.body <- quote({
                df <- J - 1
                lambda2 <- sqrt(J) * f2/sqrt(6/n + 2 * tau/sg2)
                t0 <- qt(1 - alpha, df)
                1 - pt(t0, df, lambda2)
            })
        }
    }
    
    ## Omnibus test
    if (type == "omnibus") {
        p.body <- quote({
            df1 <- 2
            df2 <- 2 * (J - 1)
            lambda1 <- sqrt(J) * f1/sqrt(4.5/n + 1.5 * tau/sg2)
            lambda2 <- sqrt(J) * f2/sqrt(6/n + 2 * tau/sg2)
            lambda3 <- lambda1^2 + lambda2^2
            f0 <- qf(1 - alpha, df1, df2)
            1 - pf(f0, df1, df2, lambda3)
        })
    }
    #### 
    if (is.null(power)) 
        power <- eval(p.body) else if (is.null(J)) 
        J <- nuniroot(function(J) eval(p.body) - power, ifelse(rep(is.null(interval), 2), c(2 - 1e-10, 1000), interval))$root else if (is.null(n)) 
        n <- nuniroot(function(n) eval(p.body) - power, ifelse(rep(is.null(interval), 2), c(3 - 1e-10, 1e+07), interval))$root else if (is.null(f1) && type != 2) 
        f1 <- nuniroot(function(f1) eval(p.body) - power, ifelse(rep(is.null(interval), 2), c(1e-07, 1e+07), interval))$root else if (is.null(f2) && type != 1) 
        f2 <- nuniroot(function(f2) eval(p.body) - power, ifelse(rep(is.null(interval), 2), c(1e-07, 1e+07), interval))$root else stop("internal error")
    
    NOTE <- "n is the number of subjects per cluster"
    METHOD <- "Multisite randomized trials with 3 arms"
    URL <- "http://psychstat.org/mrt3arm"
    structure(list(J = J, n = n, f1 = f1, f2 = f2, tau = tau, sg2 = sg2, 
        power = power, alpha = alpha, note = NOTE, method = METHOD, url = URL), 
        class = "webpower")
}

wp.regression <- function (n = NULL, p1 = NULL, p2 = 0, 
f2 = NULL, alpha = 0.05, power = NULL, 
type=c("regular", "Cohen")){
    if (sum(sapply(list(n, f2, power, alpha), is.null)) != 1) 
        stop("exactly one of n, f2, power, and alpha must be NULL")
    if (!is.null(f2) && min(f2) < 0) 
        stop("f2 must be positive")
    if (!is.null(n) && min(n) < 5) 
        stop("sample size has to be at least 5")
    if (!is.null(p1) && min(p1) < 1) 
        stop("number of predictor in the full model has to be larger than 1")
    if (!is.null(p2) && p1 < p2) 
        stop("number of predictor in the full model has to be larger than that in the reduced model")
    if (!is.null(alpha) && !is.numeric(alpha) || any(0 > alpha | 
        alpha > 1)) 
        stop(sQuote("alpha"), " must be numeric in [0, 1]")
    if (!is.null(power) && !is.numeric(power) || any(0 > power | 
        power > 1)) 
        stop(sQuote("power"), " must be numeric in [0, 1]")
    p.body <- quote({
        u <- p1 - p2
        v <- n - p1 - 1
        if (type[1]=="Cohen"){
			lambda <- f2 * (u + v + 1)
		}else{
			lambda <- f2 * n
			}
        pf(qf(alpha, u, v, lower = FALSE), u, v, lambda, lower = FALSE)
    })
    if (is.null(power)) 
        power <- eval(p.body)
    else if (is.null(n)) 
        n <- uniroot(function(n) eval(p.body) - power, c(5 + 
            p1 + 1e-10, 1e+05))$root
    else if (is.null(f2)) 
        f2 <- uniroot(function(f2) eval(p.body) - power, c(1e-07, 
            1e+07))$root
    else if (is.null(alpha)) 
        alpha <- uniroot(function(alpha) eval(p.body) - power, 
            c(1e-10, 1 - 1e-10))$root
    else stop("internal error")
    METHOD <- "Power for multiple regression"
    URL <- "http://psychstat.org/regression"
    structure(list(n = n, p1 = p1, p2 = p2, f2 = f2, alpha = alpha, 
        power = power, method = METHOD, url = URL), class = "webpower")
}

wp.logistic <- function(n = NULL, p0 = NULL, p1 = NULL, alpha = 0.05, power = NULL, 
    alternative = c("two.sided", "less", "greater"), family = c("Bernoulli", 
        "exponential", "lognormal", "normal", "Poisson", "uniform"), parameter = NULL) {
	sig.level <- alpha
    if (sum(sapply(list(n, power), is.null)) != 1) 
        stop("exactly one of n, power, and alpha must be NULL")
    
    if (is.null(p0) || !is.null(p0) && !is.numeric(p0) || any(0 > p0 | 
        p0 > 1)) 
        stop(sQuote("p0"), " must be numeric in (0,1)")
    
    if (is.null(p1) || !is.null(p1) && !is.numeric(p1) || any(0 > p1 | 
        p1 > 1)) 
        stop(sQuote("p1"), " must be numeric in (0,1)")
    
    if (!is.null(n) && min(n) < 5) 
        stop("number of observations must be at least 5")
    if (is.null(alpha) || !is.null(alpha) & !is.numeric(alpha) || any(alpha < 
        0 | alpha > 1)) 
        stop(sQuote("alpha"), " must be numeric in [0, 1]")
    if (!is.null(power) & !is.numeric(power) || any(0 > power | power > 
        1)) 
        stop(sQuote("power"), " must be numeric in [0, 1]")
    if (!(family %in% c("Bernoulli", 
        "exponential", "lognormal", "normal", "Poisson", "uniform"))) 
        stop ("family must be one of Bernoulli, exponential, lognormal, normal, Poisson, or uniform")
    
    alternative <- match.arg(alternative)
    tside <- switch(alternative, less = 1, two.sided = 2, greater = 3)
    
    
    p.body <- quote({
        s * pnorm(-qnorm(1 - alpha) - sqrt(n)/sqrt(g * v0 + (1 - g) * v1) * 
            beta1) + t * pnorm(-qnorm(1 - alpha) + sqrt(n)/sqrt(g * v0 + 
            (1 - g) * v1) * beta1)
    })
    
    
    if (family == "Bernoulli") {
        ## binomial predictor
        if (is.null(parameter)) {
            B <- 0.5
        } else {
            B <- parameter
        }
        
        g <- 0
        
        odds <- p1/(1 - p1)/(p0/(1 - p0))
        beta1 <- log(odds)  ##beta1 under H_A in the population 
        beta0 <- log(p0/(1 - p0))
        d <- B * p1 * (1 - p1) + (1 - B) * p0 * (1 - p0)
        e <- B * p1 * (1 - p1)
        f <- B * p1 * (1 - p1)
        v1 <- d/(d * f - e^2)
        ### v0
        mu1 <- B * p1 + (1 - B) * p0
        i00 <- log(mu1/(1 - mu1))
        pn <- 1/(1 + exp(-i00))
        a <- pn * (1 - pn)  #00
        b <- B * pn * (1 - pn)  #11
        c <- B * pn * (1 - pn)  #01
        v0 <- b/(a * b - c^2)
        if (tside == 1) {
            s <- 1
            t <- 0
            alpha <- alpha
        }
        if (tside == 2) {
            s <- 1
            t <- 1
            alpha <- alpha/2
        }
        if (tside == 3) {
            s <- 0
            t <- 1
            alpha <- alpha
        }
        if (is.null(power)) 
            power <- eval(p.body)
        if (is.null(n)) 
            n <- uniroot(function(n) eval(p.body) - power, c(2 + 1e-10, 
                1e+07))$root
    }
    
    if (family == "exponential") {
        ## Exponential predictor
        if (is.null(parameter)) {
            lambda <- 1
        } else {
            lambda <- parameter
        }
        
        g <- 0
        beta0 <- log(p0/(1 - p0))
        odds <- p1/(1 - p1)/(p0/(1 - p0))
        beta1 <- log(odds)  ##beta1 under H_A in the population##beta1 under H_A in the population
        d <- integrate(function(x) (1 - (1 + exp(-beta0 - beta1 * x))^(-1)) * 
            (1 + exp(-beta0 - beta1 * x))^(-1) * dexp(x, lambda), 0, Inf, 
            subdivisions = 100L)$value
        e <- integrate(function(x) x * (1 - (1 + exp(-beta0 - beta1 * x))^(-1)) * 
            (1 + exp(-beta0 - beta1 * x))^(-1) * dexp(x, lambda), 0, Inf, 
            subdivisions = 100L)$value
        f <- integrate(function(x) x^2 * (1 - (1 + exp(-beta0 - beta1 * 
            x))^(-1)) * (1 + exp(-beta0 - beta1 * x))^(-1) * dexp(x, lambda), 
            0, Inf, subdivisions = 100L)$value
        v1 <- d/(d * f - e^2)
        mu1 <- integrate(function(x) 1/(1 + exp(-beta0 - beta1 * x)) * 
            dexp(x, lambda), 0, Inf, subdivisions = 100L)$value
        
        i00 <- log(mu1/(1 - mu1))
        pn <- 1/(1 + exp(-i00))
        a <- pn * (1 - pn)  #00
        b <- 2 * lambda^(-2) * pn * (1 - pn)  #11
        c <- lambda^(-1) * pn * (1 - pn)  #01
        v0 <- a/(a * b - c^2)
        
        if (tside == 1) {
            s <- 1
            t <- 0
            alpha <- alpha
        }
        if (tside == 2) {
            s <- 1
            t <- 1
            alpha <- alpha/2
        }
        if (tside == 3) {
            s <- 0
            t <- 1
            alpha <- alpha
        }
        if (is.null(power)) 
            power <- eval(p.body)
        if (is.null(n)) 
            n <- uniroot(function(n) eval(p.body) - power, c(2 + 1e-10, 
                1e+07))$root
    }
    
    
    if (family == "lognormal") {
        
        g <- 0
        if (is.null(parameter)) {
            mu <- 0
            sigma <- 1
        } else {
            if (length(parameter) != 2) 
                stop("Both mean and standard deviation of the log-normal distribution have to be provided as a vector.")
            mu <- parameter[1]
            sigma <- parameter[2]
        }
        
        beta0 <- log(p0/(1 - p0))
        odds <- p1/(1 - p1)/(p0/(1 - p0))
        beta1 <- log(odds)  ##beta1 under H_A in the population
        d <- integrate(function(x) (1 - (1 + exp(-beta0 - beta1 * x))^(-1)) * 
            (1 + exp(-beta0 - beta1 * x))^(-1) * dlnorm(x, mu, sigma), 
            0, Inf, subdivisions = 100L)$value
        e <- integrate(function(x) x * (1 - (1 + exp(-beta0 - beta1 * x))^(-1)) * 
            (1 + exp(-beta0 - beta1 * x))^(-1) * dlnorm(x, mu, sigma), 
            0, Inf, subdivisions = 100L)$value
        f <- integrate(function(x) x^2 * (1 - (1 + exp(-beta0 - beta1 * 
            x))^(-1)) * (1 + exp(-beta0 - beta1 * x))^(-1) * dlnorm(x, 
            mu, sigma), 0, Inf, subdivisions = 100L)$value
        v1 <- d/(d * f - e^2)
        
        mu1 <- integrate(function(x) 1/(1 + exp(-beta0 - beta1 * x)) * 
            dlnorm(x, mu, sigma), 0, Inf, subdivisions = 100L)$value
        
        i00 <- log(mu1/(1 - mu1))
        pn <- 1/(1 + exp(-i00))
        a <- pn * (1 - pn)  #00
        b <- (exp(sigma * sigma) - 1) * exp(2 * mu + sigma^2) * pn * (1 - 
            pn)  #11
        c <- exp(mu + 0.5 * sigma^2) * pn * (1 - pn)  #01
        v0 <- a/(a * b - c^2)
        if (tside == 1) {
            s <- 1
            t <- 0
            alpha <- alpha
        }
        if (tside == 2) {
            s <- 1
            t <- 1
            alpha <- alpha/2
        }
        if (tside == 3) {
            s <- 0
            t <- 1
            alpha <- alpha
        }
        if (is.null(power)) 
            power <- eval(p.body)
        if (is.null(n)) 
            n <- uniroot(function(n) eval(p.body) - power, c(2 + 1e-10, 
                1e+07))$root
        
    }
    
    
    
    if (family == "normal") {
        g <- 0
        
        if (is.null(parameter)) {
            mu <- 0
            sigma <- 1
        } else {
            if (length(parameter) != 2) 
                stop("Both mean and standard deviation of the normal distribution have to be provided as a vector.")
            mu <- parameter[1]
            sigma <- parameter[2]
        }
        beta0 <- log(p0/(1 - p0))
        
        odds <- p1/(1 - p1)/(p0/(1 - p0))
        beta1 <- log(odds)  ##beta1 under H_A in the population
        
        d <- integrate(function(x) (1 - (1 + exp(-beta0 - beta1 * x))^(-1)) * 
            (1 + exp(-beta0 - beta1 * x))^(-1) * dnorm(x, mu, sigma), -Inf, 
            Inf, subdivisions = 100L)$value
        e <- integrate(function(x) x * (1 - (1 + exp(-beta0 - beta1 * x))^(-1)) * 
            (1 + exp(-beta0 - beta1 * x))^(-1) * dnorm(x, mu, sigma), -Inf, 
            Inf, subdivisions = 100L)$value
        f <- integrate(function(x) x^2 * (1 - (1 + exp(-beta0 - beta1 * 
            x))^(-1)) * (1 + exp(-beta0 - beta1 * x))^(-1) * dnorm(x, mu, 
            sigma), -Inf, Inf, subdivisions = 100L)$value
        v1 <- d/(d * f - e^2)
        
        mu1 <- integrate(function(x) 1/(1 + exp(-beta0 - beta1 * x)) * 
            dnorm(x, mu, sigma), -Inf, Inf, subdivisions = 100L)$value
        
        i00 <- log(mu1/(1 - mu1))
        pn <- 1/(1 + exp(-i00))
        a <- pn * (1 - pn)  #00
        b <- sigma^2 * pn * (1 - pn)  #11
        c <- mu * pn * (1 - pn)  #01
        v0 <- a/(a * b - c^2)
        if (tside == 1) {
            s <- 1
            t <- 0
            alpha <- alpha
        }
        if (tside == 2) {
            s <- 1
            t <- 1
            alpha <- alpha/2
        }
        if (tside == 3) {
            s <- 0
            t <- 1
            alpha <- alpha
        }
        if (is.null(power)) 
            power <- eval(p.body)
        if (is.null(n)) 
            n <- uniroot(function(n) eval(p.body) - power, c(2 + 1e-10, 
                1e+07))$root
        
    }
    
    
    if (family == "Poisson") {
        g <- 0
        
        if (is.null(parameter)) {
            lambda <- 1
        } else {
            lambda <- parameter
        }
        beta0 <- log(p0/(1 - p0))
        odds <- p1/(1 - p1)/(p0/(1 - p0))
        beta1 <- log(odds)  ##beta1 under H_A in the population
        
        d <- sum(sapply(0:1e+05, function(x) (1 - (1 + exp(-beta0 - beta1 * 
            x))^(-1)) * (1 + exp(-beta0 - beta1 * x))^(-1) * dpois(x, lambda)))
        e <- sum(sapply(0:1e+05, function(x) x * (1 - (1 + exp(-beta0 - 
            beta1 * x))^(-1)) * (1 + exp(-beta0 - beta1 * x))^(-1) * dpois(x, 
            lambda)))
        f <- sum(sapply(0:1e+05, function(x) x^2 * (1 - (1 + exp(-beta0 - 
            beta1 * x))^(-1)) * (1 + exp(-beta0 - beta1 * x))^(-1) * dpois(x, 
            lambda)))
        v1 <- d/(d * f - e^2)
        
        mu1 <- sum(sapply(0:1e+05, function(x) 1/(1 + exp(-beta0 - beta1 * 
            x)) * dpois(x, lambda)))
        
        i00 <- log(mu1/(1 - mu1))
        pn <- 1/(1 + exp(-i00))
        a <- pn * (1 - pn)  #00
        b <- lambda * pn * (1 - pn)  #11
        c <- lambda * pn * (1 - pn)  #01
        v0 <- a/(a * b - c^2)
        
        if (tside == 1) {
            s <- 1
            t <- 0
            alpha <- alpha
        }
        if (tside == 2) {
            s <- 1
            t <- 1
            alpha <- alpha/2
        }
        if (tside == 3) {
            s <- 0
            t <- 1
            alpha <- alpha
        }
        if (is.null(power)) 
            power <- eval(p.body)
        if (is.null(n)) 
            n <- uniroot(function(n) eval(p.body) - power, c(2 + 1e-10, 
                1e+07))$root
        
    }
    
    if (family == "uniform") {
        g <- 0
        if (is.null(parameter)) {
            L <- 0
            R <- 1
        } else {
            if (length(parameter) != 2) 
                stop("The lower and upper bounds have to be provided as a vector")
            L <- parameter[1]
            R <- parameter[2]
        }
        
        beta0 <- log(p0/(1 - p0))
        
        odds <- p1/(1 - p1)/(p0/(1 - p0))
        beta1 <- log(odds)  ##beta1 under H_A in the population
        
        
        d <- integrate(function(x) (1 - (1 + exp(-beta0 - beta1 * x))^(-1)) * 
            (1 + exp(-beta0 - beta1 * x))^(-1)/(R - L), L, R, subdivisions = 100L)$value
        e <- integrate(function(x) x * (1 - (1 + exp(-beta0 - beta1 * x))^(-1)) * 
            (1 + exp(-beta0 - beta1 * x))^(-1)/(R - L), L, R, subdivisions = 100L)$value
        f <- integrate(function(x) x^2 * (1 - (1 + exp(-beta0 - beta1 * 
            x))^(-1)) * (1 + exp(-beta0 - beta1 * x))^(-1)/(R - L), L, 
            R, subdivisions = 100L)$value
        v1 <- d/(d * f - e^2)
        
        mu1 <- integrate(function(x) 1/(1 + exp(-beta0 - beta1 * x))/(R - 
            L), L, R, subdivisions = 100L)$value
        
        i00 <- log(mu1/(1 - mu1))
        pn <- 1/(1 + exp(-i00))
        a <- pn * (1 - pn)  #00
        b <- (R - L)^2/12 * pn * (1 - pn)  #11
        c <- (R + L)/2 * pn * (1 - pn)  #01
        v0 <- a/(a * b - c^2)
        if (tside == 1) {
            s <- 1
            t <- 0
            alpha <- alpha
        }
        if (tside == 2) {
            s <- 1
            t <- 1
            alpha <- alpha/2
        }
        if (tside == 3) {
            s <- 0
            t <- 1
            alpha <- alpha
        }
        if (is.null(power)) 
            power <- eval(p.body)
        if (is.null(n)) 
            n <- uniroot(function(n) eval(p.body) - power, c(2 + 1e-10, 
                1e+07))$root
    }
    
    METHOD <- "Power for logistic regression"
    URL <- "http://psychstat.org/logistic"
    structure(list(p0 = p0, p1 = p1, beta0 = beta0, beta1 = beta1, n = n, 
        alpha = sig.level, power = power, alternative = alternative, family = family, 
        parameter = parameter, method = METHOD, url = URL), class = "webpower")
}

wp.mediation <- function(n = NULL, power = NULL, a = 0.5, b = 0.5, varx = 1, 
    vary = 1, varm = 1, alpha = 0.05, interval = NULL) {
    if (sum(sapply(list(n, a, b, varx, varm, vary, power, alpha), is.null)) != 
        1) 
        stop("exactly one of n, a, b, varx, varm, vary, power, and alpha must be NULL")
    if (!is.null(alpha) && !is.numeric(alpha) || any(0 > alpha | alpha > 
        1)) 
        stop(sQuote("alpha"), " must be numeric in [0, 1]")
    if (!is.null(power) && !is.numeric(power) || any(0 > power | power > 
        1)) 
        stop(sQuote("power"), " must be numeric in [0, 1]")
    
    ## power formulae
    p.body <- quote({
        numer <- sqrt(n) * a * b
        denom <- sqrt(a^2 * vary/(varm - a^2 * varx) + b^2 * (varm - a^2 * 
            varx)/varx)
        delta <- numer/denom
        alpha2 <- alpha/2
        za2 <- qnorm(1 - alpha2)
        1 - pnorm(za2 - delta) + pnorm(-za2 - delta)
    })
    
    if (is.null(power)) 
        power <- eval(p.body) else if (is.null(n)) 
        n <- uniroot(function(n) eval(p.body) - power, ifelse(rep(is.null(interval), 2), c(1e-07, 1e+07), interval))$root else if (is.null(a)) {
        astart <- varm/varx
        alow <- -sqrt(astart) + 1e-06
        aup <- sqrt(astart) - 1e-06
        a <- nuniroot(function(a) eval(p.body) - power, c(alow, aup))$root
    } else if (is.null(b)) 
        b <- nuniroot(function(b) eval(p.body) - power, ifelse(rep(is.null(interval), 2), c(-10, 10), interval))$root else if (is.null(varx)) 
        varx <- nuniroot(function(varx) eval(p.body) - power, ifelse(rep(is.null(interval), 2), c(1e-10, 1e+07), interval))$root else if (is.null(vary)) 
        vary <- nuniroot(function(vary) eval(p.body) - power, ifelse(rep(is.null(interval), 2), c(1e-10, 1e+07), interval))$root else if (is.null(varm)) 
        varm <- nuniroot(function(varm) eval(p.body) - power, ifelse(rep(is.null(interval), 2), c(1e-10, 1e+07), interval))$root else if (is.null(alpha)) 
        alpha <- nuniroot(function(alpha) eval(p.body) - power, ifelse(rep(is.null(interval), 2), c(1e-10, 1 - 1e-10), interval))$root else stop("internal error")
    
    METHOD <- "Power for simple mediation"
    URL <- "http://psychstat.org/mediation"
    
    structure(list(n = n, power = power, a = a, b = b, varx = varx, varm = varm, 
        vary = vary, alpha = alpha, method = METHOD, url = URL), class = "webpower")
}

estCRT2arm <- function(file) {
    dat <- read.table(file, header = TRUE)
    X <- dat[, 4] - 0.5
    newdata <- data.frame(dat, X)
    model <- summary(lmer(score ~ X + (1 | cluster), data = newdata))
    parest <- model$coefficients
    J <- model$ngrps
    p.value <- (1 - pt(abs(parest[, 3]), J - 2)) * 2
    test <- cbind(parest, p.value)
    cat("Fixed effects (X: Treatment main effect):\n")
    print(test)
}

estCRT3arm <- function(file) {
    dat <- read.table(file, header = T)
    X1 <- X2 <- rep(0, dim(dat)[1])
    X1[dat[, 4] == 1] <- X1[dat[, 4] == 2] <- 1/3
    X1[dat[, 4] == 0] <- -2/3
    X2[dat[, 4] == 1] <- 1/2
    X2[dat[, 4] == 2] <- -1/2
    newdata <- data.frame(dat, X1, X2)
    model <- summary(lmer(score ~ X1 + X2 + (1 | cluster), data = newdata))
    parest <- model$coefficients
    rownames(parest)[2] <- "X"
    rownames(parest)[3] <- "X1-X2"
    J <- model$ngrps
    p.value <- (1 - pt(abs(parest[, 3]), J - 3)) * 2
    test <- cbind(parest, p.value)
    # Omnibus test
    model1 <- lmer(score ~ X1 + X2 + (1 | cluster), data = newdata, REML = FALSE)
    model2 <- lmer(score ~ (1 | cluster), data = newdata, REML = FALSE)
    omnibus <- anova(model1, model2)$"Pr(>Chisq)"[2]
    cat("Fixed effects (X: Treatment main effect; X1-X2: Comparing two treatments):\n")
    print(test)
    cat("Omnibust test: p-value = ")
    cat(omnibus)
}



estMRT2arm <- function(file) {
    dat <- read.table(file, header = TRUE)
    X <- dat[, 4] - 0.5
    newdata <- data.frame(dat, X)
    model <- summary(lmer(score ~ X + (X | cluster), data = newdata))
    parest <- model$coefficients
    J <- model$ngrps
    p.value <- (1 - pt(abs(parest[, 3]), J - 1)) * 2
    test <- cbind(parest, p.value)
    cat("Fixed effects (X: Treatment main effect):\n")
    print(test)
}


estMRT3arm <- function(file) {
    dat <- read.table(file, header = T)
    X1 <- X2 <- rep(0, dim(dat)[1])
    X1[dat[, 4] == 1] <- X1[dat[, 4] == 2] <- 1/3
    X1[dat[, 4] == 0] <- -2/3
    X2[dat[, 4] == 1] <- 1/2
    X2[dat[, 4] == 2] <- -1/2
    newdata <- data.frame(dat, X1, X2)
    model <- summary(lmer(score ~ X1 + X2 + (X1 + X2 | cluster), data = newdata))
    parest <- model$coefficients
    rownames(parest)[2] <- "X"
    rownames(parest)[3] <- "X1-X2"
    J <- model$ngrps
    p.value <- (1 - pt(abs(parest[, 3]), J - 1)) * 2
    test <- cbind(parest, p.value)
    # Omnibus test
    model1 <- lmer(score ~ X1 + X2 + (X1 + X2 | cluster), data = newdata, 
        REML = FALSE)
    model2 <- lmer(score ~ (X1 + X2 | cluster), data = newdata, REML = FALSE)
    omnibus <- anova(model1, model2)$"Pr(>Chisq)"[2]
    cat("Fixed effects (X: Treatment main effect; X1-X2: Comparing two treatments):\n")
    print(test)
    cat("Omnibust test: p-value = ")
    cat(omnibus)
}

wp.correlation <- function(n = NULL, r = NULL, power = NULL, p = 0, rho0 = 0, 
    alpha = 0.05, alternative = c("two.sided", "less", "greater")) {
    if (sum(sapply(list(n, r, power, alpha), is.null)) != 1) 
        stop("exactly one of n, r, power, and alpha must be NULL")
    if (!is.null(alpha) && !is.numeric(alpha) || any(0 > alpha | alpha > 
        1)) 
        stop(sQuote("alpha"), " must be numeric in [0, 1]")
    if (!is.null(power) && !is.numeric(power) || any(0 > power | power > 
        1)) 
        stop(sQuote("power"), " must be numeric in [0, 1]")
    if (!is.null(n) && min(n) < 4) 
        stop("number of observations must be at least 4")
    alternative <- match.arg(alternative)
    tside <- switch(alternative, less = 1, two.sided = 2, greater = 3)
    if (tside == 2 && !is.null(r)) 
        r <- abs(r)
    if (tside == 3) {
        p.body <- quote({
            delta <- sqrt(n - 3 - p) * (log((1 + r)/(1 - r))/2 + r/(n - 
                1 - p)/2 * (1 + (5 + r^2)/(n - 1 - p)/4 + (11 + 2 * r^2 + 
                3 * r^4)/(n - 1 - p)^2/8) - log((1 + rho0)/(1 - rho0))/2 - 
                rho0/(n - 1 - p)/2)
            v <- (n - 3 - p)/(n - 1 - p) * (1 + (4 - r^2)/(n - 1 - p)/2 + 
                (22 - 6 * r^2 - 3 * r^4)/(n - 1 - p)^2/6)
            zalpha <- qnorm(1 - alpha)
            pnorm((delta - zalpha)/sqrt(v))
        })
    }
    if (tside == 1) {
        p.body <- quote({
            delta <- sqrt(n - 3 - p) * (log((1 + r)/(1 - r))/2 + r/(n - 
                1 - p)/2 * (1 + (5 + r^2)/(n - 1 - p)/4 + (11 + 2 * r^2 + 
                3 * r^4)/(n - 1 - p)^2/8) - log((1 + rho0)/(1 - rho0))/2 - 
                rho0/(n - 1 - p)/2)
            v <- (n - 3 - p)/(n - 1 - p) * (1 + (4 - r^2)/(n - 1 - p)/2 + 
                (22 - 6 * r^2 - 3 * r^4)/(n - 1 - p)^2/6)
            zalpha <- qnorm(1 - alpha)
            pnorm((-delta - zalpha)/sqrt(v))
        })
    }
    if (tside == 2) {
        p.body <- quote({
            delta <- sqrt(n - 3 - p) * (log((1 + r)/(1 - r))/2 + r/(n - 
                1 - p)/2 * (1 + (5 + r^2)/(n - 1 - p)/4 + (11 + 2 * r^2 + 
                3 * r^4)/(n - 1 - p)^2/8) - log((1 + rho0)/(1 - rho0))/2 - 
                rho0/(n - 1 - p)/2)
            v <- (n - 3 - p)/(n - 1 - p) * (1 + (4 - r^2)/(n - 1 - p)/2 + 
                (22 - 6 * r^2 - 3 * r^4)/(n - 1 - p)^2/6)
            zalpha <- qnorm(1 - alpha/2)
            pnorm((delta - zalpha)/sqrt(v)) + pnorm((-delta - zalpha)/sqrt(v))
        })
    }
    if (is.null(power)) 
        power <- eval(p.body) else if (is.null(n)) 
        n <- uniroot(function(n) eval(p.body) - power, c(4 + p + 1e-10, 
            1e+07))$root else if (is.null(r)) {
        if (tside == 2) {
            r <- uniroot(function(r) eval(p.body) - power, c(1e-10, 1 - 
                1e-10))$root
        } else {
            r <- uniroot(function(r) eval(p.body) - power, c(-1 + 1e-10, 
                1 - 1e-10))$root
        }
    } else if (is.null(alpha)) 
        alpha <- uniroot(function(alpha) eval(p.body) - power, c(1e-10, 
            1 - 1e-10))$root else stop("internal error")
    METHOD <- "Power for correlation"
    URL <- "http://psychstat.org/correlation"
    structure(list(n = n, r = r, alpha = alpha, power = power, alternative = alternative, 
        method = METHOD, url = URL), class = "webpower")
}

wp.poisson <- function(n = NULL, exp0 = NULL, exp1 = NULL, alpha = 0.05, power = NULL,
                        alternative = c("two.sided", "less", "greater"),
                        family = c("Bernoulli", "exponential", "lognormal",
                                   "normal", "Poisson", "uniform"), 
                        parameter = NULL, subdivisions=200L,
                       i.method=c("numerical", "MC"), mc.iter=20000)
{
  sig.level <- alpha
  if(sum(sapply(list(n, power), is.null)) != 1)
     stop("exactly one of n, power, andalpha must be NULL")
  if(is.null(exp0) || !is.null(exp0) && !is.numeric(exp0) ||
      exp0 <= 0)
    stop(sQuote("exp0"), " must be numeric in (0,Infinity)")
  if(is.null(exp1) || !is.null(exp1) && !is.numeric(exp1) ||
      exp1 <= 0)
    stop(sQuote("exp1"), " must be numeric in (0,Infinity)")
  if(!is.null(n) && min(n) < 5)
    stop("number of observations must be at least 5")
  if(is.null(alpha) || !is.null(alpha) & !is.numeric(alpha) ||
      any(alpha < 0 | alpha > 1))
    stop(sQuote("alpha"), " must be numeric in [0, 1]")
  if(!is.null(power) & !is.numeric(power) || any(0 > power |
                                                  power > 1))
    stop(sQuote("power"), " must be numeric in [0, 1]")
  alternative <- match.arg(alternative)
  tside <- switch(alternative, less = 1, two.sided = 2, greater = 3)
  p.body <- quote({
    s * pnorm(-qnorm(1 - alpha) - sqrt(n)/sqrt(v1) * beta1) +
      t * pnorm(-qnorm(1 - alpha) + sqrt(n)/sqrt(v1) *
                  beta1)
  })
  if (family == "Bernoulli") {
    if (is.null(parameter)) {
      B <- 0.5
    }
    else {
      B <- parameter
    }
    beta1 <- log(exp1)
    beta0 <- log(exp0)
    a <- (1 - B) * exp(beta0) + B * exp(beta0 + beta1)
    b <- B * exp(beta0 + beta1)
    c <- B * exp(beta0 + beta1)
    v1 <- a/(a * c - b^2)
    if (tside == 1) {
      s <- 1
      t <- 0
      alpha <- alpha
    }
    if (tside == 2) {
      s <- 1
      t <- 1
      alpha <- alpha/2
    }
    if (tside == 3) {
      s <- 0
      t <- 1
      alpha <- alpha
    }
    if (is.null(power)) {
      power <- eval(p.body)
    }
    if (is.null(n)) {
      n <- uniroot(function(n) eval(p.body) - power, c(2 +
                                                         1e-10, 1e+07))$root
    }
  }
  if (family == "exponential") {
    if (is.null(parameter)) {
      lambda <- 1
    }
    else {
      lambda <- parameter
    }
    beta1 <- log(exp1); beta0 <- log(exp0)
    
    if (i.method[1] == "MC"){
      x.s=rexp(mc.iter,lambda)
      d=mean(exp0*exp(log(exp1)*x.s))
      e=mean(x.s*exp0*exp(log(exp1)*x.s))
      f=mean(x.s^2*exp0*exp(log(exp1)*x.s)) 
    }else{
      d.integ<-function(x,beta0,beta1,lamba){#beta1!=lambda
        
        lambda*exp(beta0)/(beta1-lambda)*exp((beta1-lambda)*x)
        
      }
      
      e.integ<-function(x,beta0,beta1,lambda){
        
        lambda*exp(beta0)/(beta1-lambda)*x*exp((beta1-lambda)*x) -1/(beta1-lambda)*d.integ(x,beta0,beta1,lambda)
        
      }
      
      f.integ<-function(x,beta0,beta1,lambda){
        
        x^2*d.integ(x,beta0,beta1,lambda)-2/(beta1-lambda)*e.integ(x,beta0,beta1,lambda)
        
      }
      
      
      
      if(beta1<lambda){
        
        d=d.integ(100/(lambda-beta1),beta0,beta1,lambda)-d.integ(0,beta0,beta1,lambda)
        
        e=e.integ(100/(lambda-beta1),beta0,beta1,lambda)-e.integ(0,beta0,beta1,lambda)
        
        f=f.integ(100/(lambda-beta1), beta0,beta1,lambda)-f.integ(0,beta0,beta1,lambda)
        
        v1=d/(d*f-e^2)}
      
      if(beta1>lambda){
        
        d=d.integ(100/abs(lambda-beta1),beta0,beta1,lambda)-d.integ(0,beta0,beta1,lambda)
        
        e=e.integ(100/abs(lambda-beta1),beta0,beta1,lambda)-e.integ(0,beta0,beta1,lambda)
        
        f=f.integ(100/abs(lambda-beta1), beta0,beta1, lambda)-f.integ(0,beta0,beta1,lambda)
        
        v1=d/(d*f-e^2)
        
      }else{
        
        v1=0.00001#the actual value is 0
        
      }
    }
   
    v1 <- d/(d * f - e^2)
    if (tside == 1) {
      s <- 1
      t <- 0
      alpha <- alpha
    }
    if (tside == 2) {
      s <- 1
      t <- 1
      alpha <- alpha/2
    }
    if (tside == 3) {
      s <- 0
      t <- 1
      alpha <- alpha
    }
    if (is.null(power)) {
      power <- eval(p.body)
    }
    if (is.null(n)) {
      n <- uniroot(function(n) eval(p.body) - power, c(2 +
                                                         1e-10, 1e+07))$root
    }
  }
  if (family == "lognormal") {
    if (is.null(parameter)) {
      mu <- 0
      sigma <- 1
    }
    else {
      if (length(parameter) != 2)
        stop("Both mean and standard deviation of the log-normal distribution have to be provided as a vector.")
      mu <- parameter[1]
      sigma <- parameter[2]
    }
    beta1 <- log(exp1)
    beta0 <- log(exp0)
    
    if (i.method[1] == "MC"){
      x.s=rlnorm(mc.iter,mu,sigma)
      d=mean(exp0*exp(log(exp1)*x.s))
      e=mean(x.s*exp0*exp(log(exp1)*x.s))
      f=mean(x.s^2*exp0*exp(log(exp1)*x.s))
    }else{
      d <- integrate(function(x) exp(beta0 + beta1 * x) * dlnorm(x, mu, sigma), 0, Inf, subdivisions = subdivisions)$value
      e <- integrate(function(x) x * exp(beta0 + beta1 * x) * dlnorm(x, mu, sigma), 0, Inf, subdivisions = subdivisions)$value
      f <- integrate(function(x) x^2 * exp(beta0 + beta1 * x) * dlnorm(x, mu, sigma), 0, Inf, subdivisions = subdivisions)$value
    }
   
   
     v1 <- d/(d * f - e^2)
    if (tside == 1) {
      s <- 1
      t <- 0
      alpha <- alpha
    }
    if (tside == 2) {
      s <- 1
      t <- 1
      alpha <- alpha/2
    }
    if (tside == 3) {
      s <- 0
      t <- 1
      alpha <- alpha
    }
    if (is.null(power)) {
      power <- eval(p.body)
    }
    if (is.null(n)) {
      n <- uniroot(function(n) eval(p.body) - power, c(2 +
                                                         1e-10, 1e+07))$root
    }
  }
  if (family == "normal") {
    if (is.null(parameter)) {
      mu <- 0
      sigma <- 1
    }
    else {
      if (length(parameter) != 2)
        stop("Both mean and standard deviation of the normal distribution have to be provided as a vector.")
      mu <- parameter[1]
      sigma <- parameter[2]
    }
    beta1 <- log(exp1)
    beta0 <- log(exp0)
    
    
      d=exp0*exp(beta1^2*sigma^2/2+beta1*mu)
      
      e=exp0*exp(beta1^2*sigma^2/2+beta1*mu)*(mu+beta1*sigma^2)
      
      f=exp0*exp(beta1^2*sigma^2/2+beta1*mu)*((mu+beta1*sigma^2)^2+sigma^2)
    
    
    v1 <- d/(d * f - e^2)
    if (tside == 1) {
      s <- 1
      t <- 0
      alpha <- alpha
    }
    if (tside == 2) {
      s <- 1
      t <- 1
      alpha <- alpha/2
    }
    if (tside == 3) {
      s <- 0
      t <- 1
      alpha <- alpha
    }
    if (is.null(power)) {
      power <- eval(p.body)
    }
    if (is.null(n)) {
      n <- uniroot(function(n) eval(p.body) - power, c(2 +
                                                         1e-10, 1e+07))$root
    }
  }
  if (family == "Poisson") {
    if (is.null(parameter)) {
      lambda <- 1
    }
    else {
      lambda <- parameter
    }
    beta1 <- log(exp1)
    beta0 <- log(exp0)
    d <- sum(sapply(0:1e+05, function(x) exp(beta0 + beta1 *
                                               x) * dpois(x, lambda)))
    e <- sum(sapply(0:1e+05, function(x) x * exp(beta0 +
                                                   beta1 * x) * dpois(x, lambda)))
    f <- sum(sapply(0:1e+05, function(x) x^2 * exp(beta0 +
                                                     beta1 * x) * dpois(x, lambda)))
    v1 <- d/(d * f - e^2)
    if (tside == 1) {
      s <- 1
      t <- 0
      alpha <- alpha
    }
    if (tside == 2) {
      s <- 1
      t <- 1
      alpha <- alpha/2
    }
    if (tside == 3) {
      s <- 0
      t <- 1
      alpha <- alpha
    }
    if (is.null(power)) {
      power <- eval(p.body)
    }
    if (is.null(n)) {
      n <- uniroot(function(n) eval(p.body) - power, c(2 +
                                                         1e-10, 1e+07))$root
    }
  }
  if (family == "uniform") {
    if (is.null(parameter)) {
      L <- 0
      R <- 1
    }
    else {
      if (length(parameter) != 2)
        stop("The lower and upper bounds have to be provided as a vector")
      L <- parameter[1]
      R <- parameter[2]
    }
    beta1 <- log(exp1)
    beta0 <- log(exp0)
    d <- integrate(function(x) exp(beta0 + beta1 * x)/(R -
                                                         L), L, R, subdivisions = subdivisions)$value
    e <- integrate(function(x) x * exp(beta0 + beta1 * x)/(R -
                                                             L), L, R, subdivisions = subdivisions)$value
    f <- integrate(function(x) x^2 * exp(beta0 + beta1 *
                                           x)/(R - L), L, R, subdivisions = subdivisions)$value
    v1 <- d/(d * f - e^2)
    if (tside == 1) {
      s <- 1
      t <- 0
      alpha <- alpha
    }
    if (tside == 2) {
      s <- 1
      t <- 1
      alpha <- alpha/2
    }
    if (tside == 3) {
      s <- 0
      t <- 1
      alpha <- alpha
    }
    if (is.null(power)) {
      power <- eval(p.body)
    }
    if (is.null(n)) {
      n <- uniroot(function(n) eval(p.body) - power, c(2 +
                                                         1e-10, 1e+07))$root
    }
  }
  METHOD <- "Power for Poisson regression"
  URL <- "http://psychstat.org/poisson"
  structure(list(n = n, power = power, alpha = sig.level, exp0 = exp0,
                exp1 = exp1, beta0 = beta0, beta1 = beta1, alternative = alternative,
                family = family, method = METHOD,
                url = URL), class = "webpower")
}


wp.sem.chisq <- 
function (n = NULL, df = NULL, effect = NULL, power = NULL, alpha = 0.05) 
{
    if (sum(sapply(list(n, df, power, alpha, effect), is.null)) != 1) 
        stop("exactly one of sample size, degrees of freedom, effect size, alpha and power must be NULL")
    if (!is.null(alpha) && !is.numeric(alpha) || any(0 > 
        alpha | alpha > 1)) 
        stop(sQuote("alpha"), " must be numeric in [0, 1]")
    if (!is.null(power) && !is.numeric(power) || any(0 > power | 
        power > 1)) 
        stop(sQuote("power"), " must be numeric in [0, 1]")

    p.body <- quote({
    		ncp<-(n-1)*effect
    		#ncp<-n*effect
        	calpha<-qchisq(1-alpha, df)
        	1-pchisq(calpha, df, ncp)
        	})
            
    if (is.null(power)) 
        power <- eval(p.body)
    else if (is.null(n)) 
        n <- uniroot(function(n) eval(p.body) - power, c(10 + 
            1e-10, 1e+08))$root
    else if (is.null(df))
        df <- uniroot(function(df) eval(p.body) - power, c(1, 
                10000))$root
    else if (is.null(alpha)) 
        alpha <- uniroot(function(alpha) eval(p.body) - 
            power, c(1e-10, 1 - 1e-10))$root
    else if (is.null(effect))
    	effect <- uniroot(function(effect) eval(p.body) - power, c(0, 1))$root
    else stop("internal error")
	METHOD <- "Power for SEM (Satorra & Saris, 1985)"
	URL <- "http://psychstat.org/semchisq"
    structure(list(n = n, df = df, effect = effect, power = power, alpha = alpha, url=URL, method=METHOD), 
        class = "webpower")
}

wp.sem.rmsea <- 
function (n = NULL, df = NULL, rmsea0 = NULL, rmsea1 = NULL,  power = NULL, alpha = 0.05, type=c('close','notclose')) {
    if (sum(sapply(list(n, df, power, alpha, rmsea0, rmsea1), is.null)) != 1) 
        stop("exactly one of sample size, degrees of freedom, rmsea0, rmsea1, alpha and power must be NULL")
    if (!is.null(alpha) && !is.numeric(alpha) || any(0 > 
        alpha | alpha > 1)) 
        stop(sQuote("alpha"), " must be numeric in [0, 1]")
    if (!is.null(power) && !is.numeric(power) || any(0 > power | 
        power > 1)) 
        stop(sQuote("power"), " must be numeric in [0, 1]")
    
    type <- match.arg(type)    
    htype <- switch(type, close = 1, notclose = 2)

    if (htype==1){
    	p.body <- quote({
    		ncp0<-(n-1)*df*rmsea0^2
    		ncp1<-(n-1)*df*rmsea1^2
        	calpha<-qchisq(1-alpha, df, ncp0)
        	1-pchisq(calpha, df, ncp1)
        	})
    }else{
    	p.body <- quote({
    		ncp0<-(n-1)*df*rmsea0^2
    		ncp1<-(n-1)*df*rmsea1^2
        	calpha<-qchisq(alpha, df, ncp0)
        	pchisq(calpha, df, ncp1)
        	})
    }
        
    
    if (is.null(power)) 
        power <- eval(p.body)
    else if (is.null(n)) 
        n <- uniroot(function(n) eval(p.body) - power, c(2 + 
            1e-10, 1e+05))$root
    else if (is.null(df))
        df <- uniroot(function(df) eval(p.body) - power, c(1, 
                10000))$root
    else if (is.null(alpha)) 
        alpha <- uniroot(function(alpha) eval(p.body) - 
            power, c(1e-10, 1 - 1e-10))$root
    else if (is.null(rmsea0))
    	rmsea0 <- uniroot(function(rmsea0) eval(p.body) - 
            power, c(0, 1))$root
    else if (is.null(rmsea1))
    	rmsea1 <- uniroot(function(rmsea1) eval(p.body) - 
            power, c(0, 1))$root
    else stop("internal error")
	METHOD <- "Power for SEM based on RMSEA"
	URL <- "http://psychstat.org/rmsea"
    structure(list(n = n, df = df, rmsea0 = rmsea0, rmsea1 = rmsea1, power = power, alpha = alpha, method=METHOD, url=URL), 
        class = "webpower")
}

wp.kanova<-function (n = NULL, ndf = NULL, f = NULL, ng=NULL, alpha = 0.05, power = NULL) 
{
    if (sum(sapply(list(n, ndf, f, ng, power, alpha), is.null)) != 
        1) 
        stop("exactly one of n, ndf, f, ng, power, alpha must be NULL")
    if (!is.null(f) && min(f) < 0) 
        stop("f must be positive")
    if (!is.null(n) && min(n) < 1) 
        stop("degree of freedom u for numerator must be at least 1")
    if (!is.null(ndf) && min(ndf) < 1) 
        stop("degree of freedom v for denominator must be at least 1")
	if (!is.null(ng) && min(ng) < 1) 
        stop("number of groups must be at least 1")
    if (!is.null(alpha) && !is.numeric(alpha) || any(0 > 
        alpha | alpha > 1)) 
        stop(sQuote("alpha"), " must be numeric in [0, 1]")
    if (!is.null(power) && !is.numeric(power) || any(0 > power | 
        power > 1)) 
        stop(sQuote("power"), " must be numeric in [0, 1]")
    p.body <- quote({
        lambda <- f^2 * n
		ddf <- n - ng
        pf(qf(alpha, ndf, ddf, lower = FALSE), ndf, ddf, lambda, 
            lower = FALSE)
    })
    if (is.null(power)) 
        power <- eval(p.body)
    else if (is.null(n)){ 
        n <- uniroot(function(n) eval(p.body) - power, c(1 + 
            ng, 1e+07))$root
		ddf <- n - ng
		}
    else if (is.null(ndf)){ 
        ndf <- uniroot(function(ndf) eval(p.body) - power, c(1 + 
            1e-10, 1e+05))$root
			ddf <- n -ng
		}
	else if (is.null(ng)) {
        ng <- uniroot(function(ng) eval(p.body) - power, c(1e-07, 
            1e+07))$root
			ddf <- n -ng
		}
    else if (is.null(f)){ 
        f <- uniroot(function(f) eval(p.body) - power, c(1e-07, 
            1e+07))$root
			ddf <- n -ng
		}
    else if (is.null(alpha)){ 
        alpha <- uniroot(function(alpha) eval(p.body) - 
            power, c(1e-10, 1 - 1e-10))$root
			ddf <- n -ng
		}
    else stop("internal error")
    METHOD <- "Multiple way ANOVA analysis"
	URL <- "http://psychstat.org/kanova"
	NOTE <- "Sample size is the total sample size"
    structure(list(n = n, ndf = ndf, ddf=ddf, f = f,  ng=ng, alpha = alpha, 
        power = power, method = METHOD, url=URL, note=NOTE), class = "webpower")
}

## ANOVA with binary data
wp.anova.binary<-function(k = NULL, n = NULL, V = NULL, alpha = 0.05, power = NULL){
  if (sum(sapply(list(k, n, V, power, alpha), is.null)) != 1) 
    stop("exactly one of k, n, V, power, and alpha must be NULL")
  if (!is.null(V) && min(V) < 0) 
    stop("V must be positive")
  if (!is.null(k) && min(k) < 2) 
    stop("number of groups must be at least 2")
  if (!is.null(n) && min(n) < 2) 
    stop("number of observations in each group must be at least 2")
  if (!is.null(alpha) && !is.numeric(alpha) || any(0 > alpha | alpha > 1)) 
    stop(sQuote("alpha"), " must be numeric in [0, 1]")
  if (!is.null(power) && !is.numeric(power) || any(0 > power | power > 1)) 
    stop(sQuote("power"), " must be numeric in [0, 1]")
  p.body <- quote({
    chi <- (V)^2*n*(k-1)
    df <- k - 1
    crit.val <-  qchisq(p=1-alpha,df=df,ncp=0)
    1-pchisq(crit.val,df=df,ncp=chi)
  })
  
  if (is.null(power)) 
    power <- eval(p.body)
  else if (is.null(k)) 
    k <- uniroot(function(k) eval(p.body) - power, c(2 + 1e-10, 100))$root
  else if (is.null(n)) 
    n <- uniroot(function(n) eval(p.body) - power, c(2 + k + 1e-10, 1e+05))$root
  else if (is.null(V)) 
    V <- uniroot(function(V) eval(p.body) - power, c(1e-07, 1e+07))$root
  else if (is.null(alpha)) 
    alpha <- uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
  else stop("internal error")
  NOTE <- "n is the total sample size"
  
  METHOD <- "One-way Analogous ANOVA with Binary Data"
  URL <- "http://psychstat.org/anovabinary"
  structure(list(k = k, n = n, V = V, alpha = alpha, power = power, 
                 note = NOTE, method = METHOD, url=URL), class = "webpower")
}

wp.anova.count<-function(k = NULL, n = NULL, V = NULL, alpha = 0.05, power = NULL){
  if (sum(sapply(list(k, n, V, power, alpha), is.null)) != 1) 
    stop("exactly one of k, n, V, power, and alpha must be NULL")
  if (!is.null(V) && min(V) < 0) 
    stop("V must be positive")
  if (!is.null(k) && min(k) < 2) 
    stop("number of groups must be at least 2")
  if (!is.null(n) && min(n) < 2) 
    stop("number of observations in each group must be at least 2")
  if (!is.null(alpha) && !is.numeric(alpha) || any(0 > alpha | alpha > 1)) 
    stop(sQuote("alpha"), " must be numeric in [0, 1]")
  if (!is.null(power) && !is.numeric(power) || any(0 > power | power > 1)) 
    stop(sQuote("power"), " must be numeric in [0, 1]")
  p.body <- quote({
    chi <- (V)^2*n*(k-1)
    df <- k - 1
    crit.val <-  qchisq(p=1-alpha,df=df,ncp=0)
    1-pchisq(crit.val,df=df,ncp=chi)
  })
  
  if (is.null(power)) 
    power <- eval(p.body)
  else if (is.null(k)) 
    k <- uniroot(function(k) eval(p.body) - power, c(2 + 1e-10, 100))$root
  else if (is.null(n)) 
    n <- uniroot(function(n) eval(p.body) - power, c(2 + k + 1e-10, 1e+05))$root
  else if (is.null(V)) 
    V <- uniroot(function(V) eval(p.body) - power, c(1e-07, 1e+07))$root
  else if (is.null(alpha)) 
    alpha <- uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
  else stop("internal error")
  NOTE <- "n is the total sample size"
  
  METHOD <- "One-way Analogous ANOVA with Count Data"
  URL <- "http://psychstat.org/anovacount"
  structure(list(k = k, n = n, V = V, alpha = alpha, power = power, note = NOTE, method = METHOD, url=URL), class = "webpower")
}

### Univariate Latent Change Score Models
wp.lcsm<-function(N=100, T=5, R=1000,
	betay=0, my0=0, mys=0, varey=1, vary0=1, varys=1, vary0ys=0, alpha=0.05, ...){
	#if (sum(N < 2*T)>0) stop("The sample size has to be at least 2 times of the number of occasions")
	
	pop.model <- function(T){
	## latent y
	## Intercept
		model<-"y0 =~ 1*y1\n"
		
		## path from y(t-1) to y(t) with path 1
		for (i in 2:T){
			model<-paste(model, "y",i,"~1*y",(i-1),"\n", sep="")
		}
		
		## loading from dy(t) to y(t) with path 1
		for (i in 2:T){
			model<-paste(model, "dy",i,"=~1*y",i,"\n", sep="")
		}
		
		## path from y(t) to dy(t+1) with path betay
		for (i in 2:T){
			model<-paste(model, "dy",i,"~", betay, "*y", (i-1), "\n", sep="")
		}
		
		## latent slope ys factor model
		for (i in 2:T){
			model<-paste(model, "ys=~1*dy", i, "\n", sep="")
		}
		
		## variance for dy constraints to 0
		for (i in 2:T){
			model<-paste(model, "dy",i,"~~0*dy",i,"\n", sep="")
		}
		
		## variance for y constraints to 0
		for (i in 1:T){
			model<-paste(model, "y",i,"~~0*y",i,"\n", sep="")
		}
		
		## variance and covariance for intercept and slope
		model<-paste(model, "ys~~", vary0ys, "*y0\n", sep="")
		model<-paste(model, "y0~~", vary0, "*y0\n", sep="")
		model<-paste(model, "ys~~", varys, "*ys\n", sep="")
		
		model<-paste(model, "ys~", mys, "*1\n", sep="")	
		model<-paste(model, "y0~", my0, "*1\n", sep="")

		## constrain means of y and dy to be zero
		for (i in 1:T){
			model<-paste(model, "y",i,"~0*1\n", sep="")
		}
		for (i in 2:T){
			model<-paste(model, "dy",i,"~0*1\n", sep="")
		}
		
		
		## for observed data part
		## y(t) to Y(t)
		for (i in 1:T){
			model<-paste(model, "y",i,"=~1*", "Y",i, "\n", sep="")				
		}
		
		## set means of Y to be zero
		for (i in 1:T){
			model<-paste(model, "Y",i, "~0*1\n", sep="")		
		}
		
		## set the variance for Y
		for (i in 1:T){
			model<-paste(model, "Y",i, "~~", varey, "*", "Y",i, "\n", sep="")		
		}
		model
	}
	
	fit.model <- function(T){
	## latent y
	## Intercept
		model<-"y0 =~ 1*y1\n"
		
		## path from y(t-1) to y(t) with path 1
		for (i in 2:T){
			model<-paste(model, "y",i,"~1*y",(i-1),"\n", sep="")
		}
		
		## loading from dy(t) to y(t) with path 1
		for (i in 2:T){
			model<-paste(model, "dy",i,"=~1*y",i,"\n", sep="")
		}
		
		## path from y(t) to dy(t+1) with path betay
		for (i in 2:T){
			model<-paste(model, "dy",i,"~start(", betay, ")*y", (i-1)," + betay*y", (i-1), "\n", sep="")
		}
		
		## latent slope ys factor model
		for (i in 2:T){
			model<-paste(model, "ys=~1*dy", i, "\n", sep="")
		}
		
		## variance for dy constraints to 0
		for (i in 2:T){
			model<-paste(model, "dy",i,"~~0*dy",i,"\n", sep="")
		}
		
		## variance for y constraints to 0
		for (i in 1:T){
			model<-paste(model, "y",i,"~~0*y",i,"\n", sep="")
		}
		
		## variance and covariance for intercept and slope
		model<-paste(model, "ys~~start(", vary0ys, ")*y0 + vary0ys*y0\n", sep="")
		model<-paste(model, "y0~~start(", vary0, ")*y0 + vary0*y0\n", sep="")
		model<-paste(model, "ys~~start(", varys, ")*ys + varys*ys\n", sep="")
		
		model<-paste(model, "ys~start(", mys, ")*1 + label('mys')*1\n", sep="")	
		model<-paste(model, "y0~start(", my0, ")*1 + label('my0')*1\n", sep="")

		## constrain means of y and dy to be zero
		for (i in 1:T){
			model<-paste(model, "y",i,"~0*1\n", sep="")
		}
		for (i in 2:T){
			model<-paste(model, "dy",i,"~0*1\n", sep="")
		}
		
		
		## for observed data part
		## y(t) to Y(t)
		for (i in 1:T){
			model<-paste(model, "y",i,"=~1*", "Y",i, "\n", sep="")				
		}
		
		## set means of Y to be zero
		for (i in 1:T){
			model<-paste(model, "Y",i, "~0*1\n", sep="")		
		}
		
		## set the variance for Y
		for (i in 1:T){
			model<-paste(model, "Y",i, "~~start(", varey, ")*", "Y", i, " + varey*Y",i, "\n", sep="")		
		}
		model
	}
	
	sem.est <- function(model, data){
		temp.res <- sem(model=model, data=data)
		label <- temp.res@ParTable$label
		c(temp.res@ParTable$est[label!=""], temp.res@ParTable$se[label!=""])
	}
	## do it once for a given N and T
	fit.once <- function(N, T){
		## generate data
		pop.model.T <- pop.model(T)
		pop.model.T.res <- sem(pop.model.T, do.fit=FALSE)
		pop.model.T.cov <- inspect(pop.model.T.res, "cov.ov")
		pop.model.T.mean <- inspect(pop.model.T.res, "mean.ov")
		ynames <- row.names(pop.model.T.cov)
		gen.data <- lapply(1:R, mvrnorm, n=N, mu=pop.model.T.mean, Sigma=pop.model.T.cov)
		
		## conduct the analysis
		fit.model.T <- fit.model(T)
		fit.res <- lapply(gen.data, sem.est, model=fit.model.T)
		
		## run once to get the model information
		model.info.res <- sem(fit.model.T, gen.data[[1]])
		label <- model.info.res@ParTable$label
		label <- label[label!=""]
		label.unique <- !duplicated(label)
		label <- label[label.unique]
		npar <- length(label)
		## get the parameter estimates, sd, se, power, CI of power
		all.res <- do.call(rbind, fit.res)
		all.res <- all.res[, c(label.unique, label.unique)]
		all.res <- na.omit(all.res)
		mc.est <- colMeans(all.res[, 1:npar], na.rm=TRUE)
		mc.se <- apply(all.res[, (npar+1):(2*npar)], 2, mean, na.rm=TRUE)
		mc.sd <- apply(all.res[, 1:npar], 2, sd, na.rm=TRUE)
		
		mc.z.score <- all.res[, 1:npar]/all.res[, (npar+1):(2*npar)]
		mc.z.score.check <- abs(mc.z.score) >= qnorm(1-alpha/2)
		
		mc.power <- colMeans(mc.z.score.check, na.rm=TRUE)
		
		pop.par <- unlist(lapply(label, function(x){eval(parse(text=x))}))
		mc.output <- cbind(pop.par, mc.est, mc.sd, mc.se, mc.power, N, T)
		row.names(mc.output) <- label
		label.sort <- sort(label)
		mc.output[label.sort, ]
	}
	
	if (length(N)>1 | length(T)>1){
		all.output <- list()
		for (i in N){
			for (j in T){
				all.output [[paste('N',i,'-T',j, sep="")]]<- fit.once(i,j)
			}
		}
	}else{
		all.output <- fit.once(N,T)
	}
	class(all.output) <- "lcs.power"
	all.output 	
}

plot.lcs.power <- function(x, parameter, ...){
	## x is the output from power analysis
	power.mat <- do.call('rbind', x)
	power.par <- power.mat[rownames(power.mat)==parameter, ]
	unique.N <- unique(power.par[ ,6])
	unique.T <- unique(power.par[ ,7])
	
	if (length(unique.N)==1 & length(unique.T)==1) stop("Multiple N or T is needed for power plot.")
	
	if (length(unique.N)==1){
		## plot the power along T
		plot(power.par[, 7], power.par[, 5], type='l', xlab='Number of Occasions', ylab=paste('Power of ',parameter, sep=""), ylim=c(0,1))
		points(power.par[, 7], power.par[, 5])
	}
	
	if (length(unique.T)==1){
		plot(power.par[, 6], power.par[, 5], type='l', xlab='Sample size', ylab=paste('Power of ',parameter, sep=""), ylim=c(0,1))
		points(power.par[, 6], power.par[, 5])
	}
	
	if (length(unique.N)>1 & length(unique.T)>1){
		for (N in unique.N){
			## plot power with time for a given sample size
			temp.power <- power.par[power.par[, 6]==N, ]
			plot(temp.power[, 7], temp.power[, 5], type='l', xlab='Number of Occasions', ylab=paste('Power of ',parameter, sep=""), ylim=c(0,1))
			points(temp.power[, 7], temp.power[, 5])
			cat ("Press [enter] to continue")
			line <- readline()
		}
		
		for (T in unique.T){
			## plot power with time for a given sample size
			temp.power <- power.par[power.par[, 7]==T, ]
			plot(temp.power[, 6], temp.power[, 5], type='l', xlab='Sample size', ylab=paste('Power of ',parameter, sep=""), ylim=c(0,1))
			points(temp.power[, 6], temp.power[, 5])
			cat ("Press [enter] to continue")
			line <- readline()
		}
	}
}

wp.blcsm<-function(N=100, T=5, R=1000,
	betay=0, my0=0, mys=0, varey=1, vary0=1, varys=1, vary0ys=0, alpha=0.05,
	betax=0, mx0=0, mxs=0, varex=1, varx0=1, varxs=1, varx0xs=0, varx0y0=0,  
	varx0ys=0, vary0xs=0, varxsys=0, gammax=0, gammay=0, ...){
	
	pop.model <- function(T){
		## for y
		## latent y
		## Intercept
		model<-"y0 =~ 1*y1\n"
		
		## path from y(t-1) to y(t) with path 1
		for (i in 2:T){
			model<-paste(model, "y",i,"~1*y",(i-1),"\n", sep="")
		}
		
		## loading from dy(t) to y(t) with path 1
		for (i in 2:T){
			model<-paste(model, "dy",i,"=~1*y",i,"\n", sep="")
		}
		
		## path from y(t) to dy(t+1) with path betay
		for (i in 2:T){
			model<-paste(model, "dy",i,"~", betay, "*y", (i-1), "\n", sep="")
		}
		
		## latent slope ys factor model
		for (i in 2:T){
			model<-paste(model, "ys=~1*dy", i, "\n", sep="")
		}
		
		## variance for dy constraints to 0
		for (i in 2:T){
			model<-paste(model, "dy",i,"~~0*dy",i,"\n", sep="")
		}
		
		## variance for y constraints to 0
		for (i in 1:T){
			model<-paste(model, "y",i,"~~0*y",i,"\n", sep="")
		}
		
		## variance and covariance for intercept and slope
		model<-paste(model, "ys~~", vary0ys, "*y0\n", sep="")
		model<-paste(model, "y0~~", vary0, "*y0\n", sep="")
		model<-paste(model, "ys~~", varys, "*ys\n", sep="")
		
		model<-paste(model, "ys~", mys, "*1\n", sep="")	
		model<-paste(model, "y0~", my0, "*1\n", sep="")

		## constrain means of y and dy to be zero
		for (i in 1:T){
			model<-paste(model, "y",i,"~0*1\n", sep="")
		}
		for (i in 2:T){
			model<-paste(model, "dy",i,"~0*1\n", sep="")
		}
		
		
		## for observed data part
		## y(t) to Y(t)
		for (i in 1:T){
			model<-paste(model, "y",i,"=~1*", "Y",i, "\n", sep="")				
		}
		
		## set means of Y to be zero
		for (i in 1:T){
			model<-paste(model, "Y",i, "~0*1\n", sep="")		
		}
		
		## set the variance for Y
		for (i in 1:T){
			model<-paste(model, "Y",i, "~~", varey, "*", "Y",i, "\n", sep="")		
		}
		
		
		## for x
		## latent x
		## Intercept
		model<-paste(model, "x0 =~ 1*x1\n")
		
		## path from x(t-1) to x(t) with path 1
		for (i in 2:T){
			model<-paste(model, "x",i,"~1*x",(i-1),"\n", sep="")
		}
		
		## loading from dx(t) to x(t) with path 1
		for (i in 2:T){
			model<-paste(model, "dx",i,"=~1*x",i,"\n", sep="")
		}
		
		## path from x(t) to dx(t+1) with path betax
		for (i in 2:T){
			model<-paste(model, "dx",i,"~", betax, "*x", (i-1), "\n", sep="")
		}
		
		## latent slope xs factor model
		for (i in 2:T){
			model<-paste(model, "xs=~1*dx", i, "\n", sep="")
		}
		
		## variance for dx constraints to 0
		for (i in 2:T){
			model<-paste(model, "dx",i,"~~0*dx",i,"\n", sep="")
		}
		
		## variance for x constraints to 0
		for (i in 1:T){
			model<-paste(model, "x",i,"~~0*x",i,"\n", sep="")
		}
		
		## variance and covariance for intercept and slope
		model<-paste(model, "xs~~", varx0xs, "*x0\n", sep="")
		model<-paste(model, "x0~~", varx0, "*x0\n", sep="")
		model<-paste(model, "xs~~", varxs, "*xs\n", sep="")
		
		model<-paste(model, "xs~", mxs, "*1\n", sep="")	
		model<-paste(model, "x0~", mx0, "*1\n", sep="")

		## constrain means of x and dx to be zero
		for (i in 1:T){
			model<-paste(model, "x",i,"~0*1\n", sep="")
		}
		for (i in 2:T){
			model<-paste(model, "dx",i,"~0*1\n", sep="")
		}
		
		
		## for observed data part
		## x(t) to X(t)
		for (i in 1:T){
			model<-paste(model, "x",i,"=~1*", "X",i, "\n", sep="")				
		}
		
		## set means of X to be zero
		for (i in 1:T){
			model<-paste(model, "X",i, "~0*1\n", sep="")		
		}
		
		## set the variance for X
		for (i in 1:T){
			model<-paste(model, "X",i, "~~", varex, "*", "X",i, "\n", sep="")		
		}
	
	
	
		## coupling effects
		for (i in 2:T){
			model<-paste(model, "dy",i,"~", gammax, "*x",i-1, "\n", sep="")
		}
	
		for (i in 2:T){
			model<-paste(model, "dx",i,"~", gammay, "*y",i-1, "\n", sep="")
		}
	
		model<-paste(model, "x0~~", varx0y0, "*y0\n", sep="")	
		model<-paste(model, "x0~~", varx0ys, "*ys\n", sep="")	
		model<-paste(model, "y0~~", vary0xs, "*xs\n", sep="")
		model<-paste(model, "xs~~", varxsys, "*ys\n", sep="")
		model
	}
	
	fit.model <- function(T){
		## latent y
		## Intercept
		model<-"y0 =~ 1*y1\n"
		
		## path from y(t-1) to y(t) with path 1
		for (i in 2:T){
			model<-paste(model, "y",i,"~1*y",(i-1),"\n", sep="")
		}
		
		## loading from dy(t) to y(t) with path 1
		for (i in 2:T){
			model<-paste(model, "dy",i,"=~1*y",i,"\n", sep="")
		}
		
		## path from y(t) to dy(t+1) with path betay
		for (i in 2:T){
			model<-paste(model, "dy",i,"~start(", betay, ")*y", (i-1)," + betay*y", (i-1), "\n", sep="")
		}
		
		## latent slope ys factor model
		for (i in 2:T){
			model<-paste(model, "ys=~1*dy", i, "\n", sep="")
		}
		
		## variance for dy constraints to 0
		for (i in 2:T){
			model<-paste(model, "dy",i,"~~0*dy",i,"\n", sep="")
		}
		
		## variance for y constraints to 0
		for (i in 1:T){
			model<-paste(model, "y",i,"~~0*y",i,"\n", sep="")
		}
		
		## variance and covariance for intercept and slope
		model<-paste(model, "ys~~start(", vary0ys, ")*y0 + vary0ys*y0\n", sep="")
		model<-paste(model, "y0~~start(", vary0, ")*y0 + vary0*y0\n", sep="")
		model<-paste(model, "ys~~start(", varys, ")*ys + varys*ys\n", sep="")
		
		model<-paste(model, "ys~start(", mys, ")*1 + label('mys')*1\n", sep="")	
		model<-paste(model, "y0~start(", my0, ")*1 + label('my0')*1\n", sep="")

		## constrain means of y and dy to be zero
		for (i in 1:T){
			model<-paste(model, "y",i,"~0*1\n", sep="")
		}
		for (i in 2:T){
			model<-paste(model, "dy",i,"~0*1\n", sep="")
		}
		
		
		## for observed data part
		## y(t) to Y(t)
		for (i in 1:T){
			model<-paste(model, "y",i,"=~1*", "Y",i, "\n", sep="")				
		}
		
		## set means of Y to be zero
		for (i in 1:T){
			model<-paste(model, "Y",i, "~0*1\n", sep="")		
		}
		
		## set the variance for Y
		for (i in 1:T){
			model<-paste(model, "Y",i, "~~start(", varey, ")*", "Y", i, " + varey*Y",i, "\n", sep="")		
		}
		
		## latent x
		## Intercept
		model<-paste(model, "x0 =~ 1*x1\n")
		
		## path from x(t-1) to x(t) with path 1
		for (i in 2:T){
			model<-paste(model, "x",i,"~1*x",(i-1),"\n", sep="")
		}
		
		## loading from dx(t) to x(t) with path 1
		for (i in 2:T){
			model<-paste(model, "dx",i,"=~1*x",i,"\n", sep="")
		}
		
		## path from x(t) to dx(t+1) with path betax
		for (i in 2:T){
			model<-paste(model, "dx",i,"~start(", betax, ")*x", (i-1)," + betax*x", (i-1), "\n", sep="")
		}
		
		## latent slope xs factor model
		for (i in 2:T){
			model<-paste(model, "xs=~1*dx", i, "\n", sep="")
		}
		
		## variance for dx constraints to 0
		for (i in 2:T){
			model<-paste(model, "dx",i,"~~0*dx",i,"\n", sep="")
		}
		
		## variance for x constraints to 0
		for (i in 1:T){
			model<-paste(model, "x",i,"~~0*x",i,"\n", sep="")
		}
		
		## variance and covariance for intercept and slope
		model<-paste(model, "xs~~start(", varx0xs, ")*x0 + varx0xs*x0\n", sep="")
		model<-paste(model, "x0~~start(", varx0, ")*x0 + varx0*x0\n", sep="")
		model<-paste(model, "xs~~start(", varxs, ")*xs + varxs*xs\n", sep="")
		
		model<-paste(model, "xs~start(", mxs, ")*1 + label('mxs')*1\n", sep="")	
		model<-paste(model, "x0~start(", mx0, ")*1 + label('mx0')*1\n", sep="")

		## constrain means of x and dx to be zero
		for (i in 1:T){
			model<-paste(model, "x",i,"~0*1\n", sep="")
		}
		for (i in 2:T){
			model<-paste(model, "dx",i,"~0*1\n", sep="")
		}
		
		
		## for observed data part
		## x(t) to X(t)
		for (i in 1:T){
			model<-paste(model, "x",i,"=~1*", "X",i, "\n", sep="")				
		}
		
		## set means of X to be zero
		for (i in 1:T){
			model<-paste(model, "X",i, "~0*1\n", sep="")		
		}
		
		## set the variance for X
		for (i in 1:T){
			model<-paste(model, "X",i, "~~start(", varex, ")*", "X", i, " + varex*X",i, "\n", sep="")		
		}
		
		## coupling effects
		for (i in 2:T){
			model<-paste(model, "dy",i,"~start(", gammax, ")*x", i-1, " + gammax*x", i-1, "\n", sep="")
		}
	
		for (i in 2:T){
			model<-paste(model, "dx",i,"~start(", gammay, ")*y", i-1, " + gammay*y", i-1, "\n", sep="")
		}
		
		model<-paste(model, "x0~~start(", varx0y0, ")*y0 + varx0y0*y0\n", sep="")	
		model<-paste(model, "x0~~start(", varx0ys, ")*ys + varx0ys*ys\n", sep="")	
		model<-paste(model, "y0~~start(", vary0xs, ")*xs + vary0xs*xs\n", sep="")
		model<-paste(model, "xs~~start(", varxsys, ")*ys + varxsys*ys\n", sep="")
		
		model
	}
	
	sem.est <- function(model, data){
		temp.res <- sem(model=model, data=data)
		label <- temp.res@ParTable$label
		c(temp.res@ParTable$est[label!=""], temp.res@ParTable$se[label!=""])
	}
	## do it once for a given N and T
	fit.once <- function(N, T){
		## generate data
		pop.model.T <- pop.model(T)
		pop.model.T.res <- sem(pop.model.T, do.fit=FALSE)
		pop.model.T.cov <- inspect(pop.model.T.res, "cov.ov")
		pop.model.T.mean <- inspect(pop.model.T.res, "mean.ov")
		ynames <- row.names(pop.model.T.cov)
		gen.data <- lapply(1:R, mvrnorm, n=N, mu=pop.model.T.mean, Sigma=pop.model.T.cov)
		
		## conduct the analysis
		fit.model.T <- fit.model(T)
		fit.res <- lapply(gen.data, sem.est, model=fit.model.T)
		
		## run once to get the model information
		model.info.res <- sem(fit.model.T, gen.data[[1]])
		label <- model.info.res@ParTable$label
		label <- label[label!=""]
		label.unique <- !duplicated(label)
		label <- label[label.unique]
		npar <- length(label)
		## get the parameter estimates, sd, se, power, CI of power
		all.res <- do.call(rbind, fit.res)
		all.res <- all.res[, c(label.unique, label.unique)]
		all.res <- na.omit(all.res)
		mc.est <- colMeans(all.res[, 1:npar], na.rm=TRUE)
		
		mc.se <- apply(all.res[, (npar+1):(2*npar)], 2, mean, na.rm=TRUE)
		mc.sd <- apply(all.res[, 1:npar], 2, sd, na.rm=TRUE)
		
		mc.z.score <- all.res[, 1:npar]/all.res[, (npar+1):(2*npar)]
		mc.z.score.check <- abs(mc.z.score) >= qnorm(1-alpha/2)
		
		mc.power <- colMeans(mc.z.score.check, na.rm=TRUE)
		
		pop.par <- unlist(lapply(label, function(x){eval(parse(text=x))}))
		mc.output <- cbind(pop.par, mc.est, mc.sd, mc.se, mc.power, N, T)
		row.names(mc.output) <- label
		label.sort <- sort(label)
		mc.output[label.sort, ]
	}
	
	if (length(N)>1 | length(T)>1){
		all.output <- list()
		for (i in N){
			for (j in T){
				all.output [[paste('N',i,'-T',j, sep="")]]<- fit.once(i,j)
			}
		}
	}else{
		all.output <- fit.once(N,T)
	}
	class(all.output) <- "lcs.power"
	all.output 	
}

summary.power<-function(object, ...) {

    # Show some basic about the simulation methods

    # main part: parameter estimates
    cat("Basic information:\n\n")
	if(object$out@Options$test %in% c("satorra.bentler", "yuan.bentler",
                                  "mean.var.adjusted",
                                  "scaled.shifted") &&
       length(object$out@Fit@test) > 1L) {
        scaled <- TRUE
        if(object$out@Options$test == "scaled.shifted")
            shifted <- TRUE
        else
            shifted <- FALSE
    } else {
        scaled <- FALSE
        shifted <- FALSE
    }
	
	t0.txt <- sprintf("  %-40s", "Esimation method")
    t1.txt <- sprintf("  %10s", object$out@Options$estimator)
    t2.txt <- ifelse(scaled, 
              sprintf("  %10s", "Robust"), "")
    cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
	
	t0.txt <- sprintf("  %-40s", "Standard error")
    t1.txt <- sprintf("  %10s", object$out@Options$se)
	cat(t0.txt, t1.txt, "\n", sep="")
	
	if (!is.null(object$info$bootstrap)){
		t0.txt <- sprintf("  %-40s", "Number of requested bootstrap")
		t1.txt <- sprintf("  %10i", object$info$bootstrap)
		cat(t0.txt, t1.txt, "\n", sep="")
	}

    t0.txt <- sprintf("  %-40s", "Number of requested replications")
    t1.txt <- sprintf("  %10i", object$info$nrep)
    cat(t0.txt, t1.txt, "\n", sep="")
    t0.txt <- sprintf("  %-40s", "Number of successful replications")
    t1.txt <- sprintf("  %10i", nrow(object$results$estimates))
    cat(t0.txt, t1.txt, "\n", sep="")
		
    cat("\n")

    # local print function
    print.estimate <- function(name="ERROR", i=1, z.stat=TRUE) {       
        # cut name if (still) too long
        name <- strtrim(name, width=15L)
        txt <- sprintf("    %-15s %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
                               name, pop.value[i], est[i], mse[i], sd.est[i], power[i], cvg[i])
        cat(txt)
    }

    est <- apply(object$results$estimates, 2, mean, na.rm=TRUE)
    sd.est  <- apply(object$results$estimates, 2, sd, na.rm=TRUE)
	mse <- apply(object$results$se, 2, mean, na.rm=TRUE)
	power <- object$power
	pop.value<-object$pop.value
	cvg<-object$coverage
	
	object<-object$out
	

    for(g in 1:object@Data@ngroups) {
        ov.names <- lavNames(object, "ov")
        lv.names <- lavNames(object, "lv")

        # group header
        if(object@Data@ngroups > 1) {
            if(g > 1) cat("\n\n")
            cat("Group ", g, 
                " [", object@Data@group.label[[g]], "]:\n\n", sep="")
        }

        # estimates header
		#txt <- sprintf("%-13s %12s %12s %8s %8s %8s\n", "", "True", "Estimate", "MSE", "SD", "Power")
		#cat(txt)
        cat("                       True  Estimate      MSE      SD     Power Coverage\n")
		
        makeNames <- function(NAMES, LABELS) {
            multiB <- FALSE
            if(any(nchar(NAMES) != nchar(NAMES, "bytes")))
                multiB <- TRUE
            if(any(nchar(LABELS) != nchar(LABELS, "bytes")))
                multiB <- TRUE
            # labels?
            l.idx <- which(nchar(LABELS) > 0L)
            if(length(l.idx) > 0L) {
                if(!multiB) {
                    LABELS <- abbreviate(LABELS, 4)
                    LABELS[l.idx] <- paste(" (", LABELS[l.idx], ")", sep="")
                    MAX.L <- max(nchar(LABELS))
                    NAMES <- abbreviate(NAMES, minlength = (13 - MAX.L), 
                                        strict = TRUE)
                } else {
                    # do not abbreviate anything (eg in multi-byte locales)
                    MAX.L <- 4L
                }
                NAMES <- sprintf(paste("%-", (13 - MAX.L), "s%", MAX.L, "s",
                                       sep=""), NAMES, LABELS)
            } else {
                if(!multiB) {
                    NAMES <- abbreviate(NAMES, minlength = 13, strict = TRUE)
                } else {
                    NAMES <- sprintf(paste("%-", 13, "s", sep=""), NAMES)
                }
            }

            NAMES
        }

        NAMES <- object@ParTable$rhs

        # 1a. indicators ("=~") (we do show dummy indicators)
        mm.idx <- which( object@ParTable$op == "=~" & 
                        !object@ParTable$lhs %in% ov.names &
                         object@ParTable$group == g)
        if(length(mm.idx)) {
            cat("Latent variables:\n")
            lhs.old <- ""
            NAMES[mm.idx] <- makeNames(  object@ParTable$rhs[mm.idx],
                                       object@ParTable$label[mm.idx])
            for(i in mm.idx) {
                lhs <- object@ParTable$lhs[i]
                if(lhs != lhs.old) cat("  ", lhs, " =~\n", sep="")
                print.estimate(name=NAMES[i], i)
                lhs.old <- lhs
            }
            cat("\n")
        }

        # 1b. formative/composites ("<~")
        fm.idx <- which( object@ParTable$op == "<~" &
                         object@ParTable$group == g)
        if(length(fm.idx)) {
            cat("Composites:\n")
            lhs.old <- ""
            NAMES[fm.idx] <- makeNames(  object@ParTable$rhs[fm.idx],
                                       object@ParTable$label[fm.idx])
            for(i in fm.idx) {
                lhs <- object@ParTable$lhs[i]
                if(lhs != lhs.old) cat("  ", lhs, " <~\n", sep="")
                print.estimate(name=NAMES[i], i)
                lhs.old <- lhs
            }
            cat("\n")
        }

        # 2. regressions
        eqs.idx <- which(object@ParTable$op == "~" & object@ParTable$group == g)
        if(length(eqs.idx) > 0) {
            cat("Regressions:\n")
            lhs.old <- ""
            NAMES[eqs.idx] <- makeNames(  object@ParTable$rhs[eqs.idx],
                                        object@ParTable$label[eqs.idx])
            for(i in eqs.idx) {
                lhs <- object@ParTable$lhs[i]
                if(lhs != lhs.old) cat("  ", lhs, " ~\n", sep="")
                print.estimate(name=NAMES[i], i)
                lhs.old <- lhs
            }
            cat("\n")
        }

        # 3. covariances
        cov.idx <- which(object@ParTable$op == "~~" & 
                         !object@ParTable$exo &
                         object@ParTable$lhs != object@ParTable$rhs &
                         object@ParTable$group == g)
        if(length(cov.idx) > 0) {
            cat("Covariances:\n")
            lhs.old <- ""
            NAMES[cov.idx] <- makeNames(  object@ParTable$rhs[cov.idx],
                                        object@ParTable$label[cov.idx])
            for(i in cov.idx) {
                lhs <- object@ParTable$lhs[i]
                if(lhs != lhs.old) cat("  ", lhs, " ~~\n", sep="")
                print.estimate(name=NAMES[i], i)
                lhs.old <- lhs
            }
            cat("\n")
        }

        # 4. intercepts/means
        ord.names <- lavNames(object, "ov.ord")
        int.idx <- which(object@ParTable$op == "~1" & 
                         !object@ParTable$lhs %in% ord.names &
                         !object@ParTable$exo &
                         object@ParTable$group == g)
        if(length(int.idx) > 0) {
            cat("Intercepts:\n")
            NAMES[int.idx] <- makeNames(  object@ParTable$lhs[int.idx],
                                        object@ParTable$label[int.idx])
            for(i in int.idx) {
                print.estimate(name=NAMES[i], i)
            }
            cat("\n")
        }

        # 4b thresholds
        th.idx <- which(object@ParTable$op == "|" &
                        object@ParTable$group == g)
        if(length(th.idx) > 0) {
            cat("Thresholds:\n")
            NAMES[th.idx] <- makeNames(  paste(object@ParTable$lhs[th.idx],
                                               "|",
                                               object@ParTable$rhs[th.idx],
                                               sep=""),
                                         object@ParTable$label[th.idx])
            for(i in th.idx) {
                print.estimate(name=NAMES[i], i)
            }
            cat("\n")
        }

        # 5. (residual) variances
        var.idx <- which(object@ParTable$op == "~~" &
                         !object@ParTable$exo &
                         object@ParTable$lhs == object@ParTable$rhs &
                         object@ParTable$group == g)
        if(length(var.idx) > 0) {
            cat("Variances:\n")
            NAMES[var.idx] <- makeNames(  object@ParTable$rhs[var.idx],
                                        object@ParTable$label[var.idx])
            for(i in var.idx) {
                if(object@Options$mimic == "lavaan") {
                    print.estimate(name=NAMES[i], i, z.stat=FALSE)
                } else {
                    print.estimate(name=NAMES[i], i, z.stat=TRUE)
                }
            }
            cat("\n")
        }

        # 6. latent response scales
        delta.idx <- which(object@ParTable$op == "~*~" &
                         object@ParTable$group == g)
        if(length(delta.idx) > 0) {
            cat("Scales y*:\n")
            NAMES[delta.idx] <- makeNames(  object@ParTable$rhs[delta.idx],
                                            object@ParTable$label[delta.idx])
            for(i in delta.idx) {
                print.estimate(name=NAMES[i], i, z.stat=TRUE)
            }
            cat("\n")
        }

    } # ngroups

    # 6. variable definitions
    def.idx <- which(object@ParTable$op == ":=")
    if(length(def.idx) > 0) {
        if(object@Data@ngroups > 1) cat("\n")
        cat("Indirect/Mediation effects:\n")
        NAMES[def.idx] <- makeNames(  object@ParTable$lhs[def.idx], "")
        for(i in def.idx) {
            print.estimate(name=NAMES[i], i)
        }
        cat("\n")
    }

}

wp.popPar<-function(object){
	par.value<-object@ParTable$ustart
	for (i in 1:length(par.value)){
		if (is.na(par.value[i])){
			if (object@ParTable$op[i]=="~~"){
				par.value[i]<-1
			}else{
				par.value[i]<-0
			}
		}
	}
	names(par.value)<-object@ParTable$label
	## for indirect effec defined here
	for (i in 1:length(par.value)){
		if (object@ParTable$op[i]==":="){
			temp<-object@ParTable$rhs[i]
			temp<-gsub(' +','',temp)
			temp<-gsub('-','+', temp, fixed=TRUE)
			temp<-unlist(strsplit(temp, '+', fixed=TRUE))
			m<-length(temp)
			temp.est<-0
			par<-NULL
			for (j in 1:m){
				temp1<-unlist(strsplit(temp[j], '*', fixed=TRUE))
				par<-c(par, temp1)
			}
			ind.exp<-parse(text=object@ParTable$rhs[i])
			par.list<-as.list(par.value[par])
			par.value[i]<-eval(ind.exp, par.list)
		}
	}
	par.value
}

wp.mc.sem.basic<-function(model, indirect=NULL, nobs=100, nrep=1000, alpha=.95, skewness=NULL, kurtosis=NULL, ovnames=NULL, se="default", estimator="default", parallel="no", ncore=Sys.getenv('NUMBER_OF_PROCESSORS'), cl=NULL, ...){
	if (missing(model)) stop("A model is needed.")
	
	model.indirect<-paste(model, "\n", indirect, "\n")
	ngroups <- length(nobs)
	## Initial analysis for some model information
	newdata<-simulateData(model,sample.nobs=nobs,skewness=0,kurtosis=0, ...)			

	if (ngroups > 1){
		temp.res<-sem(model.indirect, data=newdata, se=se, estimator=estimator, group='group', ...)
	}else{
		temp.res<-sem(model.indirect, data=newdata, se=se, estimator=estimator, ...)
	}
	par.value<-wp.popPar(temp.res)
	idx <- 1:length( temp.res@ParTable$lhs )
	#cnames<-paste(temp.res@ParTable$lhs[idx], temp.res@ParTable$op[idx], temp.res@ParTable$rhs[idx])
	ov<-lavNames(temp.res,'ov')
	## dealing with skewness and kurtosis
	if (length(skewness) > 1){
		if (length(skewness)!=length(ovnames)) stop("The number of skewness is not equal to the number observed variables")	
		index<-match(ov, ovnames)
		skewness<-skewness[index]
	}
	
	if (length(kurtosis) > 1){
		if (length(kurtosis)!=length(ovnames)) stop("The number of kurtosis is not equal to the number observed variables")	
		index<-match(ov, ovnames)
		kurtosis<-kurtosis[index]
	}	
	newdata<-simulateData(model,sample.nobs=nobs,skewness=skewness,kurtosis=kurtosis, ...)		
	
	cnames<-paste(temp.res@ParTable$lhs[idx], temp.res@ParTable$op[idx], temp.res@ParTable$rhs[idx])
	out<-temp.res
	runonce<-function(i){
		## Step 1: generate data
		newdata<-simulateData(model,sample.nobs=nobs,skewness=skewness,kurtosis=kurtosis, ...)		
		## Step 2: fit the model 		
		if (ngroups > 1){
			temp.res<-try(sem(model.indirect, data=newdata, se=se, estimator=estimator, group='group', warn=FALSE, ...))
		}else{
			temp.res<-try(sem(model.indirect, data=newdata, se=se, estimator=estimator, warn=FALSE, ...))
		}
		
		## Step 3: Check significance
		
		if (!inherits(temp.res, "try-error")){
		idx <- 1:length( temp.res@ParTable$lhs )
		temp.est<-temp.res@Fit@est[idx]
		temp.se<-temp.res@Fit@se[idx]
		temp.sig<-temp.est/temp.se
		crit<-qnorm(1-(1-alpha)/2)
		temp.sig.ind<-abs(temp.sig)>crit
		temp.sig.ind[!is.finite(temp.sig)]<-NA
		
		## Step 4: Check the coverage
		ci.u<-temp.est+crit*temp.se
		ci.l<-temp.est-crit*temp.se
		temp.cvg<- (ci.u>par.value[idx]) & (ci.l<par.value[idx])
		}else{
			nna<-length(idx)
			temp.est<-rep(NA, nna)
			temp.sig.ind<-rep(NA, nna)
			temp.se<-rep(NA, nna)
			temp.cvg<-rep(NA, nna)
		}		
		return(list(temp.sig.ind=temp.sig.ind, temp.est=temp.est, temp.se=temp.se, temp.cvg=temp.cvg))
	}
	
	## run parallel or not
	# this is from the boot function in package boot
	old_options <- options(); options(warn = -1)
    have_mc <- have_snow <- FALSE
    ncore <- ncore
    if (parallel != "no" && ncore > 1L) {
        if (parallel == "multicore") have_mc <- .Platform$OS.type != "windows"
        else if (parallel == "snow") have_snow <- TRUE
        if (!have_mc && !have_snow) ncore <- 1L
    }
	
    RR <- nrep
    res <- if (ncore > 1L && (have_mc || have_snow)) {
        if (have_mc) {
            parallel::mclapply(seq_len(RR), runonce, mc.cores = ncore)
        } else if (have_snow) {
            list(...) # evaluate any promises
            if (is.null(cl)) {
                cl <- parallel::makePSOCKcluster(rep("localhost", ncore))
                if(RNGkind()[1L] == "L'Ecuyer-CMRG")
                    parallel::clusterSetRNGStream(cl)
                res <- parallel::parLapply(cl, seq_len(RR), runonce)
                parallel::stopCluster(cl)
                res
            } else parallel::parLapply(cl, seq_len(RR), runonce)
        }
    } else lapply(seq_len(RR), runonce)
	
	all.sig<-do.call(rbind, lapply(res, "[[", 'temp.sig.ind'))

	all.par<-do.call(rbind, lapply(res, "[[", 'temp.est'))
	all.se<-do.call(rbind, lapply(res, "[[", 'temp.se'))
	all.cvg<-do.call(rbind, lapply(res, "[[", 'temp.cvg'))
		
	options(old_options)	
	
	colnames(all.sig)<-cnames
	power<-apply(all.sig, 2, mean, na.rm=TRUE)
	cvg<-apply(all.cvg, 2, mean, na.rm=TRUE)
	info<-list(nobs=nobs, nrep=nrep, alpha=alpha, method="Normal", bootstrap=NULL)
	print(power)
	object<-list(power=power, coverage=cvg, pop.value=par.value, results=list(estimates=all.par, se=all.se, all=res), info=info, out=out, data=newdata)
	class(object)<-'power'
	return(object)
}

wp.mc.sem.boot<-function(model, indirect=NULL, nobs=100, nrep=1000, nboot=1000, alpha=.95, skewness=NULL, kurtosis=NULL, ovnames=NULL, se="default", estimator="default", parallel="no", ncore=Sys.getenv('NUMBER_OF_PROCESSORS'), cl=NULL, ...){		
	if (missing(model)) stop("A model is needed.")	
	
	## internal function
	coef.new<-function(x,...){
		lavaan::coef(x, type='user', ...)
	}
	model.indirect<-paste(model, "\n", indirect, "\n")
	ngroups <- length(nobs)
	## Initial analysis for some model information
	newdata<-simulateData(model,sample.nobs=nobs,skewness=0,kurtosis=0, ...)			
	if (ngroups > 1){
		temp.res<-sem(model.indirect, data=newdata, se=se, estimator=estimator, group='group', warn=FALSE, ...)
	}else{
		temp.res<-sem(model.indirect, data=newdata, se=se, estimator=estimator, warn=FALSE, ...)
	}
	par.value<-wp.popPar(temp.res)
	idx <- 1:length( temp.res@ParTable$lhs )
	#cnames<-paste(temp.res@ParTable$lhs[idx], temp.res@ParTable$op[idx], temp.res@ParTable$rhs[idx])
	ov<-lavNames(temp.res,'ov')
	ptype<-temp.res@ParTable$op
	
	## dealing with skewness and kurtosis
	if (length(skewness) > 1){
		if (length(skewness)!=length(ovnames)) stop("The number of skewness is not equal to the number observed variables")	
		index<-match(ov, ovnames)
		skewness<-skewness[index]
	}
	
	if (length(kurtosis) > 1){
		if (length(kurtosis)!=length(ovnames)) stop("The number of kurtosis is not equal to the number observed variables")	
		index<-match(ov, ovnames)
		kurtosis<-kurtosis[index]
	}
	newdata<-simulateData(model,sample.nobs=nobs,skewness=skewness,kurtosis=kurtosis, ...)		
	ci.bc1<-function(x, b, cl=.95){
		n<-length(x)
		z0<-qnorm(sum(x<b, na.rm=TRUE)/n)
		alpha<-(1-cl)/2
		alpha<-c(alpha, 1-alpha)
		alpha<-sort(alpha)
		alpha1<-alpha
		alpha<-pnorm(2*z0+qnorm(alpha))
		dig <- max(2L, getOption("digits"))
		np<-length(alpha)
		qs<-quantile(x, alpha, na.rm=TRUE)
		names(qs) <- paste(if (np < 100) 
            formatC(100 * alpha1, format = "fg", width = 1, digits = dig)
        else format(100 * alpha1, trim = TRUE, digits = dig), 
            "%", sep = "")
		qs
	}
	
	ci.bc<-function(par.boot, par0, cl=.95){
		se.boot<-apply(par.boot, 2, sd, na.rm=TRUE)
		estimate<-par0
		p<-ncol(par.boot)
		ci<-NULL
		for (i in 1:p){
			ci<-rbind(ci, ci.bc1(par.boot[,i], par0[i], cl))
		}
		cbind(estimate, se.boot, ci)
	}
	
	## repeated the following for nrep times
	runonce<-function(i){		
		## Step 1: generate data
		newdata<-simulateData(model,sample.nobs=nobs,skewness=skewness,kurtosis=kurtosis, ...)	
		
		## Step 2: fit the model 		
		#temp.res<-sem(model.indirect, data=newdata, se=se, estimator=estimator, ...)
		if (ngroups > 1){
			temp.res<-try(sem(model.indirect, data=newdata, se=se, estimator=estimator, group='group', warn=FALSE, ...))
		}else{
			temp.res<-try(sem(model.indirect, data=newdata, se=se, estimator=estimator, warn=FALSE, ...))
		}
		
		
		## Step 3: Conduct bootstrap analysis
		if (!inherits(temp.res, "try-error")){
		orig.res<-coef.new(temp.res)
		boot.res<-bootstrapLavaan(temp.res, FUN=coef.new, R=nboot, ...)
		ci.res<-ci.bc(boot.res, orig.res, cl=alpha)
		## Step 4: Check the coverage		
		 
		temp.cvg<- (ci.res[,4]>=par.value[idx]) & (ci.res[,3]<par.value[idx])
		
		## Step 5: check significance
		temp.sig<-NULL
		temp.est<-coef.new(temp.res) 
		temp.se<-ci.res[,2]
		for (jj in 1:length(ptype)){
			if (ptype[jj]=="~~"){
				crit<-qnorm(1-(1-alpha)/2)
				ci.temp<-ci.res[jj, 1] + c(-1,1)*ci.res[jj, 2]*crit
				temp.sig<-c(temp.sig, (ci.temp[1]>0 | ci.temp[2]<0))
			}else{
				temp.sig<-c(temp.sig, (ci.res[jj, 3]>0 | ci.res[jj, 4]<0))
			}
		}
		}else{
			nna<-length(idx)
			temp.est<-rep(NA, nna)
			temp.sig<-rep(NA, nna)
			temp.se<-rep(NA, nna)
			temp.cvg<-rep(NA, nna)
			boot.res<-NA
			ci.res<-NA
		}

		return(list(temp.sig=temp.sig, temp.est=temp.est, temp.se=temp.se, temp.cvg=temp.cvg, boot.res=boot.res, ci.res=ci.res))
	}
	
	## run parallel or not
	# this is from the boot function in package boot
	old_options <- options(); options(warn = -1)
    have_mc <- have_snow <- FALSE
    ncore <- ncore
    if (parallel != "no" && ncore > 1L) {
        if (parallel == "multicore") have_mc <- .Platform$OS.type != "windows"
        else if (parallel == "snow") have_snow <- TRUE
        if (!have_mc && !have_snow) ncore <- 1L
    }
	
    RR <- nrep
    res <- if (ncore > 1L && (have_mc || have_snow)) {
        if (have_mc) {
            parallel::mclapply(seq_len(RR), runonce, mc.cores = ncore)
        } else if (have_snow) {
            list(...) # evaluate any promises
            if (is.null(cl)) {
                cl <- parallel::makePSOCKcluster(rep("localhost", ncore))
                if(RNGkind()[1L] == "L'Ecuyer-CMRG")
                    parallel::clusterSetRNGStream(cl)
                res <- parallel::parLapply(cl, seq_len(RR), runonce)
                parallel::stopCluster(cl)
                res
            } else parallel::parLapply(cl, seq_len(RR), runonce)
        }
    } else lapply(seq_len(RR), runonce)
	
	all.sig<-do.call(rbind, lapply(res, "[[", 'temp.sig'))

	all.par<-do.call(rbind, lapply(res, "[[", 'temp.est'))
	all.se<-do.call(rbind, lapply(res, "[[", 'temp.se'))
	all.cvg<-do.call(rbind, lapply(res, "[[", 'temp.cvg'))
		
	options(old_options)	
	
	#colnames(all.sig)<-cnames
	power<-apply(all.sig, 2, mean, na.rm=TRUE)
	cvg<-apply(all.cvg, 2, mean, na.rm=TRUE)
	info<-list(nobs=nobs, nrep=nrep, alpha=alpha, method="Normal", bootstrap=nboot)
	## print(power)
	object<-list(power=power, coverage=cvg, pop.value=par.value, results=list(estimates=all.par, se=all.se, all=res), info=info, out=temp.res, data=newdata)
	class(object)<-'power'
	return(object)
}


wp.mc.sem.power.curve<-function(model, indirect=NULL, nobs=100, type='basic', nrep=1000, nboot=1000, alpha=.95, skewness=NULL, kurtosis=NULL, ovnames=NULL, se="default", estimator="default", parallel="no", ncore=Sys.getenv('NUMBER_OF_PROCESSORS'), cl=NULL, ...){		
	if (missing(model)) stop("A model is needed.")	
	
	allpower<-NULL
	
	## check whether nobs is a vector or not
	if (is.vector(nobs)){	
		for (N in nobs){
			if (type == 'basic'){
				## sobel test based analysis
				indpower <- wp.mc.sem.basic(model=model, indirect=indirect, nobs=N, nrep=nrep, alpha=alpha, skewness=skewness, kurtosis=kurtosis, ovnames=ovnames, se=se, estimator=estimator, parallel=parallel, ncore=ncore, cl=cl, ...)			
			}else{
				indpower <- wp.mc.sem.boot(model=model, indirect=indirect, nobs=N, nrep=nrep, nboot=nboot, alpha=alpha, skewness=skewness, kurtosis=kurtosis, ovnames=ovnames, se=se, estimator=estimator, parallel=parallel, ncore=ncore, cl=cl, ...)	
			}
			allpower<-rbind(allpower, indpower$power)
		}
	}else{
		for (k in 1:nrow(nobs)){
			N <- nobs[k]
			if (type == 'basic'){
				## sobel test based analysis
				indpower <- wp.mc.sem.basic(model=model, indirect=indirect, nobs=N, nrep=nrep, alpha=alpha, skewness=skewness, kurtosis=kurtosis, ovnames=ovnames, se=se, estimator=estimator, parallel=parallel, ncore=ncore, cl=cl, ...)			
			}else{
				indpower <- wp.mc.sem.boot(model=model, indirect=indirect, nobs=N, nrep=nrep, nboot=nboot, alpha=alpha, skewness=skewness, kurtosis=kurtosis, ovnames=ovnames, se=se, estimator=estimator, parallel=parallel, ncore=ncore, cl=cl, ...)	
			}
			allpower<-rbind(allpower, indpower$power)
		}
	}

	
	#windows(record=TRUE)
	#op <- par(ask=TRUE)
    #on.exit(par(op))
	pnames <- colnames(allpower)
	for (j in 1:ncol(allpower)){
		if (sum(is.nan(allpower[,j])) == 0){
			if (is.vector(nobs)){
			plot(nobs, allpower[, j], ylab=paste('Power of', pnames[j]), xlab='Sample size')
			lines(nobs, allpower[, j])
			}else{
			N<-apply(nobs, 1, sum)
			plot(N, allpower[, j], ylab=paste('Power of', pnames[j]), xlab='Sample size')
			lines(N, allpower[, j])
			}
		}
	}
	invisible(allpower)
}

## This function performs the sample size calculation for a mixed model of repeated measures with AR(1) correlation structure. See Lu, Luo, & Chen (2008) for parameter definitions and other details.
## power.mmrm.ar1 {longpower}

wp.mmrm.ar1 <-
function (N = NULL, rho = NULL, ra = NULL, sigmaa = NULL, rb = NULL, 
    sigmab = NULL, lambda = 1, times = 1:length(ra), delta = NULL, 
    alpha = 0.05, power = NULL, alternative = c("two.sided", 
        "one.sided")) 
{
    if (sum(sapply(list(N, rho, delta, power, alpha), is.null)) != 
        1) 
        stop("exactly one of 'N', 'rho', 'delta', 'power', and 'alpha' must be NULL")
    if (!is.null(alpha) && !is.numeric(alpha) || any(0 > 
        alpha | alpha > 1)) 
        stop("'alpha' must be numeric in [0, 1]")
    alternative <- match.arg(alternative)
    if (length(times) != length(ra)) 
        stop("ra and times should be the same length")
    phia <- phib <- NULL
    n.body <- quote({
        J <- length(ra)
        phia <- 1/ra[J] - sum(rho^(2 * (times[J] - times[-J])) * 
            (1/ra[-1] - 1/ra[-J]))
        if (!is.null(rb)) {
            phib <- 1/rb[J] - sum(rho^(2 * (times[J] - times[-J])) * 
                (1/rb[-1] - 1/rb[-J]))
        } else {
            rb <- ra
            phib <- phia
        }
        if (is.null(sigmab)) {
            sigma <- sigmaa
        } else {
            sigma <- mean(c(sigmaa, sigmab))
        }
        
        Na <-lambda*N/(1+lambda)    
        pnorm( delta/sigma*sqrt(Na/((phia + lambda * phib))) +  qnorm(ifelse(alternative == "two.sided", alpha/2, alpha)) )
    })
    
    
    if (is.null(alpha)) 
        alpha <- uniroot(function(alpha) eval(n.body) - 
            power, c(1e-10, 1 - 1e-10))$root
    else if (is.null(delta)) 
        delta <- uniroot(function(delta) eval(n.body) - 
            power, 
            c(1e-10, 1e+05))$root
    else if (is.null(N)) uniroot(function(N) eval(n.body) - 
            power, c(10, 1e+05))$root
    else if (is.null(rho)) 
        rho <- uniroot(function(rho) eval(n.body) - power, c(1e-10, 
            1 - 1e-10))$root
    else power <- eval(n.body)
	
	NOTE <- "Based power.mmrm.ar1 in longpower (Lu, Luo, & Chen ,2008)"
  
  METHOD <- "Power for Mixed Model of Repeated Measures"
  URL <- "http://psychstat.org/longar"
	
    structure(list(n=N, n1 = lambda*N/(1+lambda), n2 = N/(1+lambda), rho = rho,  phi1 = phia, phi2 = phib, 
        delta = delta, alpha = alpha, 
        power = power, alternative = alternative, note = NOTE, 
        method = METHOD, url = URL), 
        class = "webpower")
}  

## This function performs the sample size calculation for a mixed model of repeated measures with general correlation structure. See Lu, Luo, & Chen (2008) for parameter definitions and other details. This function executes Formula (3) on page 4.
## power.mmrm {longpower}

wp.mmrm<-function (N = NULL, Ra = NULL, ra = NULL, 
sigmaa = NULL, Rb = NULL, rb = NULL, sigmab = NULL, 
lambda = 1, delta = NULL, alpha = 0.05, 
power = NULL, alternative = c("two.sided", "one.sided")) 
{
    if (sum(sapply(list(N, delta, power, alpha), is.null)) != 
        1) 
        stop("exactly one of 'N', 'delta', 'power', and 'alpha' must be NULL")
    if (!is.null(alpha) && !is.numeric(alpha) || any(0 > 
        alpha | alpha > 1)) 
        stop("'alpha' must be numeric in [0, 1]")
    alternative <- match.arg(alternative)
    n.body <- quote({
        Ia <- 0
        ra <- c(ra, 0)
        for (j in 1:nrow(Ra)) {
            Raj <- matrix(0, nrow(Ra), nrow(Ra))
            Raj[1:j, 1:j] <- solve(Ra[1:j, 1:j])
            Ia <- Ia + (ra[j] - ra[j + 1]) * Raj
        }
        phia <- solve(Ia)[j, j]
        if (is.null(rb)) rb <- ra
        if (is.null(Rb)) Rb <- Ra
        Ib <- 0
        rb <- c(rb, 0)
        for (j in 1:nrow(Rb)) {
            Rbj <- matrix(0, nrow(Rb), nrow(Rb))
            Rbj[1:j, 1:j] <- solve(Rb[1:j, 1:j])
            Ib <- Ib + (rb[j] - rb[j + 1]) * Rbj
        }
        phib <- solve(Ib)[j, j]
        if (is.null(sigmab)) {
            sigma <- sigmaa
        } else {
            sigma <- mean(c(sigmaa, sigmab))
        }

        Na <-lambda*N/(1+lambda)    
        pnorm( delta/sigma*sqrt(Na/((phia + lambda * phib))) +  qnorm(ifelse(alternative == "two.sided", alpha/2, alpha)) )    
    })
    if (is.null(alpha)) 
        alpha <- uniroot(function(alpha) eval(n.body) - 
            power, c(1e-10, 1 - 1e-10))$root
    else if (is.null(delta)) 
        delta <- uniroot(function(delta) eval(n.body) - 
            power, 
            c(1e-10, 1e+05))$root
    else if (is.null(N)) uniroot(function(N) eval(n.body) - 
            power, c(10, 1e+05))$root
    else power <- eval(n.body)
    
    NOTE <- "Based power.mmrm in longpower (Lu, Luo, & Chen ,2008)"
  
  METHOD <- "Power for Mixed Model of Repeated Measures"
  URL <- "http://psychstat.org/longmmrm"
    structure(list(n=N, n1 = lambda*N/(1+lambda), n2 = N/(1+lambda), delta = delta, alpha = alpha, power = power, alternative = alternative,  note = NOTE, method = METHOD, url = URL), class = "webpower")
}

wp.mc.t <- function(n = NULL, R0 = 1e+5, R1=1e+3, mu0 = 0, mu1 = 0, sd = 1, 
    skewness = 0, kurtosis = 3, alpha = 0.05,  
    type = c("two.sample", "one.sample", "paired"), alternative = c("two.sided", 
        "less", "greater")) {
    
    if (sum(sapply(list(n, mu0, R0, R1, mean, sd, skewness, kurtosis, alpha), 
        is.null)) >= 1) 
        stop("n, mu0, mean, sd, skewness, kurtosis, and alpha cannot be NULL")
    
    if (!is.null(alpha) && !is.numeric(alpha) || any(0 > alpha | 
        alpha > 1)) 
        stop(sQuote("alpha"), " must be numeric in [0, 1]")
    
    type <- match.arg(type)
    alternative <- match.arg(alternative)
    tsample <- switch(type, one.sample = 1, two.sample = 2, paired = 1)
    tside <- switch(alternative, less = 1, two.sided = 2, greater = 3)
    
    if (tsample == 1) {
        ##### under null hypothesis
        moments0 <- c(mean = mu0, variance = sd^2, skewness = skewness, 
            kurtosis = kurtosis)
        data0 <- matrix(rpearson(n * R0, moments = moments0), ncol = n, 
            nrow = R0)
        y_bar0 <- apply(data0, 1, mean)
        var0 <- apply(data0, 1, var)
        t0 <- (y_bar0 - mu0)/sqrt(var0/n)
        ##### under alternative hypothesis
        momentsa <- c(mean = mu1, variance = sd^2, skewness = skewness, 
            kurtosis = kurtosis)
        dataa <- matrix(rpearson(n * R1, moments = momentsa), ncol = n, 
            nrow = R1)
        y_bara <- apply(dataa, 1, mean)
        vara <- apply(dataa, 1, var)
        ta <- (y_bara - mu0)/sqrt(vara/n)
    } else if (tsample == 2) {
        if (length(mu0) == 1) 
            mu0 <- rep(mu0, 2)
		if (length(mu1) == 1) 
            mu1 <- rep(mu1, 2)
        if (length(sd) == 1) 
            sd <- rep(sd, 2)
        if (length(skewness) == 1) 
            skewness <- rep(skewness, 2)
        if (length(kurtosis) == 1) 
            kurtosis <- rep(kurtosis, 2)
        if (length(n) == 1) 
            n <- rep(n, 2)
        ##### under null hypothesis
        moments01 <- c(mean = mu0[1], variance = sd[1]^2, skewness = skewness[1], 
            kurtosis = kurtosis[1])
        moments02 <- c(mean = mu0[2], variance = sd[2]^2, skewness = skewness[2], 
            kurtosis = kurtosis[2])
        data01 <- matrix(rpearson(n[1] * R0, moments = moments01), ncol = n[1], 
            nrow = R0)
        data02 <- matrix(rpearson(n[2] * R0, moments = moments02), ncol = n[2], 
            nrow = R0)
        y_bar01 <- apply(data01, 1, mean)
        var01 <- apply(data01, 1, var)
        y_bar02 <- apply(data02, 1, mean)
        var02 <- apply(data02, 1, var)
        t0 <- (y_bar01 - y_bar02)/sqrt(var01/n[1] + var02/n[2])
        ##### under alternative hypothesis
        momentsa1 <- c(mean = mu1[1], variance = sd[1]^2, skewness = skewness[1], 
            kurtosis = kurtosis[1])
        momentsa2 <- c(mean = mu1[2], variance = sd[2]^2, skewness = skewness[2], 
            kurtosis = kurtosis[2])
        dataa1 <- matrix(rpearson(n[1] * R1, moments = momentsa1), ncol = n[1], 
            nrow = R1)
        dataa2 <- matrix(rpearson(n[2] * R1, moments = momentsa2), ncol = n[2], 
            nrow = R1)
        y_bara1 <- apply(dataa1, 1, mean)
        vara1 <- apply(dataa1, 1, var)
        y_bara2 <- apply(dataa2, 1, mean)
        vara2 <- apply(dataa2, 1, var)
        ta <- (y_bara1 - y_bara2)/sqrt(vara1/n[1] + vara2/n[2])
    }
    
    
    if (tside == 2) {
        t.cri <- quantile(sort(t0), c(alpha/2, 1 - alpha/2))
        power <- mean(ta >= t.cri[2] | ta <= t.cri[1])
    }
    
    if (tside == 1) {
        t.cri <- quantile(sort(t0), alpha)
        power <- mean(ta <= t.cri)
    }
    
    if (tside == 3) {
        t.cri <- quantile(sort(t0), 1 - alpha)
        power <- mean(ta >= t.cri)
    }
    
    NOTE <- switch(type, paired = "n is number of *pairs*", two.sample = "n is number in *each* group", 
        NULL)
    METHOD <- paste(switch(type, one.sample = "One-sample", two.sample = "Two-sample", 
        paired = "Paired"), "t test power calculation")
	URL <- "http://psychstat.org/tnonnormal"
    if (tsample == 2) {
        structure(list(n1 = n[1], n2 = n[2], power = power, mean1 = mu1[1], mean2 = mu1[2], sd1 = sd[1], sd2 = sd[2], skewness1 = skewness[1], skewness2 = skewness[2], 
            kurtosis1 = kurtosis[1], kurtosis2 = kurtosis[1], alpha = alpha, 
            alternative = alternative, note = NOTE, method = METHOD, url=URL), 
            class = "webpower")
    } else if (tsample == 1) {
        structure(list(n = n, power = power, mu0=mu0, mu1 = mu1, sd = sd, skewness = skewness, 
            kurtosis = kurtosis, alpha = alpha,
            alternative = alternative, note = NOTE, method = METHOD, url=URL), class = "webpower")
    }
} 


wp.mc.chisq.diff <- function(full.model.pop, full.model, reduced.model, N=100, R=1000, alpha=0.05){
  
  chi.diff <- rep(0, R)

  for (i in 1:R){
    sim.data <- simulateData(full.model.pop, sample.nobs=N)
    full.res <- sem(full.model, data=sim.data)
    reduced.res <- sem(reduced.model, data=sim.data)
  
    chi.diff[i] <- reduced.res@Fit@test[[1]]$stat - full.res@Fit@test[[1]]$stat
  }
  df <- reduced.res@Fit@test[[1]]$df - full.res@Fit@test[[1]]$df
  power <- mean(chi.diff > qchisq(1-alpha, df))
  list(power=power, df=df, chi.dff=chi.diff)
  
  structure(list(n = N, power = power, df=df, alpha = alpha), class = "webpower")
}

sem.effect.size <- function(full.model.pop, reduced.model){
  N <- 1000
  full.res<-sem(full.model.pop, do.fit=FALSE)
  sigma.F<-fitted.values(full.res)$cov
  reduced.res<-sem(reduced.model, sample.cov=sigma.F, sample.nobs=N)
  delta <- reduced.res@Fit@test[[1]]$stat/N 
  df <- reduced.res@Fit@test[[1]]$df
  RMSEA <- fitMeasures(reduced.res, "rmsea")
  list(delta=delta, df=df, RMSEA=RMSEA)
}