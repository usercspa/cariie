### "The forecasting of menstruation based on a state-space modeling 
### of basal body temperature time series",
### 
### Supplementary R script for non-Gaussian filtering and prediction, along with a
### simulated menstrual cycle dataset

### Functions ----------------------------------------------------------

## log_sum_exp()
##      This function returns the natural logarithm of the sum of the exponentials.
##      Arguments
##          x: vector of real values
log_sum_exp <- function(x) {
    log(sum(exp(x)))
}

## normalize()
##      This function executes normalization and returns log probability densities.
##      Arguments
##          log_density: (unnormalized) log density
normalize <- function(log_density) {
    max_log_density <- max(log_density)
    log_density_scaled <- log_density - max_log_density
    log_density_scaled - log_sum_exp(log_density_scaled) + log(length(log_density))
}

## marginalize()
##      This function executes marginalization and returns log marginal probability densities.
##      Arguments
##          log_density: (joint) log conditional density p(omeaga_t, omega_{t-1})
##          direction: If "forward", p(omeaga_t, omega_{t-1}) becomes p(omeaga_t).
##              If "backward", p(omeaga_t, omega_{t-1}) becomes p(omeaga_{t-1}).
marginalize <- function(log_density, direction = c("forward", "backward")) {
    if (direction == "forward")
        log(apply(exp(log_density), 1, sum) / nrow(log_density))
    else
        log(apply(exp(log_density), 2, sum) / nrow(log_density))
}

## lerch()
##      This function evaluates Lerch's transcedental function.
##      Arguments
##          z, s, a: arguments of the function
##          limit: upper limit of summation
lerch <- function(z, s, a, limit) {
    k <- 0:limit
    len <- length(a)
    out <- vector(length = len)
    for (i in 1:len) out[i] <- sum(z^k / (a[i] + k)^s)
    out
}

## dwgamma()
##      This function returns probability density of the wrapped Gamma distribution
##          defined on [0, 1). The origin of PDF is fixed at 0.
##      Arguments
##          x: vector of locations where probability densities are evaluated
##          alpha: shape parameter
##          beta: rate parameter
##          log: logical; if TRUE log densities are returned
##          limit: upper limit of summation for the lerch() function
dwgamma <- function(x, alpha, beta, log = TRUE, limit = 1000) {
    log_pdf <- alpha * log(beta) - lgamma(alpha) - 
        beta * x + 
        log(lerch(exp(-beta), 1 - alpha, x, limit))

    # If alpha < 1 and the result of numerical integration of PDF is 
    #     significantly smaller than 1, 
    #     log_pdf that is the closest to 0 is compensated.
    if (alpha < 1) {
        if (sum(exp(log_pdf) / length(x)) < 1 - 1e-8)
            log_pdf[1] <- 
                log(exp(log_pdf[1]) + (1 - sum(exp(log_pdf) / length(x))) * length(x))
    }

    vlog <- rep(log, length(x))
    ifelse(vlog, log_pdf, exp(log_pdf))
}

## filter()
##      This function executes one-step-ahead prediction and filtering.
##          It returns log likelihood (loglik), log probability densities of the 
##          predictive distribution (log_pp), and that of filtering distribution (log_pf).
##      Arguments
##          log_pf_previous: log probability densities of filtering distribution at time t-1
##          y_new: observation of BBT at time t
##          z_new: observation of menstruation at time t
##          log_pm: system model log probability densities
##          parameters: vector of model parameters (alpha, beta, sigma, a, b, c)
##          M: order of the trigonometric series for the BBT likelihood
filter <- function(log_pf_previous, y_new, z_new, log_pm, parameters, M) {

    ## Set parameters
    sigma <- parameters[3]
    a     <- parameters[4]
    b     <- parameters[5:(4 + M)]
    c     <- parameters[(4 + M + 1):(4 + 2 * M)]

    ## Define intervals
    N <- nrow(log_pf_previous)  # Number of grid points
    omega <- 0:(N - 1) / N      # Location of grid points

    ## Expected temperature
    mu <- a
    for(m in 1:M) {
        x <- 2 * pi * m * omega
        mu <- mu + b[m] * cos(x) + c[m] * sin(x)
    }

    ## Marginalization
    log_pf_marginal <- marginalize(log_pf_previous, "forward")

    ## Prediction
    log_pp <- t(t(log_pm) + log_pf_marginal)
    log_pp <- normalize(log_pp)

    ## Log-likelihood
    # BBT
    loglik_y <- matrix(0, N, N)
    if (!is.na(y_new)) {
        loglik_y <- matrix(rep(dnorm(y_new, mu, sigma, log = TRUE), N), N, N)
    }

    # Menstruation
    loglik_z <- matrix(0, N, N)
    if(z_new) {
        loglik_z[lower.tri(loglik_z)] <- log(0)
    } else {
        loglik_z[upper.tri(loglik_z)] <- log(0)
    }

    ## Filtering
    loglik_all <- log_pp + loglik_y + loglik_z
    max_loglik_all <- max(loglik_all)
    loglik_all_scaled <- loglik_all - max_loglik_all
    loglik <- log_sum_exp(loglik_all_scaled) - 2 * log(N) + max_loglik_all
    log_pf <- normalize(loglik_all)

    ## Output
    list(loglik = loglik, log_pp = log_pp, log_pf = log_pf)
}

## smoother()
##      This function executes smoothing and returns the log probability densities of 
##          the fixed-interval type smoothed distribution.
##      Arguments
##          log_ps_previous: log probability densities of smoothed distribution at time t+1
##          log_pf: log probability densities of filtering distribution at time t
smoother <- function(log_ps_previous, log_pf) {

    ## Define intervals
    N <- nrow(log_pf)

    ## Marginalization
    log_pf_marginal <- marginalize(log_pf, "forward")
    log_ps_prev_marginal <- marginalize(log_ps_previous, "backward")

    ## Replace numerical 0 with a small real value
    log_pf_marginal[is.infinite(log_pf_marginal)] <- log(.Machine$double.xmin)
    log_ps_prev_marginal[is.infinite(log_ps_prev_marginal)] <- log(.Machine$double.xmin)

    ## Smoothing
    log_ps <- log_pf + log_ps_prev_marginal - log_pf_marginal
    log_ps <- normalize(log_ps)

    # Output
    log_ps
}

## conditional_probabilities()
##      This function returns a matrix of conditional probabilities of 
##      menstruation day, f(k | omeaga_t).
##      Arguments
##          parameters: vector of model parameters (alpha, beta, sigma, a, b, c)
##          omega: location of grid points
##          k_max: upper limit of the day for which predictive probability is evaluated
conditional_probabilities <- function(parameters, omega, k_max = 100) {

    ## Set parameters
    alpha <- parameters[1]
    beta  <- parameters[2]

    ## Number of grid points
    N <- length(omega)

    ## Conditional probabilities f(k | omega_t)
    f <- matrix(nrow = k_max, ncol = N)
    for (i in 1:N) {
        f[1, ] <- 1 - pgamma(1 - omega, alpha, beta)
        for (k in 2:k_max) {
            f[k, ] <- 
                pgamma(1 - omega, (k - 1) * alpha, beta) - 
                pgamma(1 - omega, k * alpha, beta)
        }
    }

    # Output
    f
}

## predictor()
##      This function returns the predictive probability distribution of menstruation day.
##          It returns marginal predictive probabilities (marginal_probability) and
##          a point prediction for the menstruation day (point_prediction).
##      Arguments
##          log_pf: log probability densities of filtering distribution
##          parameters: vector of model parameters (alpha, beta, sigma, a, b, c)
##          k_max: upper limit of the day for which predictive probability is evaluated
predictor <- function(log_pf, f) {

    # Number of grid points
    N <- nrow(log_pf)

    # Marginal probabilities h(k | Y_t, Z_t)
    pf_marginal <- exp(marginalize(log_pf, "forward"))
    h <- c(f %*% pf_marginal / N)

    # Output
    list(marginal_probability = h, point_prediction = which.max(h))
}

## draw_figs()
##      This function draws conditional distributions of menstrual phase as well as
##          the predictive distribution of menstruation day.
##      Arguments
##          t: the day for which figures are drawn
##          f: a matrix of conditional probabilities of menstruation day
draw_figs <- function(t, f,p) {

    ## Conditional densities
    pred <- exp(marginalize(log_pp[, , t], "forward"))
    filt <- exp(marginalize(log_pf[, , t], "forward"))
    smth <- exp(marginalize(log_ps[, , t], "forward"))
    omega <- 0:(N - 1) / N

    y.max <- max(c(pred, filt, smth))
    plot(pred ~ omega, type = "n", main = sprintf("%sth day of the menstrual cycle", t),
         xlab = "Phase", ylab = "Probability density", ylim = c(0, y.max), lty = 2)
    polygon(x = c(omega, rev(omega)), y = c(rep(0, N), rev(smth)), col = "gray", border = FALSE)
    par(new = TRUE)
    plot(pred ~ omega, type = "l", ann = FALSE, axes = FALSE, ylim = c(0, y.max), lty = 2)
    par(new = TRUE)
    plot(filt ~ omega, type = "l", ylim = c(0, y.max), ann = FALSE, axes = FALSE)

    ## Predictive distribution for menstruation day
    res <- predictor(log_pf[, , t], f)
    days_passed <- t - 1

    new_res <- vector(length = 100)
    new_res[1:days_passed] <- 0
    new_res[(days_passed + 1):100] <- res$marginal_probability[1:(100 - days_passed)]
    new_pp <- res$point_prediction + days_passed
    if (p==TRUE) {
      plot(new_res, type = "h", xlim = c(0, 100), 
           ylim = c(0, max(new_res) + 0.002),
           main = sprintf("%sth day of the menstrual cycle", t),
           xlab = "Days since the onset of previous menstruation",
           ylab = "Probability")
      # Filled circle: Menstruation day
      points(x = T - 1, y = max(new_res) + 0.001, pch = 16)
      # Open triangle: model-based prediction
      points(x = new_pp, y = max(new_res) + 0.002, pch = 6)
      # Filled triangle: conventional prediction
      points(x = T-3, y = max(new_res) + 0.003, pch = 25, bg = "black")
      # Broken vertical line: Previous menstruation day
      abline(v = 0, lty = 2)
      # Solid vertical line: Current day
      abline(v = days_passed, lty = 1)
    }
    else {
      return(which.max(new_res))
    }
}

### Functions ----------------------------------------------------------


### Data and settings --------------------------------------------------

## BBT
y <- c(36.62, 36.02, 35.99, 36.40, 36.39, 36.27, 36.18, 36.45,
       35.83, 36.05, 35.89, 35.66, 36.03, 36.8, 35.64, 36.70,
       36.26, 36.07, 36.07, 36.43, 36.19, 35.90, 36.22, 36.42,
       36.48, 36.12, 36.21, 36.75, 36.84, 36.90, 36.92, 36.77,
       36.93, 36.15, 36.87, 36.59, 36.53, 36.56, 37.10)

## Menstruation
z <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
       0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
       0, 0, 0, 0, 0, 0, 0, 0, 1)

## Length of time series
T <- length(y)

## Parameters (subject 1)
parameters <- c(0.2098523,      # alpha
                8.915085,       # beta
                0.1118375,      # sigma
                36.2991,        # a
                0.192667921,    # b1
                0.014567663,    # b2
                0.008823447,    # b3
                -0.038979880,   # b4
                0.057591311,    # b5
                -0.013692024,   # b6
                -0.023629217,   # b7
                0.037039669,    # b8
                -0.085654152,   # b9
                0.061401905,    # b10
                -0.192667655,   # c1
                -0.051462140,   # c2
                -0.089428890,   # c3
                -0.029970492,   # c4
                -0.002739703,   # c5
                -0.052827804,   # c6
                0.067939328,    # c7
                -0.050495194,   # c8
                -0.027397176,   # c9
                0.024602075     # c10
                )

## Order of the trigonometric seris
M <- 10

## Number of grid points on the state space
N <- 2^9

## Location of grid points
omega <- 0:(N - 1) / N  

## log_pm: Log of probability densities of the system model
##      log_pm[i, j] = log(p(omega_t = omega[i] | omega_{t-1} = omega[j]))
log_pm <- matrix(nrow = N, ncol = N)
rownames(log_pm) <- colnames(log_pm) <- omega
log_system <- dwgamma(omega + 1 / N, parameters[1], parameters[2])
for(i in 1:(N - 1)) {
    log_pm[, i] <- c(log_system[(N - i + 1):N], log_system[1:(N - i)])
}
log_pm[, N] <- log_system

## f: Conditional probabilities of menstruation day
f <- conditional_probabilities(parameters, omega)

### Data and settings --------------------------------------------------


### State estimation and prediction ------------------------------------

## log_pp: log probability densities of the predictive distribution
##      log_pp[i, j, t] = log(p(omega_t = omega[i], omega_{t-1} = omega[j] | Y_{t-1}, Z_{t-1}))
## log_pf: log probability densities of the filtering distribution
##      log_pf[i, j, t] = log(p(omega_t = omega[i], omega_{t-1} = omega[j] | Y_t, Z_t))
## log_ps: log probability densities of the smoothed distribution
##      log_ps[i, j, t] = log(p(omega_t = omega[i], omega_{t-1} = omega[j] | Y_T, Z_T))
log_pp <- log_pf <- log_ps <- array(dim = c(N, N, T))
dimnames(log_pp)[[1]] <- dimnames(log_pp)[[2]] <- omega
dimnames(log_pp)[[3]] <- 1:T
dimnames(log_pf)[[1]] <- dimnames(log_pf)[[2]] <- omega
dimnames(log_pf)[[3]] <- 1:T
dimnames(log_ps)[[1]] <- dimnames(log_ps)[[2]] <- omega
dimnames(log_ps)[[3]] <- 1:T

## Filtering
log_pf_init <- matrix(log(1), N, N)     ## A uniform initial distribution
loglik <- 0

res <- filter(log_pf_init, y_new = y[1], z_new = z[1],
              log_pm = log_pm, parameters = parameters, M = M)
loglik <- loglik + res$loglik
log_pp[, , 1] <- res$log_pp
log_pf[, , 1] <- res$log_pf

for(t in 2:T) {
    res <- filter(log_pf[, , t - 1], y_new = y[t], z_new = z[t],
                  log_pm = log_pm, parameters = parameters, M = M)
    loglik <- loglik + res$loglik
    log_pp[, , t] <- res$log_pp
    log_pf[, , t] <- res$log_pf
}

## Smoothing
log_ps[, , T] <- log_pf[, , T]
for (t in (T - 1):1) {
    log_ps[, , t] <- smoother(log_ps[, , t + 1], log_pf[, , t])
}

## Draw figures
par(mfrow = c(2, 1))
draw_figs(1, f, TRUE) # true = plot, false = return next potential menstrualtion day.
pred_menstrual_day = draw_figs(1, f, FALSE)
pred_menstrual_day
### State estimation and prediction ------------------------------------

