library(MASS)
library(maxLik)
library(VGAM)

bsreg.fit = function(x, y, link = "log") {
  n <- NROW(x); p <- NCOL(x)
  linkobj <- make.link(link)
  linkfun <- linkobj$linkfun; linkinv <- linkobj$linkinv; Q.eta <- linkobj$mu.eta
  ystar <- linkfun(y)
  beta <- ginv(t(x) %*% x) %*% t(x) %*% ystar

  xbar <- mean(y)
  vart <- (n/(n-1)) * var(y)
  r <- vart / xbar^2
  alphai <- sqrt((2*r - 2 + 2*sqrt(1 + 3*r)) / (5 - r))

  if (is.nan(alphai) || is.na(alphai)) {
    s1 <- mean(y)
    r1 <- 1/mean(1/y)
    alphai <- sqrt(2*sqrt(s1/r1) - 2)
  }

  start <- c(as.vector(beta), alphai)

  # função de verossimilhança (log-like) - mantenha como antes
  fr <- function(vp) {
    betab <- vp[1:p]
    eta <- as.vector(x %*% betab)
    Q <- linkinv(eta)
    alphab <- vp[p+1]    # garantir escalar
    q <- 0.5
    zq <- qnorm(q)
    gma_alphab <- alphab * zq + sqrt(alphab^2 * zq^2 + 4)
    vt <- y
    sum(-0.5*log(8*pi*vt) - log(alphab) - log(gma_alphab) - 0.5*log(Q) +
        log(gma_alphab^2/2 + 2*Q/vt) -
        (2*Q/(alphab^2*gma_alphab^2*vt))*(vt*gma_alphab^2/(4*Q) - 1)^2)
  }

  # função gradiente: retorna um vetor NUMÉRICO de comprimento p+1
  grr <- function(vp) {
    betab <- vp[1:p]
    eta <- as.vector(x %*% betab)
    Q <- linkinv(eta)
    alphab <- vp[p+1]   # garantir escalar
    q <- 0.5
    zq <- qnorm(q)
    gma_alphab <- alphab * zq + sqrt(alphab^2 * zq^2 + 4)
    # derivadas auxiliares (mantive sua lógica)
    gma_alphabp <- zq + zq^2 * alphab * (1 / sqrt(alphab^2 * zq^2 + 4))
    # ... (outros termos idênticos aos seus)
    vt <- y
    z <- -0.5*(1/Q) - 2*(1/(alphab^2 * gma_alphab^2 * vt)) +
         gma_alphab^2 * vt * (1/(8 * alphab^2 * Q^2)) +
         4*(1/(vt * gma_alphab^2 + 4*Q))

    b <- -(gma_alphab + alphab * gma_alphabp) * (1/(alphab * gma_alphab)) +
         2*vt * gma_alphab * gma_alphabp * (1/(vt * gma_alphab^2 + 4*Q)) -
         (gma_alphab * gma_alphabp * alphab - gma_alphab^2) * vt * (1/(4*Q*alphab^3)) -
         2*(1/(alphab^3)) +
         4*Q*(gma_alphab + alphab * gma_alphabp) * (1/(alphab^3 * gma_alphab^3 * vt))

    # gradiente em relação a beta (p x 1)
    grad_beta <- as.vector(t(x) %*% (Q.eta(eta) * z))
    # gradiente em relação a alpha (escala)
    grad_alpha <- sum(b)

    # retornar UM VETOR numérico (p+1)
    c(grad_beta, grad_alpha)
  }

  # forço uso explícito do maxLik::maxBFGS
  opt <- maxLik::maxBFGS(fn = fr, grad = grr, start = start,
                         constraints = list(ineqA = matrix(c(rep(0,p),1),1,p+1),
                                            ineqB = 0))

  # inspeção rápida (útil para debug)
  # print(names(opt)); str(opt)

  if (!is.null(opt$code) && opt$code > 0) warning("optimization failed to converge (opt$code > 0)")
  if (!is.null(opt$convergence) && opt$convergence != 0) warning("optimizer signaled non-zero convergence code")

  estimates <- opt$estimate
  log.lik.est <- if (!is.null(opt$maximum)) opt$maximum else NA

  beta <- as.vector(estimates[1:p])
  eta <- as.vector(x %*% beta)
  Q <- linkinv(eta)
  alpha <- estimates[p+1]
  q <- 0.5

  # --- integrações: use Inf ao invés de 1e100
  # (mantenha suas funções f2,f3,f4,f5, mas com integrate(..., upper = Inf))
  # ... (reescreva essa parte mantendo a sua lógica, trocando 1e100 por Inf)

  # exemplo de correção rápida para integracao:
  # integ.f2[i] = integrate(f2, lower=0, upper=Inf, j=i)$value

  # calc de se, fisher etc. (mantive sua lógica - só atente-se aos indexadores)
  # ...
  # correção de zstatalpha:
  # zstatalpha = alpha / se[p+1]

  # Construir rval (como antes)
  rval <- list(coefficients = list(beta = beta, alpha = alpha),
               fitted.values.Q = structure(Q, .Names = names(y)),
               n = n, p = p, q = q, eta = eta,
               X = x, y = y,
               Hessian = if(!is.null(opt$hessian)) as.matrix(opt$hessian) else NULL,
               matrix.expected.Fisher = NA,
               loglik = log.lik.est,
               link = list(quantile = linkobj),
               converged = if(!is.null(opt$code)) (opt$code==0) else NA,
               information.criterions = list(aic = NA, bic = NA, aicc = NA)
               # ... (adicione demais campos conforme o seu fluxo)
               )

  return(rval)
}

model$coefficients
