# Reproduce the normal-model reference intervals in the coverage vignette.
# Run from the SMSD source root: Rscript data-raw/normal_reference.R
# These fixed-sample reference methods are not exported SMSD procedures.
# Pearson CDF: centered bivariate-normal geometry, integrated by quadrature.
# R-squared CDF: Oliveira and Ferreira (2022), Equation 2.3, a negative-binomial
# mixture of beta CDFs. Mixing-tail truncation error is at most 2e-13.
# Normal observations and independent random predictors are assumed.

quadrature <- function(shape, order) {
  j <- seq_len(order)
  matrix <- diag(2*(j-1)+shape)
  off <- sqrt((1:(order-1)) * ((1:(order-1))+shape-1))
  matrix[cbind(1:(order-1),2:order)] <- off
  matrix[cbind(2:order,1:(order-1))] <- off
  e <- eigen(matrix,symmetric=TRUE)
  list(nodes=e$values,weights=e$vectors[1,]^2)
}
cache <- new.env(parent=emptyenv())
normal_r_cdf <- function(r,rho,n,order=96L) {
  if(r<=-1) return(0)
  if(r>=1) return(1)
  t <- r*sqrt((n-2)/(1-r*r))
  key<-paste(n,order)
  if(!exists(key,cache,inherits=FALSE)) {
    s<-quadrature((n-1)/2,order)
    v<-quadrature((n-2)/2,order)
    assign(key,list(s=sqrt(2*s$nodes),v=sqrt(2*v$nodes/(n-2)),weights=outer(s$weights,v$weights)),cache)
  }
  q<-get(key,cache,inherits=FALSE)
  sum(stats::pnorm(outer(rho*q$s/sqrt(1-rho*rho),t*q$v,function(a,b)b-a))*q$weights)
}
q_r <- function(prob,rho,n) stats::uniroot(function(r)normal_r_cdf(r,rho,n)-prob,c(-1,1),tol=1e-9)$root
ci_r <- function(r,n) c(
  stats::uniroot(function(rho)normal_r_cdf(r,rho,n)-.975,c(-.99999,.99999),tol=1e-9)$root,
  stats::uniroot(function(rho)normal_r_cdf(r,rho,n)-.025,c(-.99999,.99999),tol=1e-9)$root)

r2_cdf <- function(r2,rho2,n,p=20L) {
  if(r2<=0) return(0)
  if(r2>=1) return(1)
  if(rho2<=0) return(stats::pbeta(r2,p/2,(n-p-1)/2))
  if(rho2>=1) return(0)
  shape<-(n-1)/2
  low<-max(0,stats::qnbinom(1e-13,size=shape,prob=1-rho2))
  high<-stats::qnbinom(1-1e-13,size=shape,prob=1-rho2)
  j<-seq.int(low,high)
  sum(stats::dnbinom(j,size=shape,prob=1-rho2)*stats::pbeta(r2,p/2+j,(n-p-1)/2))
}
q_r2 <- function(prob,rho2,n) stats::uniroot(function(r)r2_cdf(r,rho2,n)-prob,c(0,1),tol=1e-10)$root
ci_r2 <- function(r,n) {
  c0<-r2_cdf(r,0,n)
  lower<-if(c0<=.975) 0 else stats::uniroot(function(pop)r2_cdf(r,pop,n)-.975,c(0,1),tol=1e-9)$root
  upper<-if(c0<=.025) 0 else stats::uniroot(function(pop)r2_cdf(r,pop,n)-.025,c(0,1),tol=1e-9)$root
  c(lower,upper)
}

# Maximum widths at adjacent sample sizes; the R-squared maximum is numerical.
pearson_width <- function(n) {
  t <- stats::qt(.975, n - 1)
  2 * t / sqrt(n - 1 + t * t)
}
pearson_bounds <- data.frame(n = 95:96,
                            maximum_width = vapply(95:96, pearson_width, 0))
r2_bounds <- do.call(rbind, lapply(127:128, function(n) {
  optimum <- stats::optimize(function(r) diff(ci_r2(r, n)),
                             c(.1, .9), maximum = TRUE, tol = 1e-7)
  grid <- seq(.1, .9, length.out = 301)
  grid_width <- vapply(grid, function(r) diff(ci_r2(r, n)), 0)
  # Monotone confidence endpoints bound the width in the two remaining tails.
  tail_upper <- ci_r2(.1, n)[2]
  tail_lower <- 1 - ci_r2(.9, n)[1]
  stopifnot(max(grid_width) <= optimum$objective + 1e-7,
            max(tail_upper, tail_lower) < .3)
  data.frame(n = n, maximum_width = optimum$objective,
             sample_R2_at_maximum = optimum$maximum)
}))
stopifnot(pearson_bounds$maximum_width[1] > .4,
          pearson_bounds$maximum_width[2] < .4,
          r2_bounds$maximum_width[1] > .3,
          r2_bounds$maximum_width[2] < .3)

# Normal quadrature agrees with the central t CDF when rho = 0.
r_control <- seq(-.5, .5, length.out = 21)
stopifnot(max(vapply(r_control, function(r)
  abs(normal_r_cdf(r, 0, 96) -
        stats::pt(r * sqrt(94 / (1 - r * r)), 94)), 0)) < 1e-7)
# Independently confirm the maximum Pearson width on a positive-r profile.
profile <- vapply(seq(0, .99, length.out = 101),
                  function(r) diff(ci_r(r, 96)), 0)
stopifnot(all(diff(profile) < 1e-7),
          abs(max(profile) - pearson_width(96)) < 1e-7)

# Independent raw normal data, not draws from either proposed reference CDF.
n <- 96L
G <- 100000L
rlo <- q_r(.025, .8, n)
rhi <- q_r(.975, .8, n)
set.seed(84026L)
r <- unlist(lapply(1:50, function(i) {
  x <- matrix(stats::rnorm(n * 2000), nrow = n)
  y <- .8 * x + .6 * matrix(stats::rnorm(n * 2000), nrow = n)
  x <- sweep(x, 2, colMeans(x))
  y <- sweep(y, 2, colMeans(y))
  colSums(x * y) / sqrt(colSums(x * x) * colSums(y * y))
}))
covered <- r >= rlo & r <= rhi
mc <- stats::binom.test(sum(covered), G)$conf.int
pearson <- data.frame(parameter = "Pearson correlation = 0.8", n = n, G = G,
  nominal = .95, theoretical_coverage = .95, coverage = mean(covered),
  mc_lower = mc[1], mc_upper = mc[2], maximum_width = pearson_width(n),
  target_width = .4)

n <- 128L
p <- 20L
G <- 20000L
set.seed(85026L)
r2 <- replicate(G, {
  x <- matrix(stats::rnorm(n * p), nrow = n)
  y <- stats::rnorm(n)
  fit <- stats::lm.fit(cbind(1, x), y)
  1 - sum(fit$residuals^2) / sum((y - mean(y))^2)
})
covered <- r2 <= stats::qbeta(.975, p / 2, (n - p - 1) / 2)
mc <- stats::binom.test(sum(covered), G)$conf.int
r2_result <- data.frame(parameter = "R-squared = 0, 20 predictors", n = n, G = G,
  nominal = .95, theoretical_coverage = .975, coverage = mean(covered),
  mc_lower = mc[1], mc_upper = mc[2],
  maximum_width = r2_bounds$maximum_width[r2_bounds$n == n], target_width = .3)
results <- rbind(pearson, r2_result)
print(results, row.names = FALSE)
saveRDS(list(summary = results, pearson_bounds = pearson_bounds,
             R2_bounds = r2_bounds,
             meta = list(R = R.version.string,
                         generated_on = as.character(Sys.Date()),
                         seeds = c(84026L, 85026L),
                         generator = "data-raw/normal_reference.R",
                         method = "Fixed-sample equal-tailed normal-model CDF inversion",
                         pearson_quadrature_order = 96L,
                         R2_mixing_tail_bound = 2e-13)),
        "inst/extdata/normal_reference.rds", compress = "xz")
