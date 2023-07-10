## Code for 'Just another Gibbs additive modeller: Interfacing JAGS and mgcv'
## Final example also requires 'sitka.jags'.  Functions 'jagam', 'sim2gam' etc
## are in R package 'mgcv' >= 1.8-4 available from https://CRAN.R-project.org/.
## To run these examples you also need to have installed SemiPar and rjags from
## CRAN, and JAGS (http://mcmc-jags.sourceforge.net/) Note that at time of
## writing there is a problem with setting the RNG and seed using rjags, making
## it difficult to produce identical simulations except via a very ugly fix.
## (See example code in ?jagam, for how it should be done in theory.)


## Section 2 example from ?jagam from R package mgcv.

library("mgcv")
library("rjags")
load.module("glm")

set.seed(1)
n <- 400

dat <- gamSim(1, n = n, dist = "normal", scale = 2)
scale <- 0.5
Ey <- exp(dat$f/2)
dat$y <- rgamma(n, shape = 1/scale, scale = Ey * scale)

jd <- jagam(y ~ s(x0) + te(x1, x2) + s(x3), data = dat, family = Gamma(link = log), 
  file = "test.jags")

jm <- jags.model("test.jags", data = jd$jags.data, inits = jd$jags.ini, n.adapt = 2000, 
  n.chains = 1)

list.samplers(jm)
sam <- jags.samples(jm, c("b", "rho", "scale"), n.iter = 10000, thin = 10)
jam <- sim2jam(sam, jd$pregam)
par(mfrow = c(1, 3))
plot(jam)


## The section 4.1 example...

require(SemiPar)  ## for the data
require(mgcv)
require(rjags)
load.module("glm")
data(trade.union)  ## from SemiPar

## setup model...
jd <- jagam(union.member ~ s(wage, k = 20), data = trade.union, family = binomial, 
  file = "union.jags")

## compile and initialize...
jm <- jags.model("union.jags", data = jd$jags.data, inits = jd$jags.ini, n.chains = 1)
list.samplers(jm)

## sample....
system.time(sam <- jags.samples(jm, c("b", "rho", "mu"), n.iter = 10000, thin = 10))

## fake a reduced gam object ('jam' object)...
jam <- sim2jam(sam, jd$pregam)

## base plot...

plot(jam, shade = TRUE, shift = coef(jam)[1], trans = binomial()$linkinv, rug = FALSE, 
  ylim = c(-100, -coef(jam)[1]), seWithMean = TRUE, xlim = c(0, 30), lwd = 3)

## add response visualization...
nu <- trade.union$union.member == 0
with(trade.union, points(wage[nu], 0 * wage[nu], pch = 3, cex = 0.5))
with(trade.union, points(wage[!nu], 0 * wage[!nu] + 0.5, pch = 3, cex = 0.5))

## add some sample curves...
set.seed(0)
ii <- 1:20 * 50
pd <- data.frame(wage = 0:300/10)
Xp <- predict(jam, type = "lpmatrix", newdata = pd)
for (i in ii) {
  p <- binomial()$linkinv(Xp %*% sam$b[, i, 1])
  lines(pd$wage, p, lty = 2)
}


## The section 4.2 example (requires the supplied file 'sitka.jags', or edit
## auto-generated 'sitka0.jags' to produce this)...
library("SemiPar")
library("mgcv")
library("rjags")
load.module("glm")
data("sitka")
head(sitka)

## produce template model file 'sitka0,jags' and data...
jd <- jagam(log.size ~ s(days) + ozone, data = sitka, file = "sitka0.jags", diagonalize = TRUE)

## result edited to become sitka.jags, requiring id (id.num), and nd (length(id))
## to be added to data...

jd$jags.data$id <- sitka$id.num
jd$jags.data$nd <- length(unique(sitka$id.num))

## compile and initialize model...
jm <- jags.model("sitka.jags", data = jd$jags.data, inits = jd$jags.ini, n.chains = 1)

list.samplers(jm)

## sample...
sam <- jags.samples(jm, c("b", "rho", "scale", "mu"), n.iter = 10000, thin = 10)

## fake a 'gam' object (class 'jam')
jam <- sim2jam(sam, jd$pregam)

par(mfrow = c(1, 3))
plot(jam)  ## the smooth in the constraint null space
hist(sam$b[2, , 1], xlab = expression(beta), main = "")  ## beta

## plot 25 sample curves (unconstrained)....
days <- 152:674
pd <- data.frame(days = days, ozone = days * 0)
Xp <- predict(jam, newdata = pd, type = "lpmatrix")
ii <- 1:25 * 20 + 500

for (i in ii) {
  fv <- Xp %*% sam$b[, i, 1]
  if (i == ii[1]) 
    plot(days, fv, type = "l", ylim = c(4, 7), ylab = "log(size)") else lines(days, fv)
}
