###################################################
### chunk number 1: 
###################################################
library("VGAM")
ps.options(pointsize = 12)
options(width = 72, digits = 4)
options(SweaveHooks = list(fig = function() par(las = 1)))
options(prompt = "R> ", continue = "+")


###################################################
### chunk number 2: 
###################################################
pneumo <- transform(pneumo, let=log(exposure.time))
fit <- vgam(cbind(normal, mild, severe) ~ s(let, df=2),
            cumulative(reverse=TRUE, parallel=TRUE), pneumo)


###################################################
### chunk number 3: 
###################################################
journal <- c("Biometrika", "Comm.Statist", "JASA", "JRSS-B")
squaremat <- matrix(c(NA, 33, 320, 284,   730, NA, 813, 276,
                      498, 68, NA, 325,   221, 17, 142, NA), 4, 4)
dimnames(squaremat) <- list(winner = journal, loser = journal)


###################################################
### chunk number 4: 
###################################################
abodat <- data.frame(A = 725, B = 258, AB = 72, O = 1073)
fit <- vglm(cbind(A, B, AB, O) ~ 1, ABO, abodat)
coef(fit, matrix = TRUE)
Coef(fit) # Estimated pA and pB


###################################################
### chunk number 5: 
###################################################
head(wffc.nc, 5)


###################################################
### chunk number 6: 
###################################################
fnc <- transform(wffc.nc,
                 finame = factor(iname),
                 fsector = factor(sector),
                 fday = factor(ceiling(session / 2)),
                 mornaft = 1 - (session %% 2),
                 fbeatboat = factor(beatboat))

fnc <- fnc[with(fnc, !is.element(comid, c(99,72,80,93,45,71,97,78))),] 
fnc <- transform(fnc,
                ordnum = ifelse(numbers <= 02, "few",
                         ifelse(numbers <= 10, "more", "most")))
fnc$ordnum <- ordered(fnc$ordnum, levels = c("few", "more", "most"))


###################################################
### chunk number 7: 
###################################################
with(fnc, table(ordnum))


###################################################
### chunk number 8: 
###################################################
fit.pom <- vglm(ordnum ~
          fsector +
          mornaft +
          fday +
          finame,
          family = cumulative(parallel = TRUE, reverse = TRUE),
          data = fnc)


###################################################
### chunk number 9: 
###################################################
head(fit.pom@y, 3)
colSums(fit.pom@y)


###################################################
### chunk number 10: 
###################################################
head(coef(fit.pom, matrix = TRUE), 10)
#head(summary(fit.pom)@coef3, 10) # Old now since 0.7-10 is nicer


###################################################
### chunk number 11: 
###################################################
head(coef(summary(fit.pom)), 10)


###################################################
### chunk number 12: 
###################################################
fit.ppom <- vglm(ordnum ~
          fsector +
          mornaft +
          fday +
          finame,
          cumulative(parallel = FALSE ~ 1 + mornaft, reverse=TRUE),
          data = fnc)
head(coef(fit.ppom, matrix = TRUE),  8)


###################################################
### chunk number 13: 
###################################################
pchisq(deviance(fit.pom) - deviance(fit.ppom),
       df = df.residual(fit.pom) - df.residual(fit.ppom), lower.tail=FALSE)


###################################################
### chunk number 14: 
###################################################
fit2.ppom <- vglm(ordnum ~
          fsector +
          mornaft +
          fday +
          finame,
          family = cumulative(parallel = FALSE ~ 1 + fday, reverse = TRUE),
          data=fnc)
head(coef(fit2.ppom, matrix = TRUE), 8)


###################################################
### chunk number 15: 
###################################################
head(fitted(fit2.ppom), 3)


###################################################
### chunk number 16: 
###################################################
head(predict(fit2.ppom), 3)


###################################################
### chunk number 17: 
###################################################
dim(model.matrix(fit2.ppom, type = "lm"))
dim(model.matrix(fit2.ppom, type = "vlm"))


###################################################
### chunk number 18: 
###################################################
constraints(fit2.ppom)[c(1, 2, 5, 6)]


###################################################
### chunk number 19: 
###################################################
head(nzmarital, 4)
summary(nzmarital)


###################################################
### chunk number 20: 
###################################################
fit.ms <- vgam(mstatus ~ s(age, df = 3), multinomial(refLevel=2),
               data = nzmarital)


###################################################
### chunk number 21: 
###################################################
head(fit.ms@y, 4)
colSums(fit.ms@y)


###################################################
### chunk number 22: 
###################################################
# Plot output
mycol <- c("red","darkgreen","blue")
 par(mfrow=c(2,2))
plot(fit.ms, se=TRUE, scale=12,
         lcol=mycol, scol=mycol)

# Plot output overlayed
#par(mfrow=c(1,1))
plot(fit.ms, se=TRUE, scale=12,
         overlay=TRUE,
         llwd=2,
         lcol=mycol, scol=mycol)


###################################################
### chunk number 23: 
###################################################
# Plot output
mycol <- c("red","darkgreen","blue")
 par(mfrow=c(2,2))
 par(mar=c(4.2,4.0,1.2,2.2)+0.1)
plot(fit.ms, se=TRUE, scale=12,
         lcol=mycol, scol=mycol)

# Plot output overlaid
#par(mfrow=c(1,1))
plot(fit.ms, se=TRUE, scale=12,
         overlay=TRUE,
         llwd=2,
         lcol=mycol, scol=mycol)


###################################################
### chunk number 24: 
###################################################
plot(fit.ms, deriv=1, lcol=mycol, scale=0.3)


###################################################
### chunk number 25: 
###################################################
# Plot output
 par(mfrow=c(1,3))
 par(mar=c(4.5,4.0,0.2,2.2)+0.1)
plot(fit.ms, deriv=1, lcol=mycol, scale=0.3)


###################################################
### chunk number 26: 
###################################################
foo <- function(x, elbow=50)
    poly(pmin(x, elbow), 2)

clist <- list("(Intercept)" = diag(3),
             "poly(age, 2)" = rbind(1, 0, 0),
             "foo(age)" = rbind(0, 1, 0),
             "age" = rbind(0, 0, 1))
fit2.ms <-
    vglm(mstatus ~ poly(age, 2) + foo(age) + age,
         family = multinomial(refLevel=2),
         constraints=clist,
         data=nzmarital)


###################################################
### chunk number 27: 
###################################################
coef(fit2.ms, matrix = TRUE)


###################################################
### chunk number 28: 
###################################################
par(mfrow=c(2,2))
plotvgam(fit2.ms, se=TRUE, scale=12,
         lcol=mycol[1], scol=mycol[1], which.term=1)
plotvgam(fit2.ms, se=TRUE, scale=12,
         lcol=mycol[2], scol=mycol[2], which.term=2)
plotvgam(fit2.ms, se=TRUE, scale=12,
         lcol=mycol[3], scol=mycol[3], which.term=3)


###################################################
### chunk number 29: 
###################################################
# Plot output
par(mfrow=c(2,2))
 par(mar=c(4.5,4.0,1.2,2.2)+0.1)
plotvgam(fit2.ms, se=TRUE, scale=12,
         lcol=mycol[1], scol=mycol[1], which.term=1)
plotvgam(fit2.ms, se=TRUE, scale=12,
         lcol=mycol[2], scol=mycol[2], which.term=2)
plotvgam(fit2.ms, se=TRUE, scale=12,
         lcol=mycol[3], scol=mycol[3], which.term=3)


###################################################
### chunk number 30: 
###################################################
deviance(fit.ms) - deviance(fit2.ms)


###################################################
### chunk number 31: 
###################################################
(dfdiff <- df.residual(fit2.ms) - df.residual(fit.ms))


###################################################
### chunk number 32: 
###################################################
1-pchisq(deviance(fit.ms) - deviance(fit2.ms), df=dfdiff)


###################################################
### chunk number 33: 
###################################################
ooo <- with(nzmarital, order(age))
with(nzmarital, matplot(age[ooo], fitted(fit.ms)[ooo,],
     type="l", las=1, lwd=2, ylim=0:1,
     ylab="Fitted probabilities",
     xlab="Age", # main="Marital status amongst NZ Male Europeans",
     col=c(mycol[1], "black", mycol[-1])))
legend(x=52.5, y=0.62, # x="topright",
       col=c(mycol[1], "black", mycol[-1]),
       lty=1:4,
       legend=colnames(fit.ms@y), lwd=2)
abline(v=seq(10,90,by=5), h=seq(0,1,by=0.1), col="gray", lty="dashed")


###################################################
### chunk number 34: 
###################################################
 par(mfrow=c(1,1))
 par(mar=c(4.5,4.0,0.2,0.2)+0.1)
ooo <- with(nzmarital, order(age))
with(nzmarital, matplot(age[ooo], fitted(fit.ms)[ooo,],
     type="l", las=1, lwd=2, ylim=0:1,
     ylab="Fitted probabilities",
     xlab="Age",
     col=c(mycol[1], "black", mycol[-1])))
legend(x=52.5, y=0.62,
       col=c(mycol[1], "black", mycol[-1]),
       lty=1:4,
       legend=colnames(fit.ms@y), lwd=2.1)
abline(v=seq(10,90,by=5), h=seq(0,1,by=0.1), col="gray", lty="dashed")


###################################################
### chunk number 35: 
###################################################
# Scale the variables? Yes; the Anderson (1984) paper did (see his Table 6).
head(backPain, 4)
summary(backPain)
backPain <- transform(backPain, sx1 = -scale(x1), sx2 = -scale(x2), sx3 = -scale(x3))


###################################################
### chunk number 36: 
###################################################
bp.rrmlm1 <- rrvglm(pain ~ sx1 + sx2 + sx3, multinomial, backPain)


###################################################
### chunk number 37: 
###################################################
Coef(bp.rrmlm1)


###################################################
### chunk number 38: 
###################################################
set.seed(123)


###################################################
### chunk number 39: 
###################################################
bp.rrmlm2 <- rrvglm(pain ~ sx1 + sx2 + sx3, multinomial, backPain, Rank=2,
                   Corner=FALSE, Uncor=TRUE)


###################################################
### chunk number 40: 
###################################################
biplot(bp.rrmlm2, Acol="blue", Ccol="darkgreen", scores=TRUE,
#      xlim=c(-1,6), ylim=c(-1.2,4), # Use this if not scaled
       xlim=c(-4.5,2.2), ylim=c(-2.2, 2.2), # Use this if scaled
       chull=TRUE, clty=2, ccol="blue")


###################################################
### chunk number 41: 
###################################################
# Plot output
 par(mfrow=c(1,1))
 par(mar=c(4.5,4.0,0.2,2.2)+0.1)

biplot(bp.rrmlm2, Acol="blue", Ccol="darkgreen", scores=TRUE,
#      xlim=c(-1,6), ylim=c(-1.2,4),  # Use this if not scaled
       xlim=c(-4.5,2.2), ylim=c(-2.2, 2.2),  # Use this if scaled
       chull=TRUE, clty=2, ccol="blue")


###################################################
### chunk number 42: 
###################################################
iam(NA, NA, M = 4, both = TRUE, diag = TRUE)


