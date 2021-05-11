# ========================================================
# R Journal
# 
# RPESE: Risk and Performance Estimators Standard Errors 
#        with Serially Dependent Data
# 
# Anthony-A. Christidis, R. Doug Martin
# 
# ========================================================

# Required libraries
library(RPESE)
library(RPEIF)

# ______________
# CODE EXAMPLES
# ______________

# -----------------------------------------------------------------
# Section: Influence Functions for Risk and Performance Estimators
# -----------------------------------------------------------------

# Loading the hedge fund data
install.packages("PerformanceAnalytics")
data(edhec, package = "PerformanceAnalytics")
class(edhec)
colnames(edhec) <- c("CA", "CTAG", "DIS", "EM", "EMN", "ED", "FIA", "GM", "LS", "MA",
                     "RV", "SS", "FoF")

# Package functions
library(RPEIF)
ls("package:RPEIF")

# Nuisance parameters function arguments
args(nuisParsFn)

# IF shape plot - alternative parameters
par(mfrow = c(2, 1))
outSD <- IF.SD(evalShape = T, IFplot = T, nuisPars = nuisParsFn(mu = 0.02, sd = 0.15))
outSR <- IF.SR(evalShape = T, IFplot = T, nuisPars = nuisParsFn(mu = 0.02, sd = 0.15))

# IF shape plot - parameters from returns series
par(mfrow = c(2, 1))
outSD <- IF.SD(returns = edhec$CA, evalShape = T, IFplot = T)
outSR <- IF.SR(returns = edhec$CA, evalShape = T, IFplot = T)

# IF plot for returns series
SDiftr <- IF.SD(returns = edhec$CA)
SRiftr <- IF.SR(returns = edhec$CA)
par(mfrow = c(3, 1))
plot(edhec$CA, lwd = 0.8, ylab = "Returns", main = "CA Hedge Fund Returns")
plot(SDiftr, lwd = 0.8, main = "IF.SD Transformed Returns")
plot(SRiftr, lwd = 0.8, main = "IF.SR Transformed Returns")

# IF plot for return series with outlier cleaning
iftrFIA <- IF.Mean(returns = edhec$FIA)
iftrFIAclean <- IF.Mean(returns = edhec$FIA, cleanOutliers = T, eff = 0.99)
par(mfrow = c(2, 1))
plot(iftrFIA, main = "FIA IF Transformed Returns", lwd = 0.8)
plot(iftrFIAclean, main = "Outlier Cleaned FIA IF Transformed Returns", lwd = 0.8)

# IF plot for return series with prewhetening
iftrFIAcl <- IF.Mean(returns = edhec$FIA, cleanOutliers = T)
PWiftrFIAcl <- IF.Mean(returns = edhec$FIA, cleanOutliers = T, prewhiten = T)
par(mfrow = c(2, 1))
plot(iftrFIAcl, main = "FIA Outlier Cleaned Returns", lwd = 0.8)
plot(PWiftrFIAcl, main = "Prewhitened FIA Outlier Cleaned Returns", lwd = 0.8)

# ------------------------------------------------------------------
# Section: Application: Hedge Funds Data Standard Errors with RPESE
# ------------------------------------------------------------------

# Package functions
library(RPESE)
ls("package:RPESE")

# Example function arguments
args(SD.SE)
args(SR.SE)

# SE computation with default arguments
SRout <- SR.SE(edhec)

# Printing the output
printSE(SRout)

# SE computation with methods specification
SRout <- SR.SE(edhec, se.method = c("IFiid","BOOTiid","IFcor","IFcorPW","BOOTcor"))
printSE(SRout)

# SE computatoin with outlier cleaning
SRout <- SR.SE(edhec, se.method = "IFcorPW", cleanOutliers = F)
SRoutClean <- SR.SE(edhec, se.method = "IFcorPW", cleanOutliers = T)
clean.compare <- data.frame(SRout$IFcorPW$se, SRoutClean$IFcorPW$se)
names(clean.compare) <- c("With Outliers", "Outliers Cleaned")
row.names(clean.compare) <- names(edhec)
round(clean.compare, 3)

# SE computation with correlations output
SRretCor <- SR.SE(edhec, corOut = c("retCor", "retIFCor"))
printSE(SRretCor)

# SE computation with Exponential-Gamma comparison
Clean.out <- SR.SE(edhec, se.method=c("IFiid","IFcor","IFcorPW"), cleanOutliers = T)
GammaClean.out <- SR.SE(edhec, se.method=c("IFiid","IFcor","IFcorPW"), cleanOutliers = T,      
                        fitting.method = "Gamma")
GammaExp.comparison <- cbind(printSE(Clean.out)[,4], printSE(GammaClean.out)[,4])
colnames(GammaExp.comparison) <- c("IFcorPW", "IFcorPW-Gamma")
rownames(GammaExp.comparison) <- names(edhec)
GammaExp.comparison

# SE computation with frequency comparison
SRall <- SR.SE(edhec, cleanOutliers = T, freq.include = "All")
SRdecimate <- SR.SE(edhec, cleanOutliers = T, freq.include = "Decimate",
                    freq.par = 0.5)
SRtruncate <- SR.SE(edhec, cleanOutliers = T, freq.include = "Truncate",
                    freq.par = 0.5)
frequency.comparison <- cbind(printSE(SRall)[,3], printSE(SRdecimate)[,3],
                              printSE(SRtruncate)[,3])
colnames(frequency.comparison) <- c("IFcorPW-All", "IFcorPW-Decimate",
                                    "IFcorPW-Truncate")
rownames(frequency.comparison) <- names(edhec)
frequency.comparison


# ________________________
# MONTE CARLO SIMULATION
# ________________________

library(nse)
library(metRology)

# Setting seed 
set.seed(0)

# Simulation parameters
N <- 5000
mean <- 0.01
scale <- 0.15
df <- 5
n <- c(60, 120, 240)
phi <- seq(0, 0.5, by=0.1)

# Matrix to store final results
final.results <- matrix(nrow=0, ncol=6)

# Normal distribution case
for(n.id in n){
  
  # Matrix to store results for phi.id
  phi.results <- matrix(nrow=6, ncol=6)
  
  for(phi.id in phi){
    
    # TS Simulation
    sim.ts <- lapply(1:N, function(t, n.id=n.id, phi.id=phi.id) arima.sim(n=n.id, model=list(order=c(1,0,0), ar=phi.id), 
                                                              innov=rnorm(n.id, mean=mean, sd=scale)), # Normal distribution innovations
                                                              # innov=rt.scaled(n.id, df=df, mean=mean, sd=scale)), # t-distribution innovations
                     n.id=n.id, phi.id=phi.id)
    
    # SE Computation
    SRout <- lapply(sim.ts, function(x) SR.SE(x, se.method=c("IFcor", "IFcorPW")))
    SRoutPoint <- sapply(SRout, function(x) x$SR, simplify=TRUE)
    SE.true <- sd(sapply(SRout, function(x) x$SR, simplify=TRUE))
    SE.IFcor <- sapply(SRout, function(x) x$IFcor$se, simplify=TRUE)
    SE.IFcorPW <- sapply(SRout, function(x) x$IFcorPW$se, simplify=TRUE)
    IFdata <- lapply(sim.ts, function(x) IF.SR(as.vector(x)))
    SE.Andrews <- unlist(lapply(IFdata, function(x) nse.andrews(x=as.vector(x), lag.prewhite=0)))
    SE.AndrewsPW <-unlist(lapply(IFdata, function(x) nse.andrews(x=as.vector(x), lag.prewhite=1)))
    SE.NeweyWest <- unlist(lapply(IFdata, function(x) nse.nw(x=as.vector(x), lag.prewhite=0)))
    SE.NeweyWestPW <- unlist(lapply(IFdata, function(x) nse.nw(x=as.vector(x), lag.prewhite=1)))
    
    # CI Computation
    CI.IFcor <- 1 - mean((mean(SRoutPoint) >= SRoutPoint - qt(0.975, n.id-1)*SE.IFcor)*
                           (mean(SRoutPoint) <= SRoutPoint + qt(0.975, n.id-1)*SE.IFcor))
    SE.IFcorPW <- sapply(SRout, function(x) x$IFcorPW$se, simplify=TRUE)
    CI.IFcorPW <- 1 - mean((mean(SRoutPoint) >= SRoutPoint - qt(0.975, n.id-1)*SE.IFcorPW)*
                             (mean(SRoutPoint) <= SRoutPoint + qt(0.975, n.id-1)*SE.IFcorPW))
    CI.Andrews <- 1 - mean((mean(SRoutPoint) >= SRoutPoint - qt(0.975, n.id-1)*SE.Andrews)*
                             (mean(SRoutPoint) <= SRoutPoint + qt(0.975, n.id-1)*SE.Andrews))
    CI.AndrewsPW <- 1 - mean((mean(SRoutPoint) >= SRoutPoint - qt(0.975, n.id-1)*SE.AndrewsPW)*
                               (mean(SRoutPoint) <= SRoutPoint + qt(0.975, n.id-1)*SE.AndrewsPW))
    CI.NeweyWest <- 1 - mean((mean(SRoutPoint) >= SRoutPoint - qt(0.975, n.id-1)*SE.NeweyWest)*
                               (mean(SRoutPoint) <= SRoutPoint + qt(0.975, n.id-1)*SE.NeweyWest))
    CI.NeweyWestPW <- 1 - mean((mean(SRoutPoint) >= SRoutPoint - qt(0.975, n.id-1)*SE.NeweyWestPW)*
                                 (mean(SRoutPoint) <= SRoutPoint + qt(0.975, n.id-1)*SE.NeweyWestPW))

    # Results for fixed phi
    phi.results[,phi.id==phi] <- c(CI.IFcor, CI.IFcorPW, CI.Andrews, CI.AndrewsPW, CI.NeweyWest, CI.NeweyWestPW)
    
  }
  
  # Aggregating results
  final.results <- rbind(final.results, phi.results)
}




