# ========================================================
# R Journal
# 
# RPESE: Risk and Performance Estimators Standard Errors 
#        with Serially Dependent Data
# 
# Anthony-A. Christidis, R. Doug Martin
# 
# ========================================================

# ------------------------
# Monte Carlo Simulation
# ------------------------

# Required libraries
library(RPESE)
library(RPEIF) 
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
    SRout <- lapply(sim.ts, function(x) SR.SE(x, se.method=c('IFcor', 'IFcorPW')))
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




