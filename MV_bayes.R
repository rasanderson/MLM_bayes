# Want to try a simple MLM using the vegan mite data

# For simplicity just use SubsDens and WatrCont (cols 1 and 2) as continuous
# Need to turn this into a function, use dplyr for quicker conversion to long
# format, and check that it's OK with factors as explanatories

library(vegan)
library(mcmcplots)
library(data.table)
library(dclone)

rm(list=ls())

extract_linpred <- function (mcmcout) {
  psi.data <- ggs(mcmcout, family="psi")
  linpred <- NULL
  psi.data.dt <- data.table(psi.data[,1], psi.data[,2], psi.data[,3], psi.data[,4])
  setkey(psi.data.dt, Parameter)
  
  for(i in 1:Nobs){
    cond <- paste0("psi[", i, "]")
    sub <- psi.data.dt[cond]$value
    sub <- log(sub/(1-sub))
    linpred[i] <- mean(sub)
  }
  return(linpred)
}



data(mite)
data(mite.env)
mite.env <- mite.env[,1:2]
mite <- decostand(mite, method="pa")
mite.env <- decostand(mite.env, method="standardize") # zero mean unit sd for env

# Traditional method via RDA
mite.rda <- rda(mite ~ SubsDens + WatrCont, data=mite.env)


# Unravel to 'long' format for MLM
mite.lng <- vector("numeric",   length=nrow(mite)*ncol(mite))
SPP      <- vector("character", length=nrow(mite)*ncol(mite))
mite.env.lng = data.frame(array(NA, dim=c(0,ncol(mite.env))))
colnames(mite.env.lng) <- colnames(mite.env)
# OK I know vectorising it or using dplyr is quicker but this is quicker to code
z = 1
for(spp in 1:ncol(mite)){
   for(site in 1:nrow(mite)){
      mite.lng[z] <- mite[site,spp]
      SPP[z] <- colnames(mite)[spp]
      mite.env.lng <- rbind(mite.env.lng, mite.env[site,])
      z <- z+1
   }
}
SPP <- as.factor(SPP)

# Now implement as Bayesian
readline("Hit return to setup Bayesian approach...")
library(rjags)
Nspecies <- ncol(mite)
Nsite <- nrow(mite)
Nobs <- Nspecies*Nsite
Ncov <- ncol(mite.env.lng)

# data:
jags_d <- list(Y=mite.lng,
               X=mite.env.lng,
               Species=SPP,
               Nspecies=Nspecies,
               Ncov=Ncov,
               Nobs=Nobs)

# parameters:
params <- c("alpha", "betas", "sd.beta", "psi")

# Try to speed things up; originall store=1000, nadap=1000; nburn=2000
store <-500
nadap <-500
nburn <-750
thin<-7
chains <- 3


print("First call to JAGS")
cl <- makePSOCKcluster(3)
parJagsModel(cl, name="mod0res", file="jags_models/MLM_0f.txt", data=jags_d, n.chains=3, n.adapt=nadap)
parUpdate(cl, "mod0res", n.iter=nburn)
out0 <- parCodaSamples(cl, "mod0res", params, n.iter=store*thin, thin=thin)
stopCluster(cl)

caterplot(out0, "betas")

print("Second call to JAGS")
cl <- makePSOCKcluster(3)
parJagsModel(cl, name="mod12res", file="jags_models/MLM_12f.txt", data=jags_d, n.chains=3, n.adapt=nadap)
parUpdate(cl, "mod12res", n.iter=nburn)
out_norand <- parCodaSamples(cl, "mod12res", params, n.iter=store*thin, thin=thin)
stopCluster(cl)

caterplot(out_norand, "betas")

################################################
# Check random vs. fixed:
################################################
# If 95 HDI of st.dev. of beta[j] overlaps zero, then fixed effect.
# (i.e. no significant variability in effect among species)
print("About to check random vs fixed")
library(ggmcmc)
source(file="HDI.R")

hdi.sd <- array(0, dim=c(Ncov, 2))

sd.beta.df <- ggs(out0, family="sd.beta") #look at the full model

for(i in 1:Ncov){
   sub <- subset(sd.beta.df, 
                 Parameter==paste("sd.beta[",i,"]", sep="")
   )$value
   hdi.sd[i, ] <- HDI(sub) #HDI of st.dev. for each covariate
}
print(hdi.sd)  # HDI does not encompass zero therefore both are rand effects



################################################
# Extract 'linear predictor' (logit(psi)) of the best model
################################################
print("About to start linear prediction extraction best")
linpred.best <- extract_linpred(out0)

################################################
# Extract 'linear predictor' of the model with NO random effects
################################################
print("About to start linear predictor extraction no random")
linpred.norand <- extract_linpred(out_norand)

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
################################################
# Conduct PCA to ordinate sites, map environmental effects
################################################
print("About to start PCA extraction")
MLM.fitted <- array(linpred.best - linpred.norand, c(Nobs/Nspecies, Nspecies))

rownames(MLM.fitted)=c(1:Nsite)
colnames(MLM.fitted)=paste0("Species", 1:Nspecies)
# standardize over spp

MLM.fitted.standard <- MLM.fitted
for(j in 1:Nspecies){
   MLM.fitted.standard[,j] <- (MLM.fitted[,j]-mean(MLM.fitted[,j]))/sd(MLM.fitted[,j])
} 

ss <- cor(MLM.fitted.standard)
U <- svd(ss)
mlm.fit <- MLM.fitted.standard %*% U$v
mlm.fit <- mlm.fit[,1:2]

# environmental variables (only those with significant random effects)
#envir.vars <- Xcov[, c(1:2)]
envir.vars <- mite.env.lng
mlm.envir <- NULL
for(j in 1:ncol(envir.vars)){
   mlm.envir <- cbind(mlm.envir, envir.vars[,j]*mlm.fit[,1],envir.vars[,j]*mlm.fit[,2])
}

envir.points <- t(array(colMeans(mlm.envir),c(2,dim(mlm.envir)[2]/2)))

# plot mlm
plot(-mlm.fit,xlab="PC1",ylab="PC2",type="n")
text(-mlm.fit,label=c(1:(Nobs/Nspecies)),cex=.75)
#points(-mlm.fit, pch=19, cex=0.5)

arrow.coordMLM <- cbind(array(0,dim(envir.points)),-envir.points)

arrows(arrow.coordMLM[,1],arrow.coordMLM[,2],arrow.coordMLM[,3],arrow.coordMLM[,4], 
       code=2, col="black", length=0.05, lwd=.8)

text(1.3*-envir.points,label=colnames(envir.vars),cex=1, font=2)



