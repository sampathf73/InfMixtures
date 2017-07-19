Alpha <- 1e-15
locus <- "D2S1338"
locnum <- 4
Power <- 3
mypath <- "./"

progsim <- function (fn, Alpha=Alpha, locus=locus, locnum=locnum,
                    Power=Power, mypath=mypath) {
Alpha <- Alpha
locus <- locus
locnum <- locnum
Power <- Power
KKvalue <- 1   # Number of active components KK
filenumber <- fn
niter <- 5000

Resultfile1 <- paste(mypath, "/P", Power, "CS-Loc-", 
                     locnum, "-Alp-", Alpha, "-fn", filenumber,
                     "A.csv", sep = "") 
Resultfile2 <- paste(mypath, "/P", Power, "LL-Loc-", 
                     locnum, "-Alp-", Alpha, "-fn", filenumber,
                     "A.csv", sep = "")
Resultfile3 <- paste(mypath, "/P", Power, "SumLL-Loc-", 
                     locnum, "-Alp-", Alpha, "-fn", filenumber,
                     "A.csv", sep = "")
library(Matrix)
library(mnormt)

#################################################################

setClass("GaussianLR", representation(hh = "list", GLR = "list"),
         contains="list")

GaussianLR <- function(hh) {
  new("GaussianLR",
      # Total number of observations N
      hh$NN <- as.integer(hh$nn), 
      # Cavariates matrix (Nx2)
      hh$XX <- as.matrix(hh$xx), 
      # A column matrix of responses
      hh$YY <- as.matrix(hh$yy), 
      # Shape parameter (a) of Inverse Gamma prior
      hh$IGA <- as.numeric(hh$IGa), 
      # Scale parameter (b) of Inverse Gamma prior
      hh$IGB <- as.numeric(hh$IGb), 
      # Mean vector of betas (2X1)
      hh$MuBeta <- as.matrix(hh$mubeta), 
      # Scale matrix of betas (2X2) 
      hh$ScaleBeta <- as.matrix(hh$ScaleBeta), 
      
      # Total number of observations N
      GLR.NN <- as.integer(hh$NN), 
      # Cavariates matrix (Nx2)
      GLR.XX <- as.matrix(hh$XX),
      # A column matrix of responses (Nx1)
      GLR.YY <- as.matrix(hh$YY), 
      # Shape parameter (a) of Inverse Gamma prior
      GLR.IGA <- as.numeric(hh$IGA),
      # Scale parameter (b) of Inverse Gamma prior
      GLR.IGB <- as.numeric(hh$IGB), 
      # Prior mean vector of betas (2X1) 
      GLR.MuBeta <- as.matrix(hh$MuBeta), 
      # Prior scale matrix of betas (2X2)
      GLR.ScaleBeta <- as.matrix(hh$ScaleBeta),  
      GLR.LL <- as.numeric(0),
      if (hh$NN != 0) {
        stop("Number of observations: nn must be zero")
      } else {
        # Posterior shape parameter (a) of Inverse Gamma 
        GLR.IGATilda <- hh$IGA
        # Posterior scale parameter (b) of Inverse Gamma 
        GLR.IGBTilda <- hh$IGB 
        # Posterior mean vector of betas (2X1)
        GLR.MuBetaTilda <- hh$MuBeta 
        # Posterior scale matrix of betas (2X2)
        GLR.ScaleBetaTilda <- hh$ScaleBeta 
      },
      
      GLR <- as.list(list(NN = GLR.NN, XX = GLR.XX, YY = GLR.YY,
                          IGA = GLR.IGA, IGB = GLR.IGB,
                          MuBeta = GLR.MuBeta, 
                          ScaleBeta = GLR.ScaleBeta,
                          IGATilda = GLR.IGATilda, 
                          IGBTilda = GLR.IGBTilda,
                          MuBetaTilda = GLR.MuBetaTilda, 
                          ScaleBetaTilda = GLR.ScaleBetaTilda,
                          LL = GLR.LL)))
  GLR
} # End of GaussianLR


logpredictive <- function(GLR, XTilda, YTilda) { 
  # Calculates the log of the predictive probability of a datum 
  # ll = logpredictive(X, A, B, Mu, Beta, Xstar)
  # log predictive probability of YTilda given XTilda 
  #   and other data items in the component
  # log p(YTilda|XTilda, (x_1, y_1),...,(x_n, y_n))
  
  # Degrees of freedom of MVSt distribution
  df <- 2 * GLR$IGATilda 
  # Mean of MVSt distribution
  Mean <- XTilda %*% GLR$MuBetaTilda 
  # Scale parameter of MVSt distribution
  # Var_Cov = Scale * df / (df - 2)
  ScaleTilda <- (GLR$IGBTilda / GLR$IGATilda)* 
    (1 + XTilda %*% GLR$ScaleBetaTilda %*% t(XTilda))
  
  ll <- dmt(YTilda, Mean, ScaleTilda, df, log=TRUE)
  ll    
} # End of logpredictive

CLL <- function(GLR) { 
  # Calculates the log likelihood of a cluster
  CLL <- 0
  for (item in 1:nrow(GLR$YY))
    CLL <- CLL + logpredictive(GLR, matrix(c(GLR$XX[item,]),1,2),
                               matrix(c(GLR$YY[item,]),1,1))
  CLL
}  # End of CLL

additem2D <- function(GLR, xx, yy) { 
  # Adds data item to a component 
  # GLR = additem2D(GLR, xx, yy)
  # Ddds datum (xx, yy) into component GLR
  # xx is a 1x2 matrix in the form (1 x)
    
  if (GLR$NN == 0) { # Add data to an empty component
    # Number of observations in the component
    GLR$NN <- nrow(yy) 
    # Observed covariates of the component
    GLR$XX <- matrix(rbind(xx), ncol=2) 
    # Observed responses of the component
    GLR$YY <- matrix(rbind(yy), ncol=1) 
  } else { # Add data to non-empty components
    # Number of observations
    GLR$NN <- GLR$NN + nrow(yy) 
    # Observed covariates of the component
    GLR$XX <- matrix(rbind(GLR$XX, xx), ncol=2) 
    # Observed responses of the component   
    GLR$YY <- matrix(rbind(GLR$YY, yy), ncol=1)   
  }

  # Posterior shape parameter (a) of Inverse Gamma
  GLR$IGATilda <- GLR$IGA + 0.5 * GLR$NN
  # Posterior scale parameter (b) of Inverse Gamma 
  GLR$IGBTilda <- GLR$IGB + 
                  0.5 * as.numeric(t(GLR$YY - 
                  GLR$XX %*% GLR$MuBeta)%*%
                  solve(Diagonal(GLR$NN) 
                  + GLR$XX %*% GLR$ScaleBeta %*% t(GLR$XX))%*%
                  (GLR$YY - GLR$XX %*% GLR$MuBeta))
  # Posterior scale matrix of betas (2X2)
  GLR$ScaleBetaTilda <- solve(solve(GLR$ScaleBeta) 
                        + t(GLR$XX) %*% GLR$XX) 
  # Posterior mean vector of betas (2X1)
  GLR$MuBetaTilda <- GLR$ScaleBetaTilda %*% 
                     (solve(GLR$ScaleBeta) %*%
                     GLR$MuBeta + t(GLR$XX) %*% GLR$YY)
  GLR$LL <- CLL(GLR)
  
  GLR
} # End of additem2D


delitem2D <- function(GLR, xx, yy) { 
  # Deletes a data item from a component
  # GLR = delitem2D(GLR, xx, yy) 
  # Deletes datum (xx, yy) from component GLR
  # xx is a 1x2 matrix in the form (1 x)
  
  if (GLR$NN == 1) { 
    # Coverts to an empty component after detetion
    # Number of observations in the component
    GLR$NN <- 0 
    # NO observed covariates in the component
    GLR$XX <- NA 
    # No observed responses in the component 
    GLR$YY <- NA 
    # Posterior shape parameter (a) of Inverse Gamma
    GLR$IGATilda <- GLR$IGA 
    # Posterior scale parameter (b) of Inverse Gamma 
    GLR$IGBTilda <- GLR$IGB 
    # Posterior mean vector of betas (2X1)
    GLR$MuBetaTilda <- GLR$MuBeta 
    # Posterior scale matrix of betas (2X2)
    GLR$ScaleBetaTilda <- GLR$ScaleBeta  
    GLR$LL <- 0
  } else { # Deletes a data item from a component
    index = -999
    j = 1
    
    while (index < 0) { # Identify the datum to delete
      if ((GLR$XX[j,2] == xx[1,2]) && (GLR$YY[j] == yy[1,1])) {
         index = j # Row number of the datum to be deleted          
      }      
      j = j + 1
    } # End of while
    
    # Number of observations in the component
    GLR$NN <- GLR$NN - 1 
    # Observed covariates of the component
    GLR$XX <- matrix(GLR$XX[-index,], ncol=2) 
    # Observed responses of the component
    GLR$YY <- matrix(GLR$YY[-index], ncol=1) 
    
    # Posterior shape parameter (a) of Inverse Gamma
    GLR$IGATilda <- GLR$IGA + 0.5 * GLR$NN
    # Posterior scale parameter (b) of Inverse Gamma 
    GLR$IGBTilda <- GLR$IGB 
                    + 0.5 * as.numeric(t(GLR$YY - GLR$XX %*% 
                      GLR$MuBeta)%*%
                      solve(Diagonal(GLR$NN) 
                      + GLR$XX %*% GLR$ScaleBeta %*% 
                      t(GLR$XX) )%*% (GLR$YY 
                      - GLR$XX %*% GLR$MuBeta))
    # Posterior scale matrix of betas (2X2)
    GLR$ScaleBetaTilda <- solve(solve(GLR$ScaleBeta) 
                          + t(GLR$XX) %*% GLR$XX) 
    # Posterior mean vector of betas (2X1)
    GLR$MuBetaTilda <- GLR$ScaleBetaTilda %*% 
                       (solve(GLR$ScaleBeta) %*%
                       GLR$MuBeta + t(GLR$XX) %*% GLR$YY)
    GLR$LL <- CLL(GLR)
    
  } # End of else
  
  GLR
} # End of delitem2D


DPMLR_Init <- function(KK,Alpha,GLR0,xx,yy,zz) {
  # Initialize DP mixture model
  # GLR0 empty GLR component with hh prior, 
  # Active mixture components   
  DPM.KK <- KK 
  # Total number of observations
  DPM.NN <- nrow(yy) 
  # Concentration parameter of DP prior
  DPM.Alpha <- Alpha
  # Mixture Components
  DPM.GLR <- vector(mode = "list", length = KK+1) 
  DPM.XX <- xx # Covarites
  DPM.YY <- yy # Responses
  DPM.ZZ <- zz # Initial cluster assignments (between 1 and KK).
  DPM.nn <- matrix(0,1,KK) # KK number of mpty clusters
  DPM <- list(KK=DPM.KK, NN=DPM.NN, Alpha=DPM.Alpha, GLR=DPM.GLR,
              XX=DPM.XX, YY=DPM.YY, ZZ=DPM.ZZ, nn=DPM.nn)
  
  # Initialize mixture components
  # Component KK+1 takes care of all inactive components
  for (kk in 1:(KK+1)) {
    # Generating KK+1 number of empty clusters
    DPM$GLR[[kk]] <- GLR0      
  }
  
  # Add data items into mixture components
  for (ii in 1:DPM$NN) {
    kk = zz[ii] # Identify the cluster index of the datum ii
    DPM$GLR[[kk]] <- additem2D(DPM$GLR[[kk]],matrix(xx[ii,],1,2),
                               matrix(yy[ii],1,1)) 
                    # Add the datum to the cluster  
    DPM$nn[[kk]] <- DPM$nn[[kk]] + 1# Cluster size increase by 1    
  }  
  DPM
  
} # End of DPMLR_Init

temp4 <- matrix(NA,nrow=500,ncol=103)
Resultfile4 <- paste(mypath, "/P", Power, "Gibs-Loc-",
                     locnum, "-Alp-", Alpha, "-fn", filenumber,
                     "A.csv", sep = "")

DPMLR_Gibbs <- function(DPMLR, niter, temp4) { 
  # Gibbs sampler for DPMLR
  KK <- DPMLR$KK # Number of active clusters
  NN <- DPMLR$NN # Total number of data items
  Alpha <- DPMLR$Alpha # Dispersion parameter of DP prior
  GLR <- DPMLR$GLR # A vector of mixture components
  XX <- DPMLR$XX # A 2-column matrix of covariates
  YY <- DPMLR$YY # A column matrix
  ZZ <- DPMLR$ZZ # Cluster indicators
  nn <- DPMLR$nn # Number of data items in each cluster
  
  for (i in 1:niter) { 
    # In each iteration, remove each data item from the model
    # Then add it back according to the conditional probabilities 
    
    for (ii in 1:NN) { # iterate over data items ii
      # Remove data item xx[ii] from component GLR[kk]
      kk <- ZZ[ii] # Current component data item ii belongs to
      # Number of data items in component kk is reduced by 1
      nn[kk] <- nn[kk] - 1 
      GLRtemp1 <- GLR[kk][[1]] # Component kk
      # Remove data item from component kk
      GLRtemp2 <- delitem2D(GLRtemp1, matrix(XX[ii,],1,2), 
                            matrix(YY[ii],1,1)) 
      # Component kk after removing the data item
      GLR[kk][[1]] <- GLRtemp2 
     
      # Delete the component if it has become empty
      if (nn[kk] == 0)
      {
        KK <- KK - 1 # Number of components is reduced by 1
        GLR <- GLR[-kk] # Delete the empty component
        # nn related to empty component is removed
        nn <- nn[-kk] 
        idx <- which(ZZ>kk)
        #  Adjust all the indicators of the components 
        #    after component kk   
        ZZ[idx] <- ZZ[idx] - 1      
      }
      
      # compute conditional probabilities pp(kk) 
      #    of data item ii belonging to each component kk
      # compute probabilities in log domain, then exponentiate
      #logpredictive(N, X, Y, A, B, MuBeta, ScaleBeta, 
      #    XTilda, YTilda)
      pp <- log(c((nn), Alpha))
      for (kkk in 1:(KK+1)) {
        GLRtemp3 <- GLR[kkk][[1]]
        pp[kkk] <- pp[kkk] + Power * logpredictive(GLRtemp3,
                  matrix(XX[ii,],1,2), matrix(YY[ii],1,1))
      }
      pp <- exp(pp - max(pp)) # -max(p) for numerical stability
      pp <- pp / sum(pp)
      
      # Select component kk by sampling from conditional 
      #    probabitilies
      uu <- runif(1)
      kk <- 1+sum(uu>cumsum(pp))
            
      # When a new active component is required
      if (kk == KK+1) {
        # Increse number of components by 1
        KK <- KK + 1 
        # Number of observations in the new component
        nn[kk] <- 0 
        # Increse the indicator of previous empty 
        #     component by 1
        GLR[kk+1] <- GLR[kk] 
      }
      
      # Add data item xx[ii] back into model (component GLR[kk])
      ZZ[ii] <- kk
      # Number of data items in component kk is reduced by 1
      nn[kk] <- nn[kk] + 1 
      GLRtemp3 <- GLR[kk][[1]] # Component kk
      # Add data item to component kk
      GLRtemp4 <- additem2D(GLRtemp3, matrix(XX[ii,],1,2),
                            matrix(YY[ii],1,1)) 
      # Component kk after adding the data item
      GLR[kk][[1]] <- GLRtemp4  
    
    } # End of iteration over data items
    if ((i > 3000)&(i%%4 == 0)) {
      r <- (i - 3000)/4
      temp4[r,1] = r
      temp4[r,2:(KK+1)]=nn
      sumt <- 0
      for (k in 1:KK) {
        temp4[r,(51+k)]=GLR[[k]]$LL
        sumt <- sumt + GLR[[k]]$LL
      }
      temp4[r,102]= KK
      temp4[r,103]= sumt
    }
    
  } # End of iteration
  
  write.csv(temp4, Resultfile4)
  
  # Update DPMLR object
  DPMLR$GLR <- GLR
  DPMLR$ZZ <- ZZ
  DPMLR$nn <- nn
  DPMLR$KK <- KK  
  DPMLR
  
} # End of Gibbs sampler

#################### APPLICATION ###########################

##### NGM DATA #####

df0 <- read.csv("NGMdata.csv",header=TRUE)
df1 <- df0[,-1]
df2 <- df1
library(plyr)
df2 <- rename(df2, c("Marker" = "locus", 
                     "Marker.code" = "marker" ,
                     "Stutter.Height" = "stheight",
                     "Allele" = "allele",
                     "LUS" = "LUS",
                     "Allele.Height" = "height",
                     "SR" = "sr"))
locname1 <- as.character(df2$locus[df2$marker==locnum][1])
df <- df2[df2$locus==locname1,]
sampsize <- nrow(df)
set.seed(as.integer(filenumber))
priordataIndex <- sort(sample(1:sampsize,5,replace=F))
sr1 <- df$sr
LUS1 <- df$LUS
y0 <- sr1[priordataIndex]
x0 <- matrix(c(rep(1,5), LUS1[priordataIndex]), 
             ncol=2, byrow=FALSE)
lm0 <- lm(y0 ~ x0[,2])

Beta0 <- summary(lm0)[['coefficients']][1,1] # Beta0
Beta1 <- summary(lm0)[['coefficients']][2,1] # Beta1
MuBeta0 <- matrix(c(Beta0, Beta1),2,1) # (Beta0, Beta1)^T
sigma0 <- summary(lm0)[['sigma']]
cov.unscaled <- matrix(summary(lm0)[['cov.unscaled']],2,2)
ScaleBeta0 <- cov.unscaled # VBeta
# sigma0 = Residual standard error of the linear model 
#        = summary(lm0)[['sigma']] = sqrt(MSE)
# Mean(Beta) = MuBeta
# Var-Cov(Beta) = VBeta * Sigmasquared
################ With n0, x0, y0 historical data: ############
# Mean(Beta)0 = MuBeta0 
#  = Beta coefficients of linear model 
#       based on the sample of size n0
# Mean(Beta)0 = (X0^TX0)^(-1)X0^TY0
# VBeta0 = cov.unscaled = (X0^TX0)^(-1)
# Var-Cov(Beta)0 = (sigma0)^2 * ScaleBeta0 
#                = (sigma0)^2 * VBeta0 
#                = sigma0^2 (X0^TX0)^(-1) 
# sigmasquare ~ IG(A,B)  
# A = (n0 - k)/2 and B = (n0 - k)* MSE /2

##### Posterior Predictive Distribution #############

# yTilda|y ~ MVSt(2ATilda)[XTilda*MuBetaTilda,
#            (BTilda/ATilda)*(I+XTilda*VBetaTilda*XTilda^T)]

A <- 1.5 # since k=2, n0 = 5
B <- 1.5*sigma0^2

hh0 <- list(nn=0, xx=NA, yy=NA, IGa=NA, IGb=NA, 
            mubeta=NA, ScaleBeta=NA)
hh0$IGa <- A
hh0$IGb <- B
hh0$mubeta <- MuBeta0
hh0$ScaleBeta <- ScaleBeta0
GLR0 <- GaussianLR(hh0)

sr2 <- sr1[-(priordataIndex)]
LUS2 <- LUS1[-(priordataIndex)]

xxall <- matrix(c(rep(1,sampsize-5), LUS2), ncol=2, 
                byrow=FALSE)
yyall <- matrix(c(sr2), ncol=1)

# Initial parameters of Gibbs sampler
KK <- KKvalue # Number of active components
#Alpha  # Dispersion parameter of Dirichlet prior
NN <- nrow(yyall) # Total number of observations N
# Initial component indicaters
zz <-  ceiling(runif(NN)*KK) 

# Initialize DPMLR object
DPMLR = DPMLR_Init(KK,Alpha,GLR0,xxall,yyall,zz)

temp1 <- matrix(NA,nrow=1,ncol=50)
temp2 <- matrix(NA,nrow=1,ncol=50)
temp3 <- matrix(NA,nrow=1,ncol=2)


RDatafile <- paste(mypath, "/P", Power, "Loc-",
                   locnum, "-Alp-", Alpha, "-fn", filenumber, 
                   "A.RData", sep = "")
DPMLRtemp <- DPMLR_Gibbs(DPMLR, niter, temp4)
save(DPMLRtemp, file = RDatafile)
temp1[1,1:length(DPMLRtemp$nn)]=DPMLRtemp$nn
sum <- 0
for (j in 1:length(DPMLRtemp$nn)) {
  temp2[1,j] <- DPMLRtemp$GLR[[j]]$LL
  sum <- sum + DPMLRtemp$GLR[[j]]$LL
}
temp3[1,1] <- length(DPMLRtemp$nn)
temp3[1,2] <- sum

write.csv(temp1, Resultfile1)
write.csv(temp2, Resultfile2)
write.csv(temp3, Resultfile3)

} # End of progsim

###########################################################

RunProg <- function(funnum) {
  if (funnum == 1) progsim (fn=1, Alpha=Alpha, locus=locus, 
                            locnum=locnum, Power=Power,
                            mypath=mypath)
  if (funnum == 2) progsim (fn=2, Alpha=Alpha, locus=locus, 
                            locnum=locnum, Power=Power, 
                            mypath=mypath)
  if (funnum == 3) progsim (fn=3, Alpha=Alpha, locus=locus, 
                            locnum=locnum, Power=Power, 
                            mypath=mypath)
  if (funnum == 4) progsim (fn=4, Alpha=Alpha, locus=locus, 
                            locnum=locnum, Power=Power, 
                            mypath=mypath)
  if (funnum == 5) progsim (fn=5, Alpha=Alpha, locus=locus, 
                            locnum=locnum, Power=Power, 
                            mypath=mypath)
} # End of RunProg

library("doParallel")
library("foreach")
cl <- makeCluster(5)
registerDoParallel(cl)
prognum <- 1:5
foreach(i=1:length(prognum)) %dopar% RunProg(prognum[i])
stopCluster(cl)
csvfile1 <- paste("P", Power, "CS-Loc-", locnum, 
                  "-Alp-", Alpha, "-fn", 1, "A.csv", sep = "") 
csvfile2 <- paste("P", Power, "LL-Loc-", locnum, 
                  "-Alp-", Alpha, "-fn", 1, "A.csv", sep = "")
csvfile3 <- paste("P", Power, "SumLL-Loc-", locnum, 
                  "-Alp-", Alpha, "-fn", 1, "A.csv", sep = "")
mydata1 = read.csv(csvfile1, header=T)
mydata2 = read.csv(csvfile2, header=T)
mydata3 = read.csv(csvfile3, header=T)

for (ff in 2:5) {
  csvfiletemp1 <- paste("P", Power, "CS-Loc-", 
                        locnum, "-Alp-", Alpha, "-fn", ff, 
                        "A.csv", sep = "") 
  csvfiletemp2 <- paste("P", Power, "LL-Loc-", 
                        locnum, "-Alp-", Alpha, "-fn", ff, 
                        "A.csv", sep = "")
  csvfiletemp3 <- paste("P", Power, "SumLL-Loc-", 
                        locnum, "-Alp-", 
                        Alpha, "-fn", ff, "A.csv", sep = "")
  mydatatemp1 = read.csv(csvfiletemp1, header=T)
  mydatatemp2 = read.csv(csvfiletemp2, header=T)
  mydatatemp3 = read.csv(csvfiletemp3, header=T)
  
  mydata1 = rbind(mydata1, mydatatemp1)
  mydata2 = rbind(mydata2, mydatatemp2)
  mydata3 = rbind(mydata3, mydatatemp3)
  
}

Fullcsvfile1 <- paste(mypath, "/P", Power, 
                      "CS-Loc-", locnum, "-Alp-", Alpha, 
                      "-Complete-", "A.csv", sep = "") 
Fullcsvfile2 <- paste(mypath, "/P", Power, 
                      "LL-Loc-", locnum, "-Alp-", Alpha, 
                      "-Complete-", "A.csv", sep = "")
Fullcsvfile3 <- paste(mypath, "/P", Power, 
                      "SumLL-Loc-", locnum, "-Alp-", Alpha, 
                      "-Complete-", "A.csv", sep = "")
write.csv(mydata1, Fullcsvfile1)
write.csv(mydata2, Fullcsvfile2)
write.csv(mydata3, Fullcsvfile3)

#########ACTUAL SIMULATION CODES  ENDS  ###############


