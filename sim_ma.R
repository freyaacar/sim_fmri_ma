################################################################################
# Read in dataset
  # A meta-analysis on pain from NeuroVault is used (http://neurovault.org/collections/1425/)
  # This consists of 21 studies, the first 10 are analysed with SPM, the last 11
  # with FSL
################################################################################
  # Preamble
    library(oro.nifti)
    library(lattice)

    nstud <- 21
    n.form <- formatC(c(1:nstud), width=2, flag="0")
    DIM <- c(91,109,91)

    df.stud <- array(data=NA,dim=nstud)
    for (n in 1:nstud) {
    	df.stud[n] <- dim(read.table(paste("21painstudies(NIDM-Results)/pain_",n.form[n],".nidm/DesignMatrix.csv",sep="")))[1] - 1
    }
    n.perstud <- df.stud + 1

  # read in t-maps
    tmaps <- array(data=NA,dim=c(nstud,DIM))
    for (n in 1:nstud) {
    	if (n < 11){
    		tmaps[n,,,]<-readNIfTI(paste("21painstudies(NIDM-Results)/pain_",n.form[n],".nidm/TStatistic.nii.gz",sep=""), verbose = FALSE, warn = -1, reorient = TRUE, call = NULL)
    	}
    	else {
    		tmaps[n,,,]<-readNIfTI(paste("21painstudies(NIDM-Results)/pain_",n.form[n],".nidm/TStatistic_T001.nii.gz",sep=""), verbose = FALSE, warn = -1, reorient = TRUE, call = NULL)
    	}
    }

  # ES-maps
    J <- 1-(3/((4*(df.stud))-1))
    ESmaps <- array(data=NA,dim=c(nstud,DIM))
    for (n in 1:nstud) {
    	ESmaps[n,,,] <- tmaps[n,,,]/sqrt(n.perstud[n])*J[n]
    }
	
  # Construct a mask for every study
    masks <- ifelse(tmaps == 0, 0, 1)
    mask <- ifelse(apply(masks,c(2,3,4),mean) == 1, 1, 0)


################################################################################
# Compute within- and between study variance
################################################################################
  # within-study variance
    varHedge <- function(g,N){
	  value <- (1/N) + (1 - (gamma((N - 2) / 2) / gamma((N - 1) / 2))^2 * (N - 3) / 2) * g^2
	    return(round(value,7))
	}

	wsvar <- array(data=NA,dim=c(nstud,DIM))
	for (n in 1:nstud) {
		wsvar[n,,,] <- varHedge(ESmaps[n,,,],n.perstud[n])
	}

  # between-study variance
	tau <- function(Y,W,k){
	  C <- sum(W)-(sum(W^2)/sum(W))
	  df <- k-1
	  Q <- sum(W*Y^2)-sum(W*Y)^2/sum(W)
	  if(Q < df){
	    T2 <- 0
	  }else{
	    T2 <- (Q-df)/C
	  }
	  return (T2)
	}

	bsvar <- array(data=NA,dim=c(nstud,DIM))
	for (n in 1:nstud) {
		bsvar[n,,,]<-tau(ESmaps[n,,,],(1/wsvar[n,,,]),nstud[n])
	}

################################################################################
# Determine parameters for the simulations
################################################################################
  # within-study variance

  # between-study variance

  # number of peaks


  # number of clusters

  # ditribution of cluster sizes

  # distribution of peak heights


################################################################################
# Simulate meta-analysis set and write out results
################################################################################
  # Simulate

  # t-maps

  # z-maps

  # ES-maps

  # FSL output

  # SPM output

  # cluster and peak locations, size and height

  # peak coordinates



################################################################################
#                               EXTRA                                          #
################################################################################
# Human Connectome Project
  # 1. Is there a between-study variance?
  # 2. How much do these differ from our parameters


