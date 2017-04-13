################################################################################
# Read in dataset
  # A meta-analysis on pain from NeuroVault is used (http://neurovault.org/collections/1425/)
  # This consists of 21 studies, the first 10 are analysed with SPM, the last 11
  # with FSL
################################################################################
  # Preamble
    library(oro.nifti)
    library(lattice)

    kstud <- 21
    n.form <- formatC(c(1:kstud), width=2, flag="0")
    DIM <- c(91,109,91)

    df.stud <- array(data=NA,dim=kstud)
    for (n in 1:kstud) {
    	df.stud[n] <- dim(read.table(paste("21painstudies(NIDM-Results)/pain_",n.form[n],".nidm/DesignMatrix.csv",sep="")))[1] - 1
    }
    n.perstud <- df.stud + 1

  # read in t-maps
    tmaps <- array(data=NA,dim=c(kstud,DIM))
    for (n in 1:kstud) {
    	if (n < 11){
    		tmaps[n,,,]<-readNIfTI(paste("21painstudies(NIDM-Results)/pain_",n.form[n],".nidm/TStatistic.nii.gz",sep=""), verbose = FALSE, warn = -1, reorient = TRUE, call = NULL)
    	}
    	else {
    		tmaps[n,,,]<-readNIfTI(paste("21painstudies(NIDM-Results)/pain_",n.form[n],".nidm/TStatistic_T001.nii.gz",sep=""), verbose = FALSE, warn = -1, reorient = TRUE, call = NULL)
    	}
    }

  # ES-maps
    J <- 1-(3/((4*(df.stud))-1))
    ESmaps <- array(data=NA,dim=c(kstud,DIM))
    for (n in 1:kstud) {
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
	    return(value)
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

	bsvar <- array(data=NA,dim=DIM)
	for (x in 1:DIM[1]) {
		for (y in 1:DIM[2]) {
			for (z in 1:DIM[3]) {
				bsvar[x,y,z] <- tau(ESmaps[,x,y,z],1/wsvar[,x,y,z],kstud)
			}
		}
	}

  # check values
	summary(bsvar[mask==1])
	summary(wsvar[1,,,][mask==1])
	bsvar[50,50,51]
	bsvar[50,50,50]
	wsvar[1,50,50,51]
	wsvar[1,50,50,50]

  # check over voxels
	bdw <- bsvar/apply(wsvar,c(2,3,4),mean)
	writeNIfTI(bdw, filename = paste("bdw",sep=''),gzipped=FALSE)
	summary(wdb[mask==1])
	levelplot(wdb[,,30])

################################################################################
# Fixend and random effects meta-analysis
################################################################################

################################################################################
# Determine parameters for the simulations
################################################################################
  # within-study variance
	for (n in 1:kstud) {
		av.wsvar[n] <- median(wsvar[n,which(mask==1,arr.ind=TRUE)])
	}

  # between-study variance
	av.bsvar <- median(bsvar[which(mask==1,arr.ind=TRUE)])

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


