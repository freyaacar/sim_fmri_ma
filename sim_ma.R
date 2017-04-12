################################################################################
# Read in dataset
  # A pain meta-analysis from NeuroVault is used (http://neurovault.org/collections/1425/)
  # This consists of 21 studies, the first 10 are analysed with SPM, the last 11
  # with FSL
################################################################################
  # Preamble
    nstud <- 21
    n.form <- formatC(c(1:nstud), width=2, flag="0")
    DIM <- c(91,109,91)
    library(oro.nifti)

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


  # tmaps

  # z-maps

  # ES-maps

################################################################################
# Compute within- and between study variance
################################################################################
  # within-study variance

  # between-study variance


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


