library(oro.nifti)
fslpath <- '/usr/local/fsl/bin/'

################################################################################
# Determine parameters of real meta-analysis
################################################################################
	kstud <- 21

	n.perstud <- array(data=NA,dim=kstud)
	for (n in 1:kstud) {
		n.perstud[n] <- dim(read.table(paste("21painstudies(NIDM-Results)/pain_",n.form[n],".nidm/DesignMatrix.csv",sep="")))[1]
	}

	n.tot.hcp <- 188
	k.hcp <- 12
	n.perst.hcp <- c(9,9,12,12,12,12,14,16,16,20,25,31)
	ind.rand <- sample(c(1:188))

################################################################################
# Read in copes and varcopes
################################################################################
	# Read in copes and varcopes
		copes <- readNIfTI("HCP/MOTOR_10_cope.nii.gz", verbose = FALSE, warn = -1, reorient = TRUE, call = NULL)
		varcopes <- readNIfTI("HCP/MOTOR_10_varcope.nii.gz", verbose = FALSE, warn = -1, reorient = TRUE, call = NULL)

	# construct mask
		masks <- ifelse(copes*varcopes == 0, 0, 1)
    	mask.hcp <- ifelse(apply(masks,c(1,2,3),mean) == 1, 1, 0)
    	table(mask.hcp)
    	setwd(home)
    	masknifti <- nifti(img=mask.hcp, dim=DIM, datatype=2)
		writeNIfTI(masknifti,filename = "HCP/mask",gzipped=FALSE)

	# Write out cope files per study
		i.temp <- 1
		for (k in 1:k.hcp) {
			cope.temp <- array(data=NA,dim=c(n.perst.hpc[k],DIM))
			varcope.temp <- array(data=NA,dim=c(n.perst.hpc[k],DIM))
			cope.temp <- copes[,,,i.temp:(i.temp+n.perst.hpc[k]-1)]
			varcope.temp <- varcopes[,,,i.temp:(i.temp+n.perst.hpc[k]-1)]
			cope.temp.n <- nifti(img=cope.temp, dim=c(DIM,n.perst.hpc[k]), datatype=16)
			varcope.temp.n <- nifti(img=varcope.temp, dim=c(DIM,n.perst.hpc[k]), datatype=16)
			writeNIfTI(cope.temp.n, filename = paste("HCP/cope",k,sep=''),gzipped=FALSE)
			writeNIfTI(varcope.temp.n, filename = paste("HCP/varcope",k,sep=''),gzipped=FALSE)
			system(paste(fslpath,"fslcpgeom HCP/MOTOR_10_cope HCP/cope",k,sep=""))
			system(paste(fslpath,"fslcpgeom HCP/MOTOR_10_cope HCP/varcope",k,sep=""))
			i.temp <- i.temp + n.perst.hpc[k]
		}

################################################################################
# Compute group level t-maps
################################################################################
	setwd(paste(home,"/HCP/",sep=""))
	for (k in 1:k.hcp) {
		fsl_model_cmd = paste(fslpath,"randomise -i cope",k," -o output",k," -1 --glm_output",sep="")
		system(fsl_model_cmd)
	}

	setwd(home)
	tmaps.hcp<-array(data=NA,dim=c(k.hcp,DIM))
	for (k in 1:k.hcp) {
    	tmaps.hcp[k,,,]<-readNIfTI(paste("HCP/output",k,"_tstat1",sep=""), verbose = FALSE, warn = -1, reorient = TRUE, call = NULL)
    }	

################################################################################
# Read in t and compute ES
################################################################################  
  	# ES-maps
	    J <- 1-(3/((4*(n.perst.hcp-1))-1))
	    ESmaps.hcp <- array(data=NA,dim=c(k.hcp,DIM))
	    for (k in 1:k.hcp) {
	    	ESmaps.hcp[k,,,] <- tmaps.hcp[k,,,]/sqrt(n.perst.hpc[k])*J[k]
	    }
	
  # Construct a mask for every study
    masks <- ifelse(tmaps == 0, 0, 1)
    mask.t.hcp <- ifelse(apply(masks,c(2,3,4),mean) == 1, 1, 0)


################################################################################
# Compute within- and between study variance
################################################################################
  # within-study variance
    varHedge <- function(g,N){
	  value <- (1/N) + (1 - (gamma((N - 2) / 2) / gamma((N - 1) / 2))^2 * (N - 3) / 2) * g^2
	    return(value)
	}

	wsvar.hcp <- array(data=NA,dim=c(k.hcp,DIM))
	for (k in 1:k.hcp) {
		wsvar.hcp[k,,,] <- varHedge(ESmaps.hcp[k,,,],n.perst.hcp[k])
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

	bsvar.hcp <- array(data=NA,dim=DIM)
	for (x in 1:DIM[1]) {
		for (y in 1:DIM[2]) {
			for (z in 1:DIM[3]) {
				bsvar.hcp[x,y,z] <- tau(ESmaps.hcp[,x,y,z],1/wsvar.hcp[,x,y,z],k.hcp)
			}
		}
	}

  # check values
	summary(bsvar.hcp[mask.hcp==1])
	summary(wsvar.hcp[1,,,][mask.hcp==1])
	bsvar[50,50,51]
	bsvar[50,50,50]
	wsvar[1,50,50,51]
	wsvar[1,50,50,50]

################################################################################
# Fixed and random effects meta-analysis
################################################################################











################################################################################
# Perform 2nd level group analysis - Credit to Han! - Later when we have time
################################################################################
	DataWrite <- "HCP/"
	

	for (k in 1:k.hcp) {
		# Write auxiliarly files to DataWrite. We need:
	    # GRCOPE in nifti
	    # GRVARCOPE in nifti
	    # 4D mask
	    # design.mat file
	    # design.grp file
	    # design.con file

	    #----- 1 ----#
	    ### Design.mat
	    fileCon <- paste(DataWrite,"/design.mat",sep="")
	    # Text to be written to the file
	    cat('/NumWaves\t1
	    /NumPoints\t',paste(n.perst.hpc[k],sep=''),'
	    /PPheights\t\t1.000000e+00
	    /Matrix
	    ',rep("1.000000e+00\n",n.perst.hpc[k]),file=fileCon)

	    #----- 2 ----#
	    ### Design.con
	    fileCon <- file(paste(DataWrite,"/design.con", sep=""))
	    	writeLines('/ContrastName1	Group Average
	    /NumWaves	1
	    /NumContrasts	1
	    /PPheights		1.000000e+00
	    /RequiredEffect		5.034
	    /Matrix
	    1.000000e+00
	    ',fileCon)
	    close(fileCon)

	      #----- 3 ----#
	      ### Design.grp
	    fileCon <- paste(DataWrite,"/design.grp",sep="")
	    # Text to be written to the file
	    cat('/NumWaves\t1
	    /NumPoints\t',paste(n.perst.hpc[k],sep=''),'
	    /Matrix
	    ',rep("1\n",n.perst.hpc[k]),file=fileCon)

	    # FSL TIME!
	    setwd(DataWrite)
	    command <- paste(fslpath, 'flameo --cope=cope',k,' --vc=varcope',k,' --mask=mask --ld=study',k,'_stats --dm=design.mat --cs=design.grp --tc=design.con --runmode=flame1', sep='')
	    Sys.setenv(FSLOUTPUTTYPE="NIFTI")
	    system(command)
	}
	setwd(home)


