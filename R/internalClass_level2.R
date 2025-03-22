#=====================================================================
# High level functions for Quantification
#=====================================================================

internalClass$set("private", "standardQuantification", function(stds_loc, samplename, thresfP=5, deconv=TRUE, verbose=1)
{
	# Filtering
	filter <- 'none'

	# Local baseline order
	obl <- 0

	# Ratio Peak/Noise : Only keep the highest peaks
	ratioPN <- 1800

	# lower limit of the ratio of negative intensities to positive intensities between -0.5 and 4.0 ppm
	thresNeg <- 0.25

	# Retrieve preprocessing parameters from PROFILE and disable MVPZTSP correction
	procPars <- get_procParams(PROFILE)
	procPars$MVPZTSP <- FALSE

	# Define deconvolution parameters for peak fitting
	Opars <- list(ratioPN=ratioPN, oneblk=1, pvoigt=1, oeta=1, oasym=1, asymmax=5000, spcv=0.001, d2cv=0.05,
				lowPeaks=0, distPeaks=0.5, addpeaks=0, sndpass=0, etamin=0.05, sigma_min=0.00025, R2limit=0.995)

	# Display sequence and parameter details if verbose mode is enabled
	if (verbose)  cat("\n================================================\n")
	if (verbose) cat('Sequence:',SEQUENCE,"\n")
	if (verbose) cat('Parameters: ratioPN =', Opars$ratioPN,', asymetric =',Opars$oasym,', lowPeaks =',Opars$lowPeaks, ', oblset =',obl,"\n")
	if (verbose) cat('filter =',paste(filters$main,collapse=','),"\n")

	# Retrieve spectra names for the given samplename
	spectra <- get_list_spectrum(QSDIR, samplename)
	M <- spectra[ spectra[,3]==SEQUENCE, , drop=F]
	spectranames <- unlist(lapply(1:nrow(M), function(x){ paste(M[x,1:2], collapse='-')}))

	# Remove TSP/TMSP/DSS from the standard compounds
	stds_loc <- stds_loc[! stds_loc$Compound %in% calib_cmpd_names, , drop=F]

	# Internal matrices to store all results
	P <- nrow(stds_loc) # Nb standards
	N <- nrow(M)           # N spectra
	fPl <- fR <- Mint <- matrix(rep(0,N*P), nrow=N)
	fK <- 0
	bad_rows <- NULL
	expno_list <- NULL

	# Iterate over each spectrum for processing
	for (k in 1:N) {
	# Read spectrum
		if (verbose) cat("\n-------------------\n")
		expno <- M[k,2]
		expno_list <- c(expno_list, expno)
		if (verbose) cat(samplename,', expno=',expno,', sequence =',SEQUENCE,': ')
		ACQDIR <- file.path(QSDIR,samplename,expno)
		spec <- Rnmr1D::readSpectrum(ACQDIR, procPars, PPM_NOISE, NULL, SCALE_INT, verbose= (verbose>1))
	# Display information if verbose
		if (verbose) cat("Path:", ACQDIR,"\n")
		if (verbose) cat("Sequence:",spec$acq$PULSE,"\n")
		if (verbose) cat("SW =",round(spec$acq$SW,4),", SI =",spec$proc$SI,"\n")
		if (verbose) cat("PW =",round(spec$acq$PULSEWIDTH,4),", NS =",spec$acq$NUMBEROFSCANS,"\n")
	# Check if the phasing is correct
		Y <- spec$int[getseq(spec,c(-0.5,4))]
		if (abs(sum(Y[Y<0]))>thresNeg*sum(Y[Y>0])) {
			if (verbose)  cat("ERROR : Phasing failed\n")
			if (verbose)  cat("\n-----------------\n\n")
			bad_rows <- c(bad_rows, k)
			next
		}
	# Check the TSP/TMSP/DSS width for calibration
		TSPwidth <- get_TSP_width(spec)
		if (verbose) cat("TSP width:", TSPwidth,"Hz\n")
		if (TSPwidth>TSPwidthMax) {
			if (verbose)  cat("ERROR : TSP width too large\n")
			if (verbose)  cat("\n-----------------\n\n")
			bad_rows <- c(bad_rows, k)
			next
		}
	# Process each standard compound
		Peaks <- NULL
		K1 <- spec$acq$SW/spec$proc$SI
		K2 <- spec$acq$PULSEWIDTH/spec$acq$NUMBEROFSCANS
		fK <- fK + K2
		for(i in 1:nrow(stds_loc)) {
			C <- as.vector(stds_loc[i, ])
			cmpd <- C$Compound
			ppmrange <- c(C$PPM1,C$PPM2)
			if (verbose) cat("-------------------\n")
			if (verbose) cat(cmpd,": PPM range = [", ppmrange[1],',',ppmrange[2],"]","\n")
			if (verbose) cat("Max Spec  = ", max(spec$int[getseq(spec,ppmrange)]),"\n")

			# Perform peak deconvolution or direct integration
			if (deconv) {
				modelF <- Rnmr1D::LSDeconv(spec, ppmrange, Opars, filters$main, obl, verbose = (verbose>1))
				if (is.null(modelF) || nrow(modelF$peaks)<1) next
				if (verbose) cat(cmpd,": Nb Peaks =",nrow(modelF$peaks)," -  R2 =",round(modelF$R2,4),"\n")
				Peaks <- rbind(Peaks, modelF$peaks)
				Iref <- sum(modelF$peaks$integral)
				if (verbose) cat("Iref :", Iref,"\n")
			} else {
				iseq <- getseq(spec,ppmrange)
				SUM <- 0
				for (l in 1:(length(iseq)-1)) SUM <- SUM + (spec$int[iseq[l+1]]+spec$int[iseq[l]])
				Iref <- 0.5*SUM*spec$dppm
				if (verbose) cat("Iref :", Iref,"\n")
			}

			# Compute concentration factor
			factor <- K1*Iref*(C$MW/C$NH)
			Mint[k,i] <- Iref
			fR[k,i] <- factor
			fPl[k,i] <- round(factor/C$MC,4)
			if (verbose) cat("fPUL =",round(factor/C$MC,4),"\n")
			if (verbose) cat("PW/NS =",round(K2,4),"\n")
		}
		if (verbose) cat("-------------------\n")
		if (deconv) {
			ratioAN <- Peaks$amp/spec$Noise
			names(ratioAN) <- 'Amp/Noise'
			Peaks <- round(cbind(Peaks, ratioAN),4)
			if (verbose) print(Peaks)
			if (verbose) cat("-------------------\n")
		}

	# Compute coefficient of variation (CV) for quality check
		CV <- sd(fPl[k,1:P])/mean(fPl[k,1:P])
		if (verbose) cat("f_PULCON mean:", mean(fPl[k,1:P]),"\n")
		if (verbose) cat("f_PULCON CV:",round(100*CV,2),"\n")
	# Check if the CV is within the acceptable threshold
		if (100*CV>thresfP) {
			if (verbose)  cat("ERROR : f_PULCON CV is higher than",thresfP,"\n")
			if (verbose)  cat("\n-----------------\n\n")
			bad_rows <- c(bad_rows, k)
			next
		}
	}

	# Remove bad rows if needed
	if (!is.null(bad_rows)) {
		fPl <- fPl[! 1:N %in% bad_rows, ]
		fR <- fR[! 1:N %in% bad_rows, ]
		Mint <- Mint[! 1:N %in% bad_rows, ]
		expno_list <- expno_list[! 1:N %in% bad_rows]
	}
	rownames(fPl) <- rownames(fR) <- rownames(Mint) <- expno_list

	if (verbose) cat("\n-------------------\n")
	if (N>1) {
		V <- apply(Mint, 2, function(v) { if (N>2) { 100*sd(v,na.rm=T)/mean(v,na.rm=T) } else { 100*abs(v[1]-v[2])/mean(v,na.rm=T) } })
		if (verbose) { cat("CV Integral(%) : "); print(round(V,2)) }
	}

	# Compute final statistics
	V <- apply(fPl,1,mean)
	fP_CV <- sd(V)/mean(V)
	fK <- fK/N
	fPUL <- list(mean=mean(V), CV=round(100*fP_CV,2))
	if (verbose) cat("f_PULCON mean:", mean(V),"\n")
	if (verbose) cat("f_PULCON CV:",round(100*fP_CV,2),"\n")
	if (verbose) cat("PW/NS:",round(fK,4),"\n")
	if (verbose) cat("-------------------\n")
	colnames(Mint) <- stds_loc[,2]
	MC <- stds_loc[,5]

	# Return results as an object
	obj <- list(sampletype='QS', fPUL=fPUL, fP=fPl, fR=fR, MC=MC, INTG=Mint, fK=fK)
	class(obj) <- 'QC-QS'
	obj
})

internalClass$set("private", "sampleQuantification", function(samplename, expno, zones=NULL, ncpu=2, verbose=1)
{
	# If no specific quantification zones are provided, use the default ones from PROFILE
	if (is.null(zones)) zones <- unique(PROFILE$quantif$zone)

	# Preprocessing: Read and process the raw spectrum (fid)
	if (verbose) cat(samplename,', expno=',expno,': ')
	ACQDIR <- file.path(RAWDIR,samplename,expno)
	spec <- applyReadSpectrum(ACQDIR)
	if (verbose) cat("Path:", ACQDIR,"\n")
	if (verbose) cat("Sequence:",spec$acq$PULSE,"\n")

	# Store sample name and experiment number in the spectrum object
	spec$samplename <- samplename
	spec$expno <- expno

	# Apply baseline correction to the spectrum
	spec <- applyBLcorrection(spec, verbose=verbose)

	# Calculate TSP width and check if it's within acceptable limits
	spec$TSPwidth <- get_TSP_width(spec)
	if (verbose) cat("TSP width:", spec$TSPwidth,"Hz\n")
	if (spec$TSPwidth>1 && verbose)
		cat("ERROR : TSP width too large\n")

	# Peak fitting: Identify peaks in the spectrum
	opars.loc <- opars
	opars.loc$peaks <- NULL
	t <- system.time({
		spec <- applyPeakFitting(spec, opars=opars.loc, zones=zones, ncpu=ncpu, verbose=verbose)
	})
	if (verbose) cat('elapsed time =', round(t[3],2),', Ended at ',format(Sys.time(), "%m/%d/%Y - %X"),"\n\n")

	# Initialize output variables
	Mquant <- SNR <- peaklist <- NULL

	if (verbose) { print(spec$fit$infos); cat("\n") }

	# If peaks were detected, proceed with quantification
	if (! is.null(spec$fit$peaks)) {
		# Perform quantification based on the fitted peaks
		Q <- applyQuantification(spec, fullPattern=TRUE, verbose=verbose)
		print(Q$quantification); cat("\n")
		print(Q$peaklist); cat("\n")

		# Store quantification results in the spectrum object
		spec$Q <- Q

		# Organize and format quantification results
		Mquant <- matrix(Q$quantification[, 4], ncol=1)
		SNR <- matrix(Q$quantification[, 5], ncol=1)
		rownames(Mquant) <- rownames(SNR) <- Q$peaklist[,1]
		colnames(Mquant) <- colnames(SNR) <- samplename
		if (verbose) { print(Mquant); cat("\n") }
	}
	if (verbose)  cat("\n-----------------\n\n")

	# Return the processed spectrum and quantification results
	list(spec=spec, quantMat=Mquant, SNR=SNR)
})

internalClass$set("private", "absSampleQuantification", function(spec, fP, quantMat, fdil=1, verbose=1)
{
	# Compute the Kx constant : we assume that all spectra are acquired with the same instrument parameters
	acq <- spec$acq    # Acquisition parameters
	proc <- spec$proc  # Processing parameters

	# Calculate constants for quantification
	K1 <- spec$acq$SW/spec$proc$SI
	K2 <- spec$acq$PULSEWIDTH/spec$acq$NUMBEROFSCANS


	# Initialize result matrix
	M <- NULL


	# Reference factor for absolute quantification
	Kref <- fP$mean*fP$fK
	if (verbose) cat("Kref =",Kref,", K1 =",K1,", K2 =",K2,", Fdilution =",fdil,"\n")

	# Compute the absolute quantification for each compound
	for (k in 1:nrow(quantMat)) {
		# Extract compound name
		cmpd <- rownames(quantMat)[k]
		# Get the number of protons (NH) contributing to the signal for the compound
		NH <- sum(PROFILE$quantif[ PROFILE$quantif$compound==cmpd, ]$np)
		# Get the molecular weight (MW) of the compound
		MW <- PROFILE$compound[PROFILE$compound$name==cmpd, ]$mw
		# Compute the absolute quantification (Ix)
		Ix <- (quantMat[k,1]/fdil)*(MW/NH)*K1*K2/Kref
		if (verbose) cat(k,"-",cmpd,": Ix =",Ix,"\n")
		# Append result to the matrix
		M <- rbind(M, Ix)
	}

	# Return the absolute quantification matrix
	rownames(M) <- rownames(quantMat)
	M
})

