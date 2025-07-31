#=====================================================================
# Miscellaneous
#=====================================================================

# Stop quietly
internalClass$set("private", "stop_quietly", function(msg='Stopped!')
{
	rlang::inform(msg)
	rlang::interrupt()
})

# Get matrix as numeric
internalClass$set("private", "get_NumMat", function(M, rownames=TRUE)
{
	M1 <- data.frame(matrix(as.numeric(M), nrow=nrow(M), ncol=ncol(M)))
	colnames(M1) <- colnames(M)
	if (rownames)
		rownames(M1) <- rownames(M)
	M1
})

# Get processing parameters from PROFILE
internalClass$set("private", "get_procParams", function(profile=NULL)
{
	if (!is.null(profile) && !is.null(profile$preprocess)) {
		procParams$LB <<- profile$preprocess$LB
		procParams$ZFFAC <<- profile$preprocess$ZFFAC
		if (profile$preprocess$ZFFAC==0) procParams$ZEROFILLING <<- FALSE
		if (!is.null(profile$preprocess$MVPZTSP)) {
			procParams$MVPZTSP <<- ifelse( profile$preprocess$MVPZTSP != 0, 1, 0)
		}
		if (!is.null(profile$preprocess$DHZPZRANGE)) {
			procParams$DHZPZRANGE <<- profile$preprocess$DHZPZRANGE
		}
		if (abs(profile$preprocess$PHC0)>0) {
			procParams$OPTPHC0 <<- procParams$OPTPHC1 <<- FALSE
			procParams$PHC0 <<- profile$preprocess$PHC0
			procParams$PHC1 <<- 0
		}
	}
	procParams
})

#=====================================================================
# Some functions about spectra under RAWDIR
#=====================================================================

# Get list of samples under RAWDIR
internalClass$set("private", "get_list_samples", function(DIR=NULL)
{
	if (is.null(DIR)) DIR <- RAWDIR
	LIST <- list.files(path = DIR, pattern = "audita.txt$",
					all.files = FALSE, full.names = TRUE, recursive = TRUE, ignore.case = FALSE, include.dirs = FALSE)
	unique(unlist(lapply(LIST, function(f) { V <- unlist(strsplit(f,'/')); V[c(length(V)-2)] })))

})

# Get matrix of spectra for each samples under DIR
internalClass$set("private", "get_list_spectrum", function(DIR, samples)
{
	M <- NULL
	if (! is.null(DIR) && dir.exists(DIR)) {
		dirs <- list.dirs(DIR, recursive = TRUE, full.names = TRUE)
		for (S in samples) {
			SDIR <- dirs[basename(dirs) == S]
			Sexp <- unlist(lapply(list.files(SDIR), function(d) { if (dir.exists(file.path(SDIR,d))) d }))
			for (id in Sexp) {
				ACQDIR <- file.path(SDIR,id)
				ACQFILE <- file.path(ACQDIR,'acqus')
				if (!file.exists(ACQFILE)) next
				ACQ <- readLines(ACQFILE)
				Pstr <- Rnmr1D:::.bruker.get_param(ACQ,"PULPROG",type="string")
				PULSE <- NULL
				if (length(grep('^zg[0-9]{0,2}$',Pstr))) PULSE <- 'zg'
				if (length(grep('^zgpr',Pstr))) PULSE <- 'zgpr'
				if (length(grep('^noesy',Pstr))) PULSE <- 'noesy'
				if (is.null(PULSE)) next
				M <- rbind(M, c(S, id, PULSE))
			}
		}
	}
	M
})

#=====================================================================
# Some functions about spectral information
#=====================================================================

internalClass$set("private", "get_TSP_width", function(spec)
{
	ppmrange <- c(-0.2,0.2)
	iseq <- getseq(spec,ppmrange)
	n <-  iseq[1] + which(spec$int[iseq]==max(spec$int[iseq])) - 1
	Vs <- iseq[1] + which(spec$int[iseq]>=spec$int[n]/2) - 1
	dv <- spec$ppm[ c( Vs[1], Vs[length(Vs)] ) ]
	TSPWIDTH <- round(sum(abs(dv))*spec$acq$SFO1,3)
	TSPWIDTH
})

# Compute Noise level - based on Bruker algo
internalClass$set("private", "get_Vnoise", function(spec, ppm_noise)
{
	I <- getseq(spec, ppm_noise, mode='interval')
	n1 <- min(I); n2 <- max(I)
	size_m <- n2-n1+1;
	size_half <- round(size_m/2);
	Som <- abs(sum(spec$int[n1:n2]))
	SQ <- sum(spec$int[n1:n2]*spec$int[n1:n2])
	SD=0.0;
	for(k in 1:size_half) {
		i1 = I[1] + size_half + k - 1;
		i2 = I[1] + size_half - k ;
		SD <- SD + (k+1)*( spec$int[i1] - spec$int[i2] );
	}
	SD <- abs(SD);
	Vnoise <- sqrt(( SQ - ( Som*Som + 3*SD*SD/(size_m*size_m-1) )/size_m )/(size_m-1) );
	Vnoise
})

# get the index sequence corresponding to the ppm range
internalClass$set("private", "getseq", function(spec, ppm, mode='seq')
{
	if (mode=='seq') { # seq
		c(which(spec$ppm>=ppm[1])[1]:length(which(spec$ppm<=ppm[2])))
	} else { # interval
		c(which(spec$ppm>=ppm[1])[1],length(which(spec$ppm<=ppm[2])))
	}
})

#=====================================================================
# Some processing functions
#=====================================================================

internalClass$set("private", "zerosSpec", function(spec, ppmranges)
{
	V <- spec$int
	if (! is.null(ppmranges) && "numeric" %in% class(ppmranges))
		ppmranges <- matrix(ppmranges, nrow=1, ncol=length(ppmranges))
	if (! is.null(ppmranges) && nrow(ppmranges)>0) {
	for( k in 1:nrow(ppmranges) )
		V[getseq(spec,ppmranges[k,])] <- 0
	}
	V
})

internalClass$set("private", "smoothSpec", function(spec, ppmrange, WS=30)
{
	# Smooth the PPM range
	iseq <- getseq(spec,ppmrange)
	V <- Smooth(spec$int[iseq], WS)
	n2 <- length(V); n1 <- n2 - WS + 1
	a <- (V[n2]-V[n1])/(n2-n1)
	for (j in n1:n2) V[j] <- a*(j-n1) + V[n1]
	spec$int[iseq] <- V
	spec$int
})

internalClass$set("private", "qnmrbc", function(spec, ppmrange, CSIG=5, NLOOP=5, WS=0)
{
	# external baseline at zero by default
	BLext <- rep(0,length(spec$int))
	if (!is.null(ppmrange)) {
		Ynoise <- spec$int[getseq(spec,PPM_NOISE)]
		specSig <- MASS::fitdistr(Ynoise, "normal")$estimate[2]
		dN <- round(0.00075/spec$dppm)
		iseq <- getseq(spec,ppmrange)
		Y <- spec$int[iseq]
		if (WS>0) Y <- Rnmr1D:::Smooth(Y,WS)
		bc <- rep(0,length(iseq))
		for (l in 1:NLOOP) {
			bci <- Rnmr1D:::C_GlobSeg(Y, dN, CSIG*specSig)
			Y <- Y - bci
			bc <- bc + bci
		}
		BLext[iseq] <- bc
   }
   BLext
})

internalClass$set("private", "getBLexternal", function(spec, blset, porder=1, ppmrange=NULL, WS=30)
{
	# external baseline at zero by default
	BLext <- rep(0,length(spec$int))
	if (blset != 0 ) {
		cmax <- switch(porder, 6, 7, 8)
		lambda <- ifelse( blset>0, cmax-blset, abs(blset) )
		# Estimation of a baseline correction
		if (is.null(ppmrange)) {
			Y <- spec$int
			if (WS>0) Y <- Rnmr1D:::Smooth(Y,WS)
			BLext <- Rnmr1D:::.airPLS(Y, lambda=10^lambda, porder=porder)
		} else {
			iseq <- getseq(spec,ppmrange)
			Y <- spec$int[iseq]
			if (WS>0) Y <- Rnmr1D:::Smooth(Y,WS)
			BLext[iseq] <- Rnmr1D:::.airPLS(Y, lambda=10^lambda, porder=porder)
		}
	}
	BLext
})

# Apply a filter on the spectrum - filter must be in daub8, symlet8, smooth0, ..., smooth3
internalClass$set("private", "filterSpectrum", function(spec, filter)
{
	specInt <- spec$int
	if (filter %in% c('daub8', 'symlet8')) {
		specInt <- filterByWT(spec$int, filter, threshold = 0.5)
	}
	if (filter %in% length(grep('smooth*', filter))) {
		fsavgol <- Rnmr1D::getDeconvParams()$flist[[filter]]
		specInt <- filterSavGol(spec$int, fsavgol$m, fsavgol$nl, fsavgol$nr)
	}
	specInt
})
