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
		procParams$TSP <<-  ifelse(profile$preprocess$TSP==1, TRUE, FALSE)
		if (profile$preprocess$ZFFAC==0) procParams$ZEROFILLING <<- FALSE
		if (!is.null(profile$preprocess$MVPZTSP)) {
			procParams$MVPZTSP <<- ifelse( profile$preprocess$MVPZTSP != 0, 1, 0)
		}
		if (!is.null(profile$preprocess$DHZPZRANGE)) {
			procParams$DHZPZRANGE <<- profile$preprocess$DHZPZRANGE
		}
		if (!is.null(profile$preprocess$ADDPARAMS)) {
			for(s in unlist(profile$preprocess$ADDPARAMS)) {
				V <- unlist(strsplit(s,'='))
				procParams[V[1]] <<- ifelse( is.na(as.numeric(V[2])), V[2] , as.numeric(V[2]) )
			}
		}
	}
	procParams
})

#=====================================================================
# Some functions about spectra under RAWDIR
#=====================================================================

# Get list of directories under DIR
internalClass$set("private", "get_list_dirs", function(DIR=NULL)
{
	LIST <- NULL
	pfile <- ifelse( procParams$VENDOR=='bruker', "audita.txt$", "procpar$")
	if (DIR != RAWDIR || is.null(RAWDIR_SLIST)) {
		LIST <- unique(dirname(list.files(path = DIR, pattern = pfile, all.files = FALSE, full.names = TRUE, recursive = TRUE, ignore.case = FALSE, include.dirs = FALSE)))
	} else {
		LIST <- RAWDIR_SLIST
	}
	if (DIR == RAWDIR && is.null(RAWDIR_SLIST))
		RAWDIR_SLIST <<- LIST
	LIST
})

# Get list of samples under RAWDIR
internalClass$set("private", "get_list_samples", function(DIR=NULL)
{
	LIST <- get_list_dirs(DIR)
	if (procParams$VENDOR=='bruker')
		LIST <- unique(basename(dirname(LIST)))
	else
		LIST <- unique(basename(LIST))
	LIST
})

# Get matrix of spectra for each samples under DIR
internalClass$set("private", "get_list_spectrum", function(DIR, samples, sequence=NULL)
{
	bruker_spectra_list <- function(S)
	{
		L <- grep(paste0('/',S,'/'), dirs, value = TRUE)
		Sexp <- unique(basename(L), dirs, value = TRUE)
		if (sum(is.na(as.numeric(Sexp)))>0) next # Non-Bruker Spectra
		SDIR <- unique(dirname(L))[1]
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
			if (!is.null(sequence) && sequence != PULSE) next
			cnt <<- cnt + 1
			M[cnt, 1:5] <<- c(S, paste(S,id, sep="-"),id, PULSE, ACQDIR)
			if (cnt==nrow(M)) break
		}
	}

	varian_spectra_list <- function(S)
	{
		L <- grep(paste0('/',S), dirs, value = TRUE)
		ACQDIR <- unique(L)[1]
		ACQFILE <- file.path(ACQDIR,'procpar')
		if (file.exists(ACQFILE)) {
			ACQ <- readLines(ACQFILE)
			PULSE <- Rnmr1D:::.varian.get_param(ACQ,"pslabel",type="string")
			if (!is.null(PULSE) && (is.null(sequence) || (!is.null(sequence) && sequence == PULSE)))  {
				cnt <<- cnt + 1
				M[cnt, 1:5] <<- c(S, paste(S,1, sep="-"),1, PULSE, ACQDIR)
			}
		}
	}

	M <- NULL
	if (! is.null(DIR) && dir.exists(DIR))
	{
		dirs <- get_list_dirs(DIR)
		M <- matrix(0, nrow=length(dirs), ncol=5)
		cnt <- 0
		for (S in samples) {
			if (procParams$VENDOR=='bruker')
				bruker_spectra_list(S)
			else
				varian_spectra_list(S)
			if (cnt==nrow(M)) break
		}
		if (cnt>0) { M <- M[1:cnt, , drop=F] }
		else       { M <- NULL }
		if (!is.null(M)) {
			colnames(M) <- c('Spectrum', 'Samplename', 'expno', 'sequence', 'path')
			M <- as.data.frame(M)
		}
	}
	M
})

# Get the samplecode from SAMPLES based on the ID in the matrix of spectra
internalClass$set("private", "get_samplecode_by_ID", function(Slist, ID)
{
	if (procParams$VENDOR=='bruker')
		samplecode <- SAMPLES[which(paste(SAMPLES[,1], SAMPLES[,3],sep="-") == Slist[ID,2]), 2]
	else
		samplecode <- SAMPLES[SAMPLES[,1] == Slist[ID,1], 2]
	samplecode
})

# Filter the matrix of spectra based on SAMPLES
internalClass$set("private", "get_spectralist_by_SAMPLES", function(Slist)
{
	Slist <- Slist[ Slist[,1] %in% SAMPLES[,1], , drop=F]
	if (procParams$VENDOR=='bruker')
		Slist <- Slist[Slist[,2] %in% paste(SAMPLES[,1], SAMPLES[,3],sep="-"), , drop=F]
	Slist
})

#=====================================================================
# Some functions about spectral information
#=====================================================================

internalClass$set("private", "get_negRatio", function(spec, ppmrange, dppm=0.1)
{
	is <- getseq(spec,ppmrange)
	S <- spec$int
	a <- (S[is[1]] - S[is[length(is)]])/(ppmrange[1]-ppmrange[2])
	b <- S[is[1]] - a*ppmrange[1]
	k <- which(S[is]==min(S[is]))+is[1]
	x0 <- spec$ppm[k]
	Y0 <- 0.9*(a*x0 + b)
	Ymin <- S[k]
	is <- getseq(spec,c(x0-dppm,x0+dppm))
	Ymax <- S[which(S[is]==max(S[is]))+is[1]]
	ifelse( (Ymax/abs(Ymin))>10, 100*(Y0-Ymin)/(Ymax-Y0), 0 )

})

internalClass$set("private", "get_TSP_width", function(spec)
{
	# Local baseline correction 
	BLext <- rep(0,length(spec$int))
	iseq <- getseq(spec,c(-0.5,0.5))
	BLext[iseq] <- Rnmr1D:::.airPLS(spec$int[iseq], lambda=10^9, porder=2)
	spec$int <- spec$int - BLext
	# Find the two points for which the spectrum intensity is equal to A/2 knowing that A is the max amplitude of the TSP peak
	iseq <- getseq(spec,c(-0.05,0.05))
	n <-  iseq[1] + which(spec$int[iseq]==max(spec$int[iseq])) - 1
	Vs <- iseq[1] + which(spec$int[iseq]>=spec$int[n]/2) - 1
	dv <- spec$ppm[ c( Vs[1], Vs[length(Vs)] ) ]
	round(sum(abs(dv))*spec$acq$SFO1,3)
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
	V <- Rnmr1D::Smooth(spec$int[iseq], WS)
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
		BLext[iseq] <- bc # + 2*specSig
   }
   BLext
})

internalClass$set("private", "getBLexternal", function(spec, blset, porder=1, ppmrange=NULL, WS=0)
{
	# external baseline at zero by default
	BLext <- rep(0,length(spec$int))
	if (blset != 0 ) {
		cmax <- switch( min(max(round(porder),1),3), 6, 7, 8 )
		lambda <- ifelse( blset>0, cmax-blset, abs(blset) )
		# Estimation of a baseline correction
		if (is.null(ppmrange)) {
			Y <- spec$int
			if (WS>0) Y <- Rnmr1D:::Smooth(Y,WS)
			BLext <- Rnmr1D:::.airPLS(Y, lambda=10^lambda, porder=porder) - 2*spec$B
		} else {
			iseq <- getseq(spec,ppmrange)
			Y <- spec$int[iseq]
			if (WS>0) Y <- Rnmr1D:::Smooth(Y,WS)
			BLext[iseq] <- Rnmr1D:::.airPLS(Y, lambda=10^lambda, porder=porder) + spec$B
		}
	}
	BLext
})

# Apply a filter on the spectrum - filter must be in daub8, symlet8, smooth0, ..., smooth3
internalClass$set("private", "filterSpectrum", function(spec, filter)
{
	specInt <- spec$int
	if (filter %in% c('daub8', 'symlet8')) {
		specInt <- Rnmr1D::filterByWT(spec$int, filter, threshold = 0.5)
	}
	if (filter %in% length(grep('smooth*', filter))) {
		fsavgol <- Rnmr1D::getDeconvParams()$flist[[filter]]
		specInt <- Rnmr1D::filterSavGol(spec$int, fsavgol$m, fsavgol$nl, fsavgol$nr)
	}
	specInt
})
