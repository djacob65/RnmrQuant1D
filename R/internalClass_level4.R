#=====================================================================
#  User function  for Integration computing
#=====================================================================

internalClass$set("public", "proc_Integrals", function(zones, ncpu=2, verbose=1)
{
	# Check the quantification profile
	check_profile(zones)

	pkfit <- PROFILE$fitting
	if (!is.null(zones)) pkfit <- pkfit[ pkfit$zone %in% zones, , drop=F]
	if (nrow(pkfit)==0)
		stop_quietly(paste0("Error: there no zone(s) '",paste0(collapse=","),"' defined in the fitting section in the quantification profile"))
	ppm_range <<- c(min(pkfit[,1]), max(pkfit[,2]))

	# Function of combining results for parallelization
	combine_list <- function(LL1, LL2) {
		getList <- function(L) { list(id=L$id, spec=L$spec, quantif=L$quantif, peaklist=L$peaklist, infos=L$infos) }
		mylist <- list()
		for ( i in 1:length(LL1) ) mylist[[LL1[[i]]$id]] <- getList(LL1[[i]])
		for ( i in 1:length(LL2) ) mylist[[LL2[[i]]$id]] <- getList(LL2[[i]])
		return(mylist)
	}

	# Start Cluster
	rq1d <<- self
	if (ncpu>1) {
		cltype <- ifelse(.Platform$OS.type == "unix", "FORK", "PSOCK")
		ncpu <- min(nrow(SAMPLES), max(ncpu,2))
		LOGFILE <- file.path(TMPDIR,"proc_Integrals_par.log")
		if (file.exists(LOGFILE)) unlink(LOGFILE)
		cl <- parallel::makeCluster(ncpu, type=cltype, outfile=LOGFILE)
		doParallel::registerDoParallel(cl)
		parallel::clusterExport(cl=cl, varlist=c("rq1d"), envir=globalenv())
		on.exit(parallel::stopCluster(cl))
	} else {
		foreach::registerDoSEQ()
	}
	Sys.sleep(1)

	Slist <- get_list_spectrum(RAWDIR,get_list_samples(RAWDIR))
	Slist <- Slist[ Slist[,1] %in% SAMPLES[,1] & Slist[,3] == SEQUENCE, 1:2]
	Slist <- Slist[paste0(Slist[,1], Slist[,2],sep="-") %in% paste0(SAMPLES[,1], SAMPLES[,3],sep="-"), ]

# ,.export=c("RnmrQuant1D")
	# Process all samples in parallel
	out <- foreach::foreach(ID=1:nrow(Slist),
		.combine=combine_list,
		.packages=c("Rnmr1D")) %dopar%
	{
		priv <- rq1d$.__enclos_env__$private

		# Directory of the raw spectrum
		samplename <- paste(Slist[ID,],collapse="/")
		ACQDIR <- file.path(rq1d$RAWDIR,samplename)
		print(paste(ID,"/",nrow(Slist),"..."))
		sink(file = paste0(rq1d$TMPDIR,"/log-",Slist[ID,1],"-",ID,".txt"), append = FALSE, type = c("output", "message"), split = FALSE)

		# Read the preprocessed spectrum (1r)
		cat(samplename,': ')
		spec <- priv$applyReadSpectrum(ACQDIR, verbose=verbose)
		if (verbose) cat("Path:", ACQDIR,"\n")
		if (verbose) cat("Sequence:",spec$acq$PULSE,"\n")
		spec$samplename <- Slist[ID,1]
		spec$expno <- Slist[ID,2]
		SID <- paste0(spec$samplename, spec$expno, sep="-")
		spec$samplecode <- rq1d$SAMPLES[which(paste0(rq1d$SAMPLES[,1], rq1d$SAMPLES[,3],sep="-") %in% SID ),2]

		# Baseline correction
		spec <- priv$applyBLcorrection(spec, verbose=verbose)

		# TSP width
		spec$TSPwidth <- priv$get_TSP_width(spec)
		if (verbose) cat("TSP width:", spec$TSPwidth,"Hz\n")
		if (spec$TSPwidth>rq1d$TSPwidthMax && verbose)
			cat("ERROR : TSP width too large\n")

		# Peak fitting
		Opars <- rq1d$opars
		Opars$peaks <- NULL
		spec <- priv$applyPeakFitting(spec, Opars, zones=zones, ncpu=1, verbose=verbose)

		# Quantification
		if (!is.null(spec$fit$peaks)) {
			if (verbose) { print(spec$fit$infos); cat("\n") }
			Q <- priv$applyQuantification(spec, verbose=verbose)
			print(Q$quantification); cat("\n")
			print(Q$peaklist); cat("\n")
		}
		sink()

		# Accumulating results
		mylist <- list()
		if (!is.null(spec$fit$peaks)) {
			Tinfos <- cbind( rep(spec$samplecode, nrow(spec$fit$infos)), spec$fit$infos )
			colnames(Tinfos)[1] <- 'Samplecode'
			mylist[[paste0('S',ID)]] <- list( id=paste0('S',ID), spec=spec, quantif=Q$quantification, peaklist=Q$peaklist, infos=Tinfos)
		} else {
			mylist[[paste0('S',ID)]] <- list( id=paste0('S',ID), spec=spec, quantif=NULL, peaklist=NULL, infos=NULL )
		}
		return(mylist)
	}

	# Clean memory
	invisible(gc())

	# Results proc-process
	res <<- list(allquantifs=NULL, peaklist=NULL, infos=NULL, zones=zones, ncpu=ncpu, proctype='integration')
	specList <<- list()

	# Merging results
	for(k in 1:length(out)) {
		spec <- out[[k]]$spec
		if (!is.null(out[[k]]$quantif)) {
			M <- matrix(rep('.',nrow(out[[k]]$quantif)*2), ncol=2)
			M[,1] <- spec$samplename; M[,2] <- spec$samplecode
			res$allquantifs <<- rbind(res$allquantifs, cbind(M, out[[k]]$quantif))
			res$peaklist <<- rbind(res$peaklist, cbind(M, out[[k]]$peaklist))
			res$infos <<- rbind(res$infos, out[[k]]$infos)
		}
		specList[[k]] <<- spec
	}
})


#=====================================================================
# User functions for quantification computing
#=====================================================================

internalClass$set("public", "proc_fPULCON", function(QSname, thresfP=5, deconv=TRUE, verbose=1)
{
	if (verbose) cat("Compute the PULCON factor ...\n")
	check_calibration(QS=QSname)

	LogFile <- paste0(TMPDIR,'/stds_',QSname,'-',SEQUENCE,'.txt')
	stds <- CALIBRATION
	stds <- stds[ stds$Type==QStype, , drop=F]
	stds <- stds[! stds$Compound %in% calib_cmpd_names, , drop=F]
	unlink(LogFile)
	t <- system.time({
		sink(LogFile)
		calib <- standardQuantification(stds, QSname, thresfP, deconv, verbose)
		if (is.nan(calib$fPUL$mean)) {
			fP <<- list()
			stop_quietly(paste0("Error: the spectrum '",QSname,"' does not appear to contain the correct quantification standards."))
		}
		if (verbose) print(calib)
		Mat <- get_factor_table(calib)
		sink()
	})
	fP <<- list(QSname=QSname, Mat=Mat, values=calib$fP, CV=calib$fPUL$CV, mean=calib$fPUL$mean,
			fK=calib$fK, elapsed=round(as.numeric(t[3]),2))
	if (verbose) cat('fP =', round(fP$mean),",  CV =",fP$CV,"\n")
	if (verbose) cat('elapsed time =', round(t[3],2),"\n\n")
})

internalClass$set("public", "proc_Quantification", function(cmpdlist=NULL, zones=NULL, ncpu=2, reset=FALSE, CR=FALSE, verbose=1)
{
	check_all()
	if (is.null(cmpdlist) && ! is.null(zones)) check_profile(zones)

	if (is.null(fP) || length(fP)==0 || is.nan(fP$mean))
		stop_quietly(paste0("Error: the PULCON factor must be computed before"))

	if (verbose) cat("Do quantification ... \n")
	if (verbose && CR) cat("\n")

	# Depending on input, return the correspond zones
	if (!is.null(cmpdlist)) {
		L <-  cmpdlist %in% PROFILE$quantif$compound
		if (length(L) != sum(L)) {
			bad <- paste(cmpdlist[which(L==FALSE)],sep=',')
			stop_quietly(paste0("Error: Some compounds (",bad,") are not defined in the quantification profile\n"))
		}
		zones <- unique(PROFILE$quantif[ PROFILE$quantif$compound %in% cmpdlist, ]$zone)
	}
	if (is.null(zones)) zones <- unique(PROFILE$quantif$zone)
	if (is.null(cmpdlist))
		cmpdlist <- unique(PROFILE$quantif[PROFILE$quantif$zone %in% zones, ]$compound)

	# delete all Rdata and log files if required
	if (reset) {
		unlink(paste0(RDATADIR,"/*.RData"))
		unlink(paste0(TMPDIR,"/output_*.txt"))
	}

	# For each samples
	samplelist <- SAMPLES[,1]
	expnolist <<- character(0)
	cnt <- 0
	tottime <- 0
	totcnt <- length(samplelist)

	for (samplename in unique(samplelist))
	{
		slist <- get_list_spectrum(RAWDIR, samplename)
		if (!SEQUENCE %in% unique(slist[,3])) next
		slist <- slist[ slist[,3] == SEQUENCE, , drop=F]
		slist <- slist[paste0(slist[,1],slist[,2],sep="-") %in% paste0(SAMPLES[,1],SAMPLES[,3],sep="-"), ]
		sampleid <- which(samplename==samplelist)
		fdil <- SAMPLES[FDILfield][sampleid,1]
		# For each expno
		for (k in 1:nrow(slist))
		{
			cnt <- cnt + 1
			expno <- slist[k,2]
			expnolist <<- c(expnolist, expno)
			RDataFile <- paste0(RDATADIR, '/',samplename,'-',expno,'.RData')
			if (file.exists(RDataFile)) { totcnt <- totcnt - 1; next }
			if (verbose) cat(sprintf("%3d out of %3d - ", cnt, length(samplelist)))
			if (verbose) cat(sprintf("%20s expno = %2s ...",samplename, expno));
			LogFile <- paste0(TMPDIR, '/output_',samplename,'_',expno,'.txt')
			t <- system.time({
				# Compute the (absolute) quantification matrix (metabolites X repetitions)
				sink(file = LogFile, append = FALSE, type = c("output", "message"), split = FALSE)
				# Deconvolution / Quantification
				if (verbose)  cat("\n\n================================================\n\n")
				out <- sampleQuantification(samplename, expno, zones, ncpu=ncpu, verbose=1)
				quantMat <- out$quantMat
				SNR <- out$SNR
				spec <- out$spec
				spec$sampleidx <- sampleid[k]
				spec$peaklist <- spec$Q$peaklist
				spec$Q <- NULL
				absQuantMat <- NULL
				if (!is.null(quantMat)) {
					if (!is.null(cmpdlist)) quantMat <- quantMat[rownames(quantMat) %in% cmpdlist,,drop=F]
					if (!is.null(cmpdlist)) SNR <- SNR[rownames(SNR) %in% cmpdlist,,drop=F]
					absQuantMat <- absSampleQuantification(spec, fP, quantMat, fdil=fdil[k], verbose=1)
				}
				rm(out)
				sink()
			})
			tottime <- tottime + t[3]
			remaintime <- (totcnt-cnt)*tottime/cnt
			if (verbose) cat(sprintf("OK - duration = %2.2f, remaining time = %4d sec - Ended at %-30s\n",
									round(t[3],2),round(remaintime), as.character(Sys.time()+remaintime)))
			if (verbose && CR) cat("\n")
			quantParams <- list(samplename=samplename, expno=expno, type=TYPE, field=FIELD, sampleidx=spec$sampleidx,
								cmpdlist=cmpdlist, zones=zones)
			save(quantParams, quantMat, absQuantMat, SNR, spec, file = RDataFile)
		}
	}
	quantpars <<- list(cmpdlist=cmpdlist, zones=zones, tottime=as.numeric(tottime), ncpu=ncpu)
	res$proctype <<- 'quantification'
})

internalClass$set("private", "get_results_for_all_samples", function()
{
	if (! dir.exists(RDATADIR) )
		stop_quietly(paste0("ERROR : the directory '",RDATADIR,"' does not exist !\n"))

	if (res$proctype != 'quantification')
		stop_quietly(paste0("ERROR : Quantifications must be computed with the proc_Quantification() method before !\n"))

	Rdata <- stringr::str_replace_all(list.files(path = RDATADIR, pattern = "*.RData"), setNames('','.RData'))
	Msnr <- Mint <- Mquant <- NULL
	samplenames <- sampleinfos <- sampletypes <- fitinfos <- NULL
	cmpdlist <- NULL
	SN1 <- SN2 <- NULL
	repeat {
		if (length(Rdata)<1) break
		for (sample in Rdata) {
			DIR <- file.path(RAWDIR,gsub('-[0-9]+$','', sample))
			if (sample == 'RQ1D' || !dir.exists(DIR)) next
			load(file = paste0(file.path(RDATADIR,sample),'.RData'))
			samplenames <- c(samplenames, quantParams[['samplename']])
			infos <- c(spec$TSPwidth, spec$acq$PULSEWIDTH, spec$acq$NUMBEROFSCANS, spec$acq$SW, spec$proc$SI)
			sampleinfos <- rbind(sampleinfos, infos)
			sampletypes <- c(sampletypes, quantParams[['type']])
			fitinfos <- rbind(fitinfos, cbind(rep(quantParams[['samplename']], nrow(spec$fit$infos)),spec$fit$infos))
			if (!is.null(quantMat)) {
				Mint <- rbind(Mint, t(quantMat))
				Msnr <- rbind(Msnr, t(SNR))
				Mquant <- rbind(Mquant, t(absQuantMat))
				cmpdlist <- rownames(quantMat)
				SN1 <- c(SN1, sample)
			} else {
				SN2 <- c(SN2, sample)
			}
		}

		# Add missing samples if needed
		if (!is.null(SN2)) {
			M <- matrix(rep(NA,length(SN2)*ncol(Mquant)), nrow=length(SN2), ncol=ncol(Mquant))
			Mquant <- rbind(Mquant, M)
			Mint <- rbind(Mint, M)
			Msnr <- rbind(Msnr, M)
		}
		Mquant <- matrix(as.numeric(Mquant), nrow=nrow(Mquant))

		# Add rownames and colnames
		SN <- c(SN1, SN2)
		rownames(Mint) <- rownames(Msnr) <- rownames(Mquant) <- SN
		colnames(Mint) <- colnames(Msnr) <- colnames(Mquant) <- cmpdlist

		# Reorder samples if needed
		if (!is.null(SN2)) {
			V <- simplify2array(lapply(Rdata, function(x){which(SN==x)}))
			Mquant <- Mquant[V, , drop=F]
			Mint <- Mint[V, , drop=F]
			Msnr <- Msnr[V, , drop=F]
		}

		colnames(sampleinfos) <- c('TSPwidth', 'PULSEWIDTH', 'NUMBEROFSCANS', 'SW', 'SI')
		rownames(sampleinfos) <- samplenames
		sampletypes <- as.matrix(sampletypes, ncol=1)
		colnames(sampletypes) <- 'Sample Type'
		rownames(sampletypes) <- samplenames
		colnames(fitinfos)[1] <- 'Spectrum'
		break
	}
	list(sampleinfos=sampleinfos, sampletypes=sampletypes, fitinfos=fitinfos, quantif=Mquant, SNR=Msnr, Int=Mint)
})

internalClass$set("public", "get_output_results", function()
{
	out <- get_results_for_all_samples()
	colnames(out$Int) <- colnames(out$SNR) <- colnames(out$quantif) <-
			sapply(colnames(out$quantif), function(x) gsub('-','_',gsub(' ','_',x)))
	out$Int <- get_NumMat(out$Int, rownames=FALSE)
	out$SNR <- get_NumMat(out$SNR, rownames=FALSE)
	out$quantif <- get_NumMat(out$quantif, rownames=FALSE)

	# Samples
	M <- cbind(SAMPLES[,c(1:3)], SAMPLES[, ncol(SAMPLES)], out$sampletypes, out$sampleinfos)
	colnames(M)[4] <- FDILfield
	out$samples <- M
	out$sampleinfos <- out$sampletypes <- NULL

	# Calibration
	out$fP_Mat <- fP$Mat

	# Standard Profile
	out$stds <- CALIBRATION

	# Quantification Profile
	cmpdlist <- quantpars$cmpdlist
	zones <- quantpars$zones
	if (!is.null(cmpdlist)) {
		profil_quantif <- PROFILE$quantif[ PROFILE$quantif$compound %in% cmpdlist, , drop=F ]
	} else {
		profil_quantif <- PROFILE$quantif[ PROFILE$quantif$zone %in% zones, , drop=F ]
	}
	out$profil_preprocess <- t(PROFILE$preprocess)
	out$profil_quantif <- profil_quantif
	out$profil_fitting <- PROFILE$fitting[PROFILE$fitting$zone %in% unique(profil_quantif$zone), , drop=F ]
	out$profil_compound <- PROFILE$compound[PROFILE$compound$name %in% unique(profil_quantif$compound), ]
	out
})

internalClass$set("public", "get_spectra_data", function()
{
	if (! dir.exists(RDATADIR) )
		stop_quietly(paste0("ERROR : the directory '",RDATADIR,"' does not exist !\n"))

	if (res$proctype != 'quantification')
		stop_quietly(paste0("ERROR : Quantifications must be computed with the proc_Quantification() method before !\n"))

	# Results proc-process
	res <<- list(allquantifs=NULL, peaklist=NULL, infos=NULL, proctype='quantification')
	specList <<- list()

	Rdata <- stringr::str_replace_all(list.files(path = RDATADIR, pattern = "*.RData"), setNames('','.RData'))
	if (length(Rdata)>1) for (sample in Rdata) {
		DIR <- file.path(RAWDIR,gsub('-[0-9]+$','', sample))
		if (sample == 'RQ1D' || !dir.exists(DIR)) next
		load(file = paste0(file.path(RDATADIR,sample),'.RData'))
		PL <- spec$peaklist
		if (!is.null(quantParams$cmpdlist)) PL <- PL[PL[,1] %in% quantParams$cmpdlist, , drop=F]
		n <- nrow(PL)
		M <- matrix(rep(0,n*4), nrow=n, ncol=4)
		M[,1:2] <- cbind(rep(SAMPLES[spec$sampleidx,1],n), rep(SAMPLES[spec$sampleidx,2],n))
		M[,3:4] <- PL[,1:2]
		res$peaklist <<- rbind(res$peaklist, M)
		res$infos <<- rbind(res$infos, spec$fit$infos)
		specList[[spec$sampleidx]] <<- spec
		pkfit <- PROFILE$fitting[PROFILE$fitting$zone %in% quantParams$zones, , drop=F]
        ppm_range <<- c( min(pkfit$ppm1), max(pkfit$ppm2))
	}
})

internalClass$set("public", "save_Results", function(file, filelist=NULL)
{
	# Styles
	styBH <- openxlsx::createStyle(fgFill = "#0070C0", halign = "CENTER", textDecoration = "Bold", border = "Bottom", fontColour = "white")
	styBOLD2 <- openxlsx::createStyle(textDecoration = "Bold")
	styWrap <- openxlsx::createStyle(wrapText = TRUE)
	centerStyle <- openxlsx::createStyle(halign = "center")
	styWarn <- openxlsx::createStyle(fontColour = "#FFFFFF", bgFill = "#FF0000")
	styFmt0 <- openxlsx::createStyle(numFmt = "0")
	styFmt2 <- openxlsx::createStyle(numFmt = "0.00")
	styFmt3 <- openxlsx::createStyle(numFmt = "0.000")
	styFmt4 <- openxlsx::createStyle(numFmt = "0.0000")

	if (res$proctype != 'quantification')
		stop_quietly(paste0("ERROR : Quantifications must be computed with the proc_Quantification() method before !\n"))

	# Create Workbook
	wb <- openxlsx::createWorkbook()

	# Get all result tables
	results <- get_output_results()

	# Create tabs
	tabs <- c( "Samples", "Integrals", "SNR", "Quantifications", "Calibration", "Profile", "About")
	for (i in 1:length(tabs))  openxlsx::addWorksheet(wb = wb, sheetName = tabs[i], gridLines = TRUE)

	# Write tabs
	Tid <- 1 # Samples
	Spectrum <- results$samples[,2]
	openxlsx::writeData(wb, Tid, x = results$samples, colNames=TRUE, rowNames=FALSE, withFilter = FALSE)
	openxlsx::addStyle(wb, Tid, style = styBH, rows = 1, cols = c(1:ncol(results$samples)), gridExpand = TRUE)
	openxlsx::conditionalFormatting(wb, Tid, rows = c(2:nrow(results$samples)), cols = ncol(results$samples)-4,
						rule=">1", style = styWarn)
	openxlsx::setColWidths(wb, Tid, cols=1, widths=30,  ignoreMergedCells = FALSE)

	Tid <- Tid + 1 # Integrals
	M <- cbind(Spectrum, results$Int)
	openxlsx::writeData(wb, Tid, x = M, colNames=TRUE, rowNames=FALSE, withFilter = FALSE)
	openxlsx::addStyle(wb, Tid, style = styBH, rows = 1, cols = c(1:ncol(M)), gridExpand = TRUE)
	openxlsx::setColWidths(wb, Tid, cols=1, widths=30,  ignoreMergedCells = FALSE)
	SR <- 2; ER <- nrow(results$Int)+1; SC <- 2; EC <- ncol(results$Int)+1
	openxlsx::addStyle(wb, Tid, style = styFmt0, rows = c(SR:ER), cols = c(SC:EC), gridExpand = TRUE)

	Tid <- Tid + 1 # SNR
	M <- cbind(Spectrum, results$SNR)
	openxlsx::writeData(wb, Tid, x = M, colNames=TRUE, rowNames=FALSE, withFilter = FALSE)
	openxlsx::addStyle(wb, Tid, style = styBH, rows = 1, cols = c(1:ncol(M)), gridExpand = TRUE)
	openxlsx::setColWidths(wb, Tid, cols=1, widths=30,  ignoreMergedCells = FALSE)

	Tid <- Tid + 1 # Quantification
	M <- cbind(Spectrum, results$quantif)
	openxlsx::writeData(wb, Tid, x = M, colNames=TRUE, rowNames=FALSE, withFilter = FALSE)
	openxlsx::addStyle(wb, Tid, style = styBH, rows = 1, cols = c(1:ncol(M)), gridExpand = TRUE)
	openxlsx::setColWidths(wb, Tid, cols=1, widths=30,  ignoreMergedCells = FALSE)
	SR <- 2; ER <- nrow(results$quantif)+1; SC <- 2; EC <- ncol(results$quantif)+1
	openxlsx::addStyle(wb, Tid, style = styFmt4, rows = c(SR:ER), cols = c(SC:EC), gridExpand = TRUE)

	Tid <- Tid + 1 # Calibration profile
	openxlsx::writeData(wb, Tid, x = results$stds, colNames=TRUE, rowNames=FALSE, withFilter = FALSE)
	openxlsx::addStyle(wb, Tid, style = styBH, rows = 1, cols = c(1:ncol(results$stds)), gridExpand = TRUE)
	startRow <- nrow(results$stds) + 3
	openxlsx::writeData(wb, Tid, x = results$fP_Mat, startRow=startRow, colNames=TRUE, rowNames=TRUE, withFilter = FALSE)
	openxlsx::addStyle(wb, Tid, style = styBH, rows = startRow, cols = c(1:(ncol(results$fP_Mat)+1)), gridExpand = TRUE)
	SR <- startRow+1; ER <- startRow+nrow(results$fP_Mat); SC <- 2; EC <- ncol(results$fP_Mat)
	openxlsx::addStyle(wb, Tid, style = styFmt0, rows = c(SR:ER), cols = c(SC:EC), gridExpand = TRUE)
	openxlsx::addStyle(wb, Tid, style = styFmt3, rows = c(SR:ER), cols = ncol(results$fP_Mat)+1, gridExpand = TRUE)


	Tid <- Tid + 1 # Quantification profile
	# Preprocess
	startRow <- 1
	preprocess <- results$profil_preprocess
	openxlsx::writeData(wb, Tid, x = preprocess, startRow=startRow, colNames=TRUE, rowNames=FALSE, withFilter = FALSE)
	openxlsx::addStyle(wb, Tid, style = styBH, rows = startRow, cols = c(1:ncol(preprocess)), gridExpand = TRUE)
	# fitting
	startRow <- 4
	fitting <- results$profil_fitting
	fitting$obl <- as.numeric(fitting$obl)
	openxlsx::writeData(wb, Tid, x = fitting, startRow=startRow, colNames=TRUE, rowNames=FALSE, withFilter = FALSE)
	openxlsx::addStyle(wb, Tid, style = styBH, rows = startRow, cols = c(1:ncol(fitting)), gridExpand = TRUE)
	# quantif
	startRow <- nrow(fitting) + 3
	quantif <- results$profil_quantif
	openxlsx::writeData(wb, Tid, x = quantif, startRow=startRow, colNames=TRUE, rowNames=FALSE, withFilter = FALSE)
	SR <- startRow; ER <- startRow+nrow(quantif)
	openxlsx::addStyle(wb, Tid, style = centerStyle, rows = c(SR:ER), cols = c(2:7), gridExpand = TRUE)
	openxlsx::addStyle(wb, Tid, style = styBH, rows = startRow, cols = c(1:ncol(quantif)), gridExpand = TRUE)
	openxlsx::setColWidths(wb, Tid, cols=1, widths=25,  ignoreMergedCells = FALSE)
	openxlsx::setColWidths(wb, Tid, cols=5, widths=17, ignoreMergedCells = FALSE)
	# compound
	startRow <- startRow + nrow(quantif) + 3
	compounds <- results$profil_compound
	openxlsx::writeData(wb, Tid, x = compounds, startRow=startRow, colNames=TRUE, rowNames=FALSE, withFilter = FALSE)
	openxlsx::addStyle(wb, Tid, style = styBH, rows = startRow, cols = c(1:ncol(compounds)), gridExpand = TRUE)

	# About
	infos <- rbind(
		c("Raw Data",RAWDIR),
		c("SampleFile", ifelse(is.list(filelist) && !is.null(filelist$SAMPLEFILE), filelist$SAMPLEFILE,'-')),
		c("Processing Profile",ifelse(is.list(filelist) && !is.null(filelist$PROFILE), filelist$PROFILE,'-')),
		c("Calibration Profile",ifelse(is.list(filelist) && !is.null(filelist$CALIBRATION), filelist$CALIBRATION,'-')),
		c("Instrument Field",FIELD),
		c("Wine Type",TYPE),
		c("Pulse Sequence",SEQUENCE),
		c("",""),
		c("",""),
		c("Calibration",""),
		c("PULCON Factor", fP$mean),
		c("CV (%)", fP$CV),
		c("Elapsed time (s)",round(fP$elapsed,2)),
		c("",""),
		c("",""),
		c("Processing",""),
		c("Running date", date()),
		c("Nb samples", nrow(SAMPLES)),
		c("PPM range for Noise estimation ", sprintf("%2.1f - %2.1f", PPM_NOISE[1], PPM_NOISE[2])),
		c("CPU numbers",quantpars$ncpu),
		c("Elapsed time (s)",round(quantpars$tottime,2)),
		c("",""),
		c("","")
	)
	NL <- nrow(infos)

	# Environment
	V <- sessionInfo()
	p <- ls(V$loadedOnly)
	packages <- NULL
	for (i in 1:length(p))
		packages <- c( packages, paste0(V$loadedOnly[[p[i]]]$Package,'_',V$loadedOnly[[p[i]]]$Version) )
	p <- ls(V$otherPkgs)
	others <- NULL
	for (i in 1:length(p))
		others <- c( others, paste0(V$otherPkgs[[p[i]]]$Package,'_',V$otherPkgs[[p[i]]]$Version) )
	infos <- rbind(infos,
		c("Environment",""),
		c("R version", gsub("R version ","", V$R.version$version.string)),
		c("Running under", V$running),
		c("Platform", V$platform),
		c("Blas", V$BLAS),
		c("Lapack", V$LAPACK),
		c("Locale", gsub(';',', ', V$locale)),
		c("Base Packages", paste(V$basePkgs, collapse=', ')),
		c("Loarded Packages", paste(packages, collapse=', ')),
		c("Others Packages", paste(others, collapse=', ')) )
	colnames(infos) <- c("Label","Value")
	Tid <- Tid + 1
	openxlsx::writeData(wb, Tid, x = infos,  colNames=TRUE, rowNames=FALSE, withFilter = FALSE)
	openxlsx::addStyle(wb, Tid, style = styBH,    rows = 1, cols = c(1:ncol(infos)), gridExpand = TRUE)
	openxlsx::addStyle(wb, Tid, style = styBOLD2, rows = c(11,17, NL+2), cols = 1, gridExpand = TRUE)
	openxlsx::addStyle(wb, Tid, style = styWrap,  rows = c(NL+9,NL+10,NL+11), cols = 2, gridExpand = TRUE)
	openxlsx::setColWidths(wb, Tid, cols=1, widths=30,  ignoreMergedCells = FALSE)
	openxlsx::setColWidths(wb, Tid, cols=2, widths=125, ignoreMergedCells = FALSE)

	# Save Workbook
	openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
})

#=====================================================================
# User functions for Spectra visualisation
#=====================================================================

internalClass$set("public", "plot_spectra", function(id, compound, ...)
{
	if (sum(compound %in% PROFILE$compound$name)==0)
		stop_quietly(paste("Error:",compound,"not known as a compound"))

	if (res$proctype != 'quantification')
		stop_quietly(paste0("ERROR : Quantifications must be computed with the proc_Quantification() method before !\n"))

	if (is.null(res$peaklist))
		stop_quietly(paste0("ERROR : Spectra must be retrieved in the RnmrQuant1D instance with the get_spectra_data() method before !\n"))

	zones <- unique(PROFILE$quantif[ PROFILE$quantif$compound %in% compound, ]$zone)
	pkfit <- PROFILE$fitting[PROFILE$fitting$zone %in% zones, ]
	ppm_range <<- c( min(pkfit$ppm1), max(pkfit$ppm2))
	view_spectra(id, ...)
})

