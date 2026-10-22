#=====================================================================
# User function to generate the sample table
#=====================================================================

# Get the sample table under RAWDIR
internalClass$set("public", "get_samples_table", function(sequence=NULL, infos=FALSE)
{
	seq <- ifelse(is.null(sequence), self$SEQUENCE, sequence)
	tbl <- get_list_spectrum(RAWDIR, get_list_samples(RAWDIR), sequence=seq)
	tbl <- cbind(tbl[,1:3], rep(0,nrow(tbl)), tbl[,4])
	colnames(tbl) <- c('Spectrum', 'Samplecode', 'EXPNO', 'PROCNO', 'PULSE')
	if (infos) {
		M <- NULL
		for (i in 1:nrow(tbl)) {
			ACQDIR <- file.path(RAWDIR,tbl[i,1],tbl[i,3])
			spec <- applyReadSpectrum(ACQDIR, verbose=0)
			M <- rbind(M,c(get_TSP_width(spec),spec$acq$PULSEWIDTH, spec$acq$NUMBEROFSCANS, spec$acq$SW, spec$proc$SI))	
		}
		colnames(M) <- c('TMSPWIDTH', 'PULSEWIDTH', 'NUMBEROFSCANS', 'SW', 'SI')
		tbl <- cbind(tbl,M)
	}
	tbl
})

#=====================================================================
# User functions for calibration using QC-QS
#=====================================================================

internalClass$set("public", "get_response_factors", function(sampletype,  samplename, thresfP=5, deconv=TRUE, qbl=FALSE, verbose=1)
{
	L2 <- get_list_samples(QSDIR)
	if (length(samplename)==1) {
		samplelist <- L2[ grep(pattern=paste0('^',samplename), L2, value=FALSE) ]
	} else {
		samplelist <- samplename
	}
	if (is.null(get_list_spectrum(QSDIR, samplelist)))
		stop_quietly("Error: ",paste(samplename,collapse=','),"is(are) not a valid spectrum(spectra) name(s).\n")

	if (!sampletype %in% c(QCtype, QStype))
		stop_quietly("Error: sampletype must be either",QCtype,"or",QStype,".\n")

	check_calibration(sequence=SEQUENCE, verbose=(verbose>1))

	if (verbose) cat(samplename, "/",SEQUENCE,"...\n");
	stds_profil_sub <- CALIBRATION[ CALIBRATION$Type==sampletype, , drop=F]

	out <- standardQuantification(stds_profil_sub, samplelist, thresfP, deconv, qbl=qbl, verbose=(verbose>1))
	V <- apply(out$fP,1,mean)
	fP_CV <- sd(V)/mean(V)
	fPUL <- list(mean=mean(V), CV=round(100*fP_CV,2))
	obj <- list(sampletype=sampletype, fPUL=fPUL, fP=out$fP, fR=out$fR, MC=out$MC, INTG=out$INTG, fK=out$fK)
	class(obj) <- 'QC-QS'
	obj
})

internalClass$set("public", "get_factor_table", function(QS)
{
	if (! 'QC-QS' %in% class(QS))
		stop_quietly("Error: You must provide a QC-QS type object\n")

	check_calibration(sequence=SEQUENCE, verbose=FALSE)
	if (!QS$sampletype %in% c(QCtype, QStype))
		stop_quietly("Error: sampletype must be either",QCtype,"or",QStype,".\n")

	stds_loc <- CALIBRATION[! CALIBRATION$Compound %in% calib_cmpd_names, , drop=F]
    stds_loc <- stds_loc[ stds_loc$Type==QS$sampletype, , drop=F]
    compounds <- stds_loc[, 2]
	fPl <- QS$fP
	if (is.null(fPl))
		stop_quietly(paste("Error: all QS spectra appear to have failed.\n"))
    M <- cbind(fPl, apply(fPl,1,mean), 100*apply(fPl,1,sd)/apply(fPl,1,mean))
    colnames(M) <- c(compounds, 'Mean', 'CV%')
    rownames(M) <- rownames(fPl)
    M
})

internalClass$set("public", "get_QC_estimation", function(QC, QS, merge=TRUE, verbose=1)
{
	if (! 'QC-QS' %in% class(QC))
		stop_quietly("Error: You must provide a QC-QS type object\n")
	if (! 'QC-QS' %in% class(QS))
		stop_quietly("Error: You must provide a QC-QS type object\n")

	check_calibration(sequence=SEQUENCE, verbose=FALSE)
	if (verbose) cat("Quality Control for",SEQUENCE,"sequence\n")

	Cest <- apply(QC$fR,2,mean)/QS$fPUL$mean
	Creal <- QC$MC
	Diff <- round(100*(Cest/Creal-1),2)
	stds_profil_QC <- CALIBRATION[ CALIBRATION$Type==QCtype, , drop=F]
	stds_profil_QC <- stds_profil_QC[! stds_profil_QC$Compound %in% calib_cmpd_names, , drop=F]
	M <- cbind(Creal, Cest, Diff)
	rownames(M) <- stds_profil_QC$Compound
	colnames(M) <- c('Real','Estimated', '%Diff')
	if (merge) {
		cmpds <- unique(simplify2array(lapply(rownames(M), function(s) gsub("[0-9]$","",s))))
		M2 <- t(simplify2array(lapply(1:length(cmpds),
				function(k) apply(M[grepl(cmpds[k], rownames(M)),, drop=F],2,mean))))
		rownames(M2) <- cmpds
		M <- M2
	}
	cat('R-squared =',round(stats::cor(Cest,Creal),3),"\n")
	M
})

internalClass$set("public", "plot_QC_estimation", function(QCest)
{
	if (! 'matrix' %in% class(QCest) || ncol(QCest)!=3)
		stop_quietly("Error: You must provide a suitable matrix\n")

	df1 <- as.data.frame(QCest)
	Yest = lm(data=df1, Estimated~Real)
	df1$model <- coef(Yest)[2]*df1$Real + coef(Yest)[1]
	df1$model2 <- coef(Yest)[2]*df1$Real
	cat ('Rate =',round(coef(Yest)[2],4), ', Intercept =',round(coef(Yest)[1],4),"\n")
	fig1 <- plotly::plot_ly() |>
		plotly::add_trace(data = df1, x = ~Real, y = ~Estimated, name = "Estimated", type = 'scatter', mode='markers') |>
		plotly::add_trace(data = df1, x = ~Real, y = ~model, name = "Model", type = 'scatter', mode='lines') |>
		plotly::add_trace(data = df1, x = ~Real, y = ~model2, name = "Model wo Intercept", type = 'scatter', mode='lines') |>
		plotly::layout(
			colorway = c('blue','red','green'),
			title = paste0("QC estimation - ",SEQUENCE)
		)
	fig1
})

#=====================================================================
# User functions using integration and SNR tables
#=====================================================================

# Get the integration matrix
internalClass$set("public", "get_Matrix_Integrals", function()
{
	if (res$proctype != 'integration')
		stop_quietly(paste0("ERROR : Integrals must be computed with the proc_Integrals() method before !\n"))

	if (is.null(res$allquantifs))
		stop_quietly(paste0("ERROR : No integration / quantification !\n"))

	S <-unique(res$allquantifs[,1])
	cmpds <- unique(res$allquantifs[,3])

	# Merge all integration into a matrix
	M <- NULL
	for(k in 1:length(S)) {
		long <- res$allquantifs[res$allquantifs[,1]==S[k],c(2,3,6)]
		colnames(long)[1] <- "name"
		M <- rbind(M, reshape(long, timevar="Compound", direction="wide", idvar="name"))
	}
	M2 <- matrix(round(as.numeric(as.matrix(M[,2:ncol(M)])),2), nrow=nrow(M))
	colnames(M2) <- cmpds
	rownames(M2) <- M[,1]

	# Add missing samples
	V <- SAMPLES[which(!SAMPLES[,2] %in% rownames(M2)),2]
	if (length(V)>0) {
		M3 <- matrix(rep(NA,length(V)*ncol(M2)), nrow=length(V), ncol=ncol(M2))
		M4 <- rbind(M2,M3)
		rownames(M4)[(nrow(M2)+1):(nrow(M2)+nrow(M3))] <- V
		# Reorder rownames in the same order than samples
		V <- simplify2array(lapply(SAMPLES[,2], function(x){which(rownames(M4)==x)}))
		MatInt <- M4[V, , drop=F]
	} else {
		MatInt <- M2
	}
	MatInt
})

# Get the SNR matrix
internalClass$set("public", "get_Matrix_SNR", function()
{
	if (res$proctype != 'integration')
		stop_quietly(paste0("ERROR : Integrals must be computed with the proc_Integrals() method before !\n"))

	if (is.null(res$allquantifs))
		stop_quietly(paste0("ERROR : No integration / quantification !\n"))

	S <-unique(res$allquantifs[,1])
	cmpds <- unique(res$allquantifs[,3])
	
	# Merge all SNR into a matrix
	M <- NULL
	for(k in 1:length(S)) {
		long <- res$allquantifs[res$allquantifs[,1]==S[k],c(2,3,7)]
		colnames(long)[1] <- "name"
		M <- rbind(M, reshape(long, timevar="Compound", direction="wide", idvar="name"))
	}
	M2 <- matrix(round(as.numeric(as.matrix(M[,2:ncol(M)])),2), nrow=nrow(M))
	colnames(M2) <- cmpds
	rownames(M2) <- M[,1]

	# Add lacked samples
	V <- SAMPLES[which(!SAMPLES[,2] %in% rownames(M2)),2]
	if (length(V)>0) {
		M3 <- matrix(rep(NA,length(V)*ncol(M2)), nrow=length(V), ncol=ncol(M2))
		M4 <- rbind(M2,M3)
		rownames(M4)[(nrow(M2)+1):(nrow(M2)+nrow(M3))] <- V
		# Reorder rownames in the same order than samples
		V <- simplify2array(lapply(SAMPLES[,2], function(x){which(rownames(M4)==x)}))
		MatSNR <- M4[V, , drop=F]
	} else {
		MatSNR <- M2
	}
	MatSNR
})

# Get the integration matrix
internalClass$set("public", "get_Matrix_Stats", function(MatInt)
{
    M <- NULL
    M <- rbind(M, as.numeric(apply(MatInt,2, function(v){ min(v, na.rm=T)})))
    M <- rbind(M, as.numeric(apply(MatInt,2, function(v){ max(v, na.rm=T)})))
    M <- rbind(M, as.numeric(apply(MatInt,2, function(v){ mean(v, na.rm=T)})))
    M <- rbind(M, as.numeric(100*apply(MatInt,2, function(v){
                    sd(v, na.rm=T)})/apply(MatInt,2, function(v){ mean(v, na.rm=T)})))
    colnames(M) <- colnames(MatInt)
    rownames(M) <- c('Min', 'Max', 'Mean', 'CV%')
    M
})

# Calculate the CV% for each sample
internalClass$set("public", "get_Matrix_CV", function(MatInt=NULL)
{
	if (is.null(MatInt))
		MatInt <- get_Matrix_Integrals()
	SL <- unique(SAMPLES[,1])
	MatCV <- matrix(0, nrow=length(SL), ncol=ncol(MatInt))
	for (k in 1:length(SL)) {
		SC <- SAMPLES[SAMPLES[,1] == SL[k],2]
		M1 <- MatInt[rownames(MatInt) %in% SC, , drop=FALSE]
		if (is.null(M1) || nrow(M1)==0) next
		MatCV[k, 1:ncol(M1)] <- sapply(1:ncol(M1), function(x){
				V <-M1[ !is.na(M1[,x]), x ]
				ifelse( length(V)>0, 100*sd(V)/mean(V), NA )
		})
	}
	colnames(MatCV) <- colnames(MatInt)
	rownames(MatCV) <- SL
	MatCV
})

# Get a summary 
internalClass$set("public", "get_Matrix_Summary", function()
{
	if (is.null(res$allquantifs))
		stop_quietly(paste0("ERROR : No integration / quantification !\n"))

	M1 <- res$infos[,c(1,3,4,5,6,11,12)]
	M2 <- res$allquantifs[,c(2,3,6,7)]
	CZ <- unique(PROFILE$quantif[PROFILE$quantif$compound %in% M2$Compound, c('compound','zone'), drop=F])
	M <- NULL
	for(S in unique(M1$Samplecode)) {
		M1b <- M1[M1[,1] == S,,drop=F ]
		M2b <- M2[M2[,1] == S,,drop=F ]
		for (k in 1:nrow(M1b)) {
			n <- sum(CZ$zone %in% M1b$zone[k])
			M <- rbind(M, cbind( matrix(rep(M1b[k,],n), ncol=ncol(M1b), nrow=n, byrow=T), 
					M2b[ M2b$Compound %in% unique(CZ[CZ$zone %in% M1b$zone[k], 1]), -1, drop=F ] ))
		}
	}
	colnames(M) <- c( colnames(M1), colnames(M2)[-1] )
	M$ppm1 <- unlist(M$ppm1)
	M$ppm2 <- unlist(M$ppm2)
	M$zone <- unlist(M$zone)
	M$nbpeaks <- unlist(M$nbpeaks)
	M$R2 <- unlist(M$R2)
	M['Int. Diff. %'] <- unlist(M['Int. Diff. %'])
	M$Integral <- unlist(M$Integral)
	M$SNR <- unlist(M$SNR)
	rownames(M) <- NULL
	M
})

internalClass$set("public", "save_Matrices", function(file, filelist=NULL)
{
	# Styles
	styBH <- openxlsx::createStyle(fgFill = "#0070C0", halign = "CENTER", textDecoration = "Bold", border = "Bottom", fontColour = "white")
	styBOLD2 <- openxlsx::createStyle(textDecoration = "Bold")
	styWrap <- openxlsx::createStyle(wrapText = TRUE)
	centerStyle <- openxlsx::createStyle(halign = "center")
	styWarn <- openxlsx::createStyle(fontColour = "#FFFFFF", bgFill = "#FF0000")

	if (res$proctype != 'integration')
		stop_quietly(paste0("ERROR : Integrals must be computed with the proc_Integrals() method before !\n"))

	# Get all result tables
	results <- list(Int=get_NumMat(get_Matrix_Integrals()), SNR=get_NumMat(get_Matrix_SNR()), infos=res$infos)

	# Create Workbook
	wb <- openxlsx::createWorkbook()

	# Create tabs
	tabs <- c( "Samples", "Integrals", "SNR", "Profile", "Infos", "About")
	for (i in 1:length(tabs))  openxlsx::addWorksheet(wb = wb, sheetName = tabs[i], gridLines = TRUE)

	# Write Samples
	Tid <- 1
	samples <- SAMPLES[which(SAMPLES[,2] %in% rownames(results$Int)), ]
	infos <- t(sapply(1:length(specList), function(x) {
		spec <- specList[[x]]
		c(spec$TSPwidth, spec$acq$PULSEWIDTH, spec$acq$NUMBEROFSCANS, spec$acq$SW, spec$proc$SI)
	}))
	colnames(infos) <- c('TSPwidth', 'PULSEWIDTH', 'NUMBEROFSCANS', 'SW', 'SI' )
	samples <- cbind(samples, infos)
	openxlsx::writeData(wb, Tid, x = samples, colNames=TRUE, rowNames=FALSE, withFilter = FALSE)
	openxlsx::addStyle(wb, Tid, style = styBH, rows = 1, cols = c(1:ncol(samples)), gridExpand = TRUE)
	openxlsx::conditionalFormatting(wb, Tid, rows = c(2:(nrow(samples)+1)), cols = which(colnames(samples)=='TSPwidth'),
						rule=paste0(">",TSPwidthMax), style = styWarn)
	openxlsx::setColWidths(wb, Tid, cols=1, widths=30,  ignoreMergedCells = FALSE)

	# Write Integration
	Tid <- Tid + 1
	M <- cbind(rownames(results$Int), results$Int)
	colnames(M)[1] <- 'Samplecode'
	openxlsx::writeData(wb, Tid, x = M, colNames=TRUE, rowNames=FALSE, withFilter = FALSE)
	openxlsx::addStyle(wb, Tid, style = styBH, rows = 1, cols = c(1:ncol(M)), gridExpand = TRUE)

	# Write SNR
	Tid <- Tid + 1
	M <- cbind(rownames(results$SNR), results$SNR)
	colnames(M)[1] <- 'Samplecode'
	openxlsx::writeData(wb, Tid, x = M, colNames=TRUE, rowNames=FALSE, withFilter = FALSE)
	openxlsx::addStyle(wb, Tid, style = styBH, rows = 1, cols = c(1:ncol(M)), gridExpand = TRUE)

	# Write Quantification profile
	Tid <- Tid + 1
	# Preprocess
	startRow <- 1
	preprocess <- t(PROFILE$preprocess)
	openxlsx::writeData(wb, Tid, x = preprocess, startRow=startRow, colNames=TRUE, rowNames=FALSE, withFilter = FALSE)
	openxlsx::addStyle(wb, Tid, style = styBH, rows = startRow, cols = c(1:ncol(preprocess)), gridExpand = TRUE)
	# fitting
	startRow <- 5
	fitting <- PROFILE$fitting[PROFILE$fitting$zone %in% res$zones, ]
	fitting$obl <- as.numeric(fitting$obl)
	openxlsx::writeData(wb, Tid, x = fitting, startRow=startRow, colNames=TRUE, rowNames=FALSE, withFilter = FALSE)
	openxlsx::addStyle(wb, Tid, style = styBH, rows = startRow, cols = c(1:ncol(fitting)), gridExpand = TRUE)
	# quantif
	startRow <- startRow + nrow(fitting) + 3
	quantif <- as.data.frame(PROFILE$quantif[PROFILE$quantif$zone %in% res$zones, ])
	openxlsx::writeData(wb, Tid, x = quantif, startRow=startRow, colNames=TRUE, rowNames=FALSE, withFilter = FALSE)
	SR <- startRow; ER <- startRow+nrow(quantif)
	openxlsx::addStyle(wb, Tid, style = centerStyle, rows = c(SR:ER), cols = c(2:10), gridExpand = TRUE)
	openxlsx::addStyle(wb, Tid, style = styBH, rows = startRow, cols = c(1:ncol(quantif)), gridExpand = TRUE)
	openxlsx::setColWidths(wb, Tid, cols=1, widths=25,  ignoreMergedCells = FALSE)
	openxlsx::setColWidths(wb, Tid, cols=5:6, widths=17, ignoreMergedCells = FALSE)

	# Write Infos
	Tid <- Tid + 1
	openxlsx::writeData(wb, Tid, x = results$infos, colNames=TRUE, rowNames=FALSE, withFilter = FALSE)
	openxlsx::addStyle(wb, Tid, style = styBH, rows = 1, cols = c(1:ncol(results$infos)), gridExpand = TRUE)

	# About
	infos <- rbind(
		c("Raw Data", RAWDIR),
		c("SampleFile", ifelse(is.list(filelist) && !is.null(filelist$SAMPLEFILE), filelist$SAMPLEFILE,'-')),
		c("Profile", ifelse(is.list(filelist) && !is.null(filelist$PROFILE), filelist$PROFILE,'-')),
		c("Instrument Field", FIELD),
		c("Sample Type", TYPE),
		c("Pulse Sequence", SEQUENCE),
		c("",""),
		c("",""),
		c("Processing",""),
		c("Running date", date()),
		c("Nb samples", nrow(SAMPLES)),
		c("Zone numbers", paste(res$zones,collapse=',')),
		c("CPU numbers",res$ncpu),
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
	openxlsx::addStyle(wb, Tid, style = styBOLD2, rows = c(10,NL+2), cols = 1, gridExpand = TRUE)
	openxlsx::addStyle(wb, Tid, style = styWrap,  rows = c(NL+9,NL+10,NL+11), cols = 2, gridExpand = TRUE)
	openxlsx::setColWidths(wb, Tid, cols=1, widths=30,  ignoreMergedCells = FALSE)
	openxlsx::setColWidths(wb, Tid, cols=2, widths=125, ignoreMergedCells = FALSE)

	# Save Workbook
	openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
})
