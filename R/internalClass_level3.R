#=====================================================================
# User function to generate the sample table
#=====================================================================

# Get the sample table under RAWDIR
internalClass$set("public", "get_samples_table", function(sequence=NULL, infos=FALSE)
{
	tbl <- get_list_spectrum(RAWDIR, get_list_samples(RAWDIR))
	seq <- ifelse(is.null(sequence), self$SEQUENCE, sequence)
	tbl <- tbl[ tbl[,3] %in% seq, ]
	tbl <- cbind(tbl[,1], sapply(1:nrow(tbl), function(x) { paste(tbl[x,1:2], collapse='-') }), tbl[,2:3])
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

internalClass$set("public", "get_response_factors", function(sampletype,  samplename, thresfP=5, deconv=TRUE, verbose=1)
{
	if (is.null(get_list_spectrum(QSDIR, samplename)))
		stop_quietly("Error: ",samplename,"is not a valid spectrum name.\n")

	if (!sampletype %in% c(QCtype, QStype))
		stop_quietly("Error: sampletype must be either",QCtype,"or",QStype,".\n")

	if (sampletype == QCtype) check_calibration(QC=samplename, sequence=SEQUENCE, verbose=(verbose>1))
	if (sampletype == QStype) check_calibration(QS=samplename, sequence=SEQUENCE, verbose=(verbose>1))
	
	if (verbose) cat(samplename, "/",SEQUENCE,"...\n");
	stds_profil_sub <- CALIBRATION[ CALIBRATION$Type==sampletype, , drop=F]

	out <- standardQuantification(stds_profil_sub, samplename, thresfP, deconv=deconv, verbose=(verbose>1))
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

internalClass$set("public", "save_Matrices", function(file, filelist=NULL)
{
	# Styles
	styBH <- openxlsx::createStyle(fgFill = "#0070C0", halign = "CENTER", textDecoration = "Bold", border = "Bottom", fontColour = "white")
	styBOLD2 <- openxlsx::createStyle(textDecoration = "Bold")
	styWrap <- openxlsx::createStyle(wrapText = TRUE)

	if (res$proctype != 'integration')
		stop_quietly(paste0("ERROR : Integrals must be computed with the proc_Integrals() method before !\n"))

	# Get all result tables
	results <- list(Int=get_NumMat(get_Matrix_Integrals()), SNR=get_NumMat(get_Matrix_SNR()), infos=res$infos)

	# Create Workbook
	wb <- openxlsx::createWorkbook()

	# Create tabs
	tabs <- c( "Integrals", "SNR", "Infos", "About")
	for (i in 1:length(tabs))  openxlsx::addWorksheet(wb = wb, sheetName = tabs[i], gridLines = TRUE)

	# Write Integration
	Tid <- 1
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
		c("Wine Type", TYPE),
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

#=====================================================================
# User functions for Spectra Visualisation
#=====================================================================

# View spectra along with models & compounds
internalClass$set("public", "view_spectra", function (id, plotmodel=TRUE, plotTrueSpec=TRUE, plotresidus=FALSE, plotzones=TRUE, tags='none', lw=2, showlegend=TRUE, legendhoriz=FALSE, showgrid=TRUE, title=NULL, colspecs=NULL, colcpmds=NULL, verbose=FALSE)
{
	if (! res$proctype %in% c('integration', 'quantification') || length(specList)==0)
		stop_quietly(paste0("ERROR : Integrals or quantification must be computed before !\n"))

	S <- ifelse(is.numeric(id), SAMPLES[id,2], id)
	idx <- which(sapply(1:nrow(SAMPLES), function(k){ specList[[k]]$samplecode==S }))
	spec <- specList[[idx]]
	peaklist <- res$peaklist

	# Compound colors - see https://derekogle.com/NCGraphing/resources/colors
	if (is.null(colcpmds))
		colcpmds <- c('mediumorchid1','palegreen4','deepskyblue3','lightsalmon2','steelblue4',
					'lightpink3','purple','blue','magenta','green','chocolate','chartreuse',
					'bisque3','firebrick3','slateblue2')
	colcpmds <- c( colcpmds, colcpmds, colcpmds )

	ppmview <- ppm_range
	ppm <- spec$ppm

	if (plotTrueSpec) {
		ycurves <- cbind(spec$int, spec$fit$Ymodel)
	} else {
		ycurves <- cbind(spec$fit$Y, spec$fit$Ymodel)
	}
	if (is.null(colspecs) || length(colspecs)<2)
		colspecs <- c('grey60',ifelse(spec$TSPwidth>TSPwidthMax,'violetred','lightslateblue'));
	ynames <- c( S, 'model' )

	if (plotresidus) {
		ycurves <- cbind(ycurves, spec$fit$Y - spec$fit$Ymodel)
		ynames <- c(ynames, 'residus')
		if (length(colspecs)<3) colspecs <- c(colspecs, 'pink')
		colspecs <- colspecs[1:3]
	} else {
		colspecs <- colspecs[1:2]
	}

	names(colspecs) <- ynames
	if (is.null(title)) title <- S
	p <- Rnmr1D::plotSpec(ppmview, ppm, ycurves,  ynames, ycolors=colspecs, lw=lw, title=title)
	arrColors <- colspecs

  # Get the fitting zones
	fit <- PROFILE$fitting
	fit <- fit[ fit$ppm1>=ppmview[1], ]
	fit <- fit[ fit$ppm2<=ppmview[2], ]

  # Get the peaklist corresponding to a quantif zone
	sok <- TRUE
	PL <- cmpdlist <- NULL
	if (sok && (plotmodel || tags != 'none')) {
		cmpdlist <- peaklist[peaklist[,2]==S,3:4, drop=F]
		cmpdlist <- cmpdlist[!is.na(cmpdlist[,2]), , drop=F]
		cmpd_profil <- PROFILE$quantif[PROFILE$quantif$zone %in% fit$zone, ,drop=F]$compound
		cmpdlist <- cmpdlist[cmpdlist[,1] %in% cmpd_profil, , drop=F]
		sok <- FALSE
		if (nrow(cmpdlist)>0) {
			strlist <- paste0(cmpdlist[,2], collapse=",")
			idpeaks <- sort(as.numeric(unlist(strsplit(strlist,","))))
			PL <- spec$fit$peaks[idpeaks, ]
			if (!is.null(PL)) sok <- TRUE
		}
	}
	if (verbose && ! is.null(cmpdlist) && nrow(cmpdlist)>0) { print(cmpdlist); cat("\n") }
	if (verbose && ! is.null(PL) && nrow(PL)>0) { print(PL); cat("\n") }

	if (sok && nrow(fit)>0 && (plotmodel||plotzones)) {

		if (plotmodel) {
			for (k in 1:nrow(cmpdlist)) {
				strlist <- paste0(cmpdlist[k,2], collapse=",")
				idpeaks <- sort(as.numeric(unlist(strsplit(strlist,","))))
				PLk <- spec$fit$peaks[idpeaks, , drop=F]
				if (nrow(PLk)==0) next
				M <- PROFILE$fitting[PROFILE$fitting$zone %in% fit$zone, 1:2]
				ppm_zone <- c(min(M[,1]), max(M[,2]))
				iseq <- getseq(spec, ppm_zone)
				V <- simplify2array(lapply(1:nrow(PLk), function(i) {
						Rnmr1D::PVoigt(ppm[iseq], PLk$amp[i], PLk$ppm[i], PLk$sigma[i], PLk$asym[i], PLk$eta[i])}))
				fmodel <- rep(0, length(ppm))
				fmodel[iseq] <- apply(V,1,sum)
				if (sum(fmodel)>0 && min(ppm[iseq])>=ppmview[1] && max(ppm[iseq])<=ppmview[2]) {
					df <- data.frame(x=ppm[iseq], y=fmodel[iseq])
					p <- plotly::add_trace(p, data=df, x = ~x, y = ~y, name=cmpdlist[k,1], mode = 'lines', fill = 'tozeroy', fillcolor=colcpmds[k])
					arrColors <- c(arrColors, colcpmds[k])
					names(arrColors)[length(arrColors)] <- cmpdlist[k,1]
				}
			}
		}

		if (plotzones) {
			iseq <- getseq(spec,ppmview)
			Ymax <- max(spec$int[iseq])
			lshapes <- list()
			for (k in 1:nrow(fit)) {
					lshapes[[k]] <- list(type="rect", fillcolor="blue", line=list(color="blue"), opacity=0.2,
									x0 = fit$ppm1[k], x1 = fit$ppm2[k], y0 = 0, y1 = Ymax)
			}
			p <- plotly::layout(p, shapes=lshapes)
		}

		if (tags %in% c('id','name','auto') && nrow(cmpdlist)>0 && nrow(PL)>0) {
			if (tags=='auto')
				tags <- ifelse( (ppmview[2]-ppmview[1])>1.5, 'id', 'name' )
			M <- NULL
			for (k in 1:nrow(PL)) {
				pid <- rownames(PL)[k]
				V <- simplify2array(lapply(1:nrow(cmpdlist),
						function(i) { which(pid %in% as.numeric(unlist(strsplit(cmpdlist[i,2],",")))) }))
				M <- rbind(M, c(PL[k,2],1.05*PL[k,3], cmpdlist[which(as.numeric(V)==1),1]))
			}
			data <- data.frame(lab=M[,3], x=M[,1], y=M[,2], tags=sapply(1:nrow(M), function(k) {which(unique(M[,3]) == M[k,3])}))
			if (tags=='name')
				p <- p |> plotly::add_annotations(x = data$x, y = data$y, text = as.character(data$lab),
					showarrow = TRUE, arrowcolor='red', textangle=-30,
					font = list(color = 'black', family = 'sans serif', size = 16))
			if (tags=='id')
				p <- p |> plotly::add_annotations(x = data$x, y = data$y, text = as.character(data$tags),
					showarrow = TRUE, arrowcolor='red', hovertext=data$lab,
					hoverlabel=list(font = list(color = 'blue', family = 'sans serif', size = 18)),
					font = list(color = 'black', family = 'sans serif', size = 12))
		}
	}

	if (legendhoriz)
		p <- p |> plotly::layout(legend = list(orientation = 'h', xanchor = "center",  x = 0.5))
	if (!showgrid)
		p <- p |> plotly::layout(xaxis = list(showgrid = F), yaxis = list(showgrid = F))

	p <- p |> plotly::layout(colorway = arrColors, showlegend=showlegend)
	p
})

