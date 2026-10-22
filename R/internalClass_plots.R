# Defaut Graphic settings
internalClass$set("private", "def_plot_settings", list(
	mode = 'html',
	IMGDIR = 'tmp',
	width = 600,
	height = 600,
	plotTrueSpec = FALSE,
	plotmodel = TRUE,
	plotresidus = FALSE,
	plotzones = FALSE,
	tags = 'none',
	showgrid = TRUE,
	legendhoriz = TRUE,
	legendtop = FALSE,
	showlegend = TRUE,
	yaxis = TRUE,
	ylabel = '', # 'Intensity (a.u)'
	xlabel = '', # 'Shift (ppm)'
	verbose = TRUE,
	lw = 2,
	title = '',
	opacity = 0.7,
	yview = NULL,
	ppmview = NULL,
	colspecs = c('gray70','#86c1db','deeppink4'),
	colcpmds = c('#5ba8c9','dodgerblue1','#5b75c9', 'slateblue2', '#8334b8', 'steelblue1'),
	font = list(family = "Arial", size = 20),
	font2 = list(family = "Arial", size = 18),
	filename = NULL
))

#=====================================================================
# Main function for Spectra Visualisation
#=====================================================================

# View spectra along with models & compounds
internalClass$set("public", "view_spectra", function (id, plotmodel=TRUE, plotTrueSpec=TRUE, plotresidus=FALSE, plotzones=TRUE, tags='none', lw=2, showlegend=TRUE, legendhoriz=FALSE, legendtop=FALSE, showgrid=TRUE, title=NULL, colspecs=NULL, colcpmds=NULL, opacity=0.7, verbose=FALSE)
{
	if (! res$proctype %in% c('integration', 'quantification') || length(specList)==0)
		stop_quietly(paste0("ERROR : Integrals or quantification must be computed before !\n"))

	S <- ifelse(is.numeric(id), SAMPLES[id,2], id)
	idx <- which(sapply(1:nrow(SAMPLES), function(k){ specList[[k]]$samplecode==S }))
	spec <- specList[[idx]]
	peaklist <- res$peaklist

	# default colors - see https://derekogle.com/NCGraphing/resources/colors
	defcols <- c('mediumorchid1','palegreen4','deepskyblue3','lightsalmon2','steelblue4',
					'lightpink3','purple','blue','magenta','green','chocolate','chartreuse',
					'bisque3','firebrick3','slateblue2')

	# Compound colors - see https://derekogle.com/NCGraphing/resources/colors
	if (is.null(colcpmds))
		colcpmds <- defcols
	colcpmds <- c( colcpmds, colcpmds, colcpmds )

	ppmview <- ppm_range
	ppm <- spec$ppm

	# Select the spectrum with or without baseline correction
	if (plotTrueSpec) {
		ycurves <- cbind(spec$int, spec$fit$Ymodel)
	} else {
		ycurves <- cbind(spec$fit$Y, spec$fit$Ymodel)
	}

	# Color with a different color if TSPwith > TSPwidthMax
	if (is.null(colspecs) || length(colspecs)<2)
		colspecs <- c('grey60',ifelse(spec$TSPwidth>TSPwidthMax,'violetred','lightslateblue'));
	ynames <- c( S, 'model' )

	# Add residus if selected
	if (plotresidus) {
		ycurves <- cbind(ycurves, spec$fit$Y - spec$fit$Ymodel)
		ynames <- c(ynames, 'residus')
		if (length(colspecs)<3) colspecs <- c(colspecs, 'pink')
		colspecs <- colspecs[1:3]
	} else {
		colspecs <- colspecs[1:2]
	}

    # Plot spectrum & model
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

	# Plot compounds / zones
	pkcmpds <- NULL
	if (sok && nrow(fit)>0 && (plotmodel||plotzones)) {

		if (plotmodel) {
			for (k in 1:nrow(cmpdlist)) {
				strlist <- paste0(cmpdlist[k,2], collapse=",")
				idpeaks <- sort(as.numeric(unlist(strsplit(strlist,","))))
				PLk <- spec$fit$peaks[idpeaks, , drop=F]
				if (nrow(PLk)==0) next
				pkcmpds <- c(pkcmpds, idpeaks)
				M <- PROFILE$fitting[PROFILE$fitting$zone %in% fit$zone, 1:2]
				ppm_zone <- c(min(M[,1]), max(M[,2]))
				iseq <- getseq(spec, ppm_zone)
				V <- simplify2array(lapply(1:nrow(PLk), function(i) {
						Rnmr1D::PVoigt(ppm[iseq], PLk$amp[i], PLk$ppm[i], PLk$sigma[i], PLk$asym[i], PLk$eta[i])}))
				fmodel <- rep(0, length(ppm))
				fmodel[iseq] <- apply(V,1,sum)
				if (sum(fmodel)>0 && min(ppm[iseq])>=ppmview[1] && max(ppm[iseq])<=ppmview[2]) {
					df <- data.frame(x=ppm[iseq], y=fmodel[iseq])
					fcol <- paste0('rgba(',paste(c(as.vector(col2rgb(colcpmds[k])), opacity), collapse=','),')')
					p <- plotly::add_trace(p, data=df, x = ~x, y = ~y, name=cmpdlist[k,1], mode = 'lines', fill = 'tozeroy', fillcolor = fcol)
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
					lshapes[[k]] <- list(type="rect", fillcolor="gray", line=list(color="gray"), opacity=0.1,
									x0 = fit$ppm1[k], x1 = fit$ppm2[k], y0 = 0, y1 = Ymax)
			}
			p <- plotly::layout(p, shapes=lshapes)
		}
	}

	if (tags == 'allpeaks') {
		M <- PROFILE$fitting[PROFILE$fitting$zone %in% fit$zone, 1:2]
		ppm_zone <- c(min(M[,1]), max(M[,2]))
		iseq <- getseq(spec, ppm_zone)
		PLall <- spec$fit$peaks[ spec$fit$peaks$ppm>=ppmview[1] & spec$fit$peaks$ppm<=ppmview[2], , drop=F]
		for (k in 1:nrow(PLall)) {
			if (!is.null(pkcmpds) && rownames(PLall[k,]) %in% pkcmpds) next
			V <- Rnmr1D::PVoigt(ppm[iseq], PLall$amp[k], PLall$ppm[k], PLall$sigma[k], PLall$asym[k], PLall$eta[k])
			df <- data.frame(x=ppm[iseq], y=V)
			ic <- k %% length(defcols) + 1
			fcol <- paste0('rgba(',paste(c(as.vector(col2rgb(defcols[ic])), opacity), collapse=','),')')
			p <- plotly::add_trace(p, data=df, x = ~x, y = ~y, mode = 'lines', fill = 'tozeroy', fillcolor = fcol, showlegend = FALSE)
		}
	}

	# Add Tags based on selection
	if (tags=='auto')
		tags <- ifelse( (ppmview[2]-ppmview[1])>1.5, 'peak', 'name' )

	if (tags == 'allpeaks' || (tags %in% c('id','name','peak') && !is.null(PL) && nrow(PL)>0)) {

		if (tags %in% c('id','peak','name')) {
			M <- NULL
			for (k in 1:nrow(PL)) {
				pid <- rownames(PL)[k]
				V <- simplify2array(lapply(1:nrow(cmpdlist),
						function(i) { which(pid %in% as.numeric(unlist(strsplit(cmpdlist[i,2],",")))) }))
				M <- rbind(M, c(PL[k,2], 1.05*PL[k,3], cmpdlist[which(as.numeric(V)==1),1], pid))
			}
			data <- data.frame(lab=M[,3], x=M[,1], y=M[,2], pkid=M[,4],
						tags=sapply(1:nrow(M), function(k) {which(unique(M[,3]) == M[k,3])}))
		}

		if (tags=='allpeaks') {
			PL <- spec$fit$peaks[ spec$fit$peaks$ppm>=ppmview[1] & spec$fit$peaks$ppm<=ppmview[2], , drop=F]
			data <- data.frame(x=PL[,2], y=1.05*PL[,3], pkid=rownames(PL))
		}

		if (tags %in% c('peak','allpeaks'))
			p <- p |> plotly::add_annotations(x = data$x, y = data$y, text = as.character(data$pkid),
				showarrow = TRUE, arrowcolor='red', hovertext=data$lab,
				hoverlabel=list(font = list(color = 'blue', family = 'sans serif', size = 18)),
				font = list(color = 'black', family = 'sans serif', size = 12))

		if (tags=='id')
				p <- p |> plotly::add_annotations(x = data$x, y = data$y, text = as.character(data$tags),
				showarrow = TRUE, arrowcolor='red', hovertext=data$lab,
				hoverlabel=list(font = list(color = 'blue', family = 'sans serif', size = 18)),
				font = list(color = 'black', family = 'sans serif', size = 12))

		if (tags=='name') {
			L <- lapply(unique(M[,3]), function(m){
				M0 <- M[ M[,3]==m, , drop=F]; 
				M0[which(as.numeric(M0[,2])==max(as.numeric(M0[,2]))), , drop=F]
			})
			M2 <- NULL
			for(k in 1:length(L)) M2 <- rbind(M2, L[[k]])
			data <- data.frame(lab=M2[,3], x=M2[,1], y=M2[,2])
			p <- p |> plotly::add_annotations(x = data$x, y = data$y, text = as.character(data$lab),
				showarrow = TRUE, arrowcolor='red', textangle=-30,
				font = list(color = 'black', family = 'sans serif', size = 12))
		}

	}

	# legend position
	if (legendhoriz)
		p <- p |> plotly::layout(legend = list(orientation = 'h', xanchor = "center",
				x = 0.5, y=ifelse(legendtop, 1.1, -0.2)))
	# Show / hide grid
	if (!showgrid)
		p <- p |> plotly::layout(xaxis = list(showgrid = F), yaxis = list(showgrid = F))

	# Final plot
	p <- p |> plotly::layout(colorway = arrColors, showlegend=showlegend)
	p
})


#=====================================================================
# Wrappers for plotting functions
#=====================================================================

internalClass$set("private", "getPlotParams", function(params=NULL)
{
	g <- def_plot_settings
	if ("list" %in% class(params))
		for (p in ls(params)) g[[p]] <- params[[p]]
	g
})

internalClass$set("private", "plotFinish", function(fig, params=NULL)
{
	g <- getPlotParams(params)

	# Remove Y-axis
	axcfg <- list(title="", zeroline=FALSE, showline=FALSE, showticklabels=FALSE, showgrid=FALSE)
	if (!g$yaxis)
        fig <- fig |> plotly::layout(yaxis = axcfg, font=g$font)
	else
        fig <- fig |> plotly::layout(yaxis = list(title=g$ylabel, rangemode="nonnegative", tickfont=g$font2, showgrid=g$showgrid))

	fig <- fig |> plotly::layout(xaxis = list(title=g$xlabel, tickfont=g$font2, showgrid=g$showgrid), font=g$font)
	if (!is.null(g$ppmview))
        fig <- fig |> plotly::layout(xaxis = list(autorange=FALSE, range=g$ppmview))
	if (!is.null(g$yview))
        fig <- fig |> plotly::layout(yaxis = list(autorange=FALSE, range=g$yview))

	OUTTYPE <<- g$mode
	displayWidget(fig, tmpdir=g$IMGDIR, width=g$width, height=g$height, filename=g$filename)
})

internalClass$set("public", "plotZones", function(ID, zones, params=NULL)
{
	g <- getPlotParams(params)

	zones <- c(min(zones):max(zones))
	if (sum(zones %in% res$zones) != length(zones))
		stop('Error: some zones are not included in the spectra set',"\n")

	S <- ifelse( is.numeric(ID), specList[[ID]]$samplecode, ID )
	if (is.null(g$filename)) {
		if (length(zones)>1)
			g$filename <- paste(S, min(zones), max(zones), sep='_')
		else
			g$filename <- paste(S, zones, sep='_')
	}

	pkfit <- PROFILE$fitting
	pkfit <- pkfit[pkfit$zone %in% zones, ]
	ppm_range <<- c( min(pkfit$ppm1), max(pkfit$ppm2))
	if (g$verbose) cat('ppm range =',min(ppm_range),'-',max(ppm_range),"\n")

	fig <- view_spectra(ID, plotmodel=g$plotmodel, plotTrueSpec=g$plotTrueSpec,
				plotresidus=g$plotresidus, plotzones=g$plotzones,
				tags=g$tags, lw=g$lw, showlegend=g$showlegend,
				legendhoriz=g$legendhoriz,  legendtop=g$legendtop, showgrid=g$showgrid, 
				title=g$title, colspecs=g$colspecs, colcpmds=g$colcpmds,
				opacity=g$opacity, verbose=g$verbose)
    plotFinish(fig, g)
})

internalClass$set("public", "plotCmpds", function(ID, cpmd, params=NULL)
{
	g <- getPlotParams(params)

	S <- ifelse( is.numeric(ID), specList[[ID]]$samplecode, ID )
	if (is.null(g$filename)) {
		if (length(cpmd)>1)
			g$filename <- paste(S, paste(sapply(cpmd, function(s) { substring(s,1,1) }), collapse='_'), sep='_')
		else
			g$filename <- paste(S, cpmd, sep='_')
	}

	fig <- plot_spectra(ID, cpmd, plotmodel=g$plotmodel, plotTrueSpec=g$plotTrueSpec, 
						plotresidus=g$plotresidus, plotzones=g$plotzones,
						tags=g$tags, lw=g$lw, showlegend=g$showlegend, legendhoriz=g$legendhoriz, legendtop=g$legendtop,
						showgrid=g$showgrid, title=g$title, colspecs=g$colspecs, colcpmds=g$colcpmds,
						verbose=g$verbose)
    plotFinish(fig, g)
})

