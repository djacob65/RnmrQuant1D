#=====================================================================
# Wrappers for plotting functions
#=====================================================================

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

	if (is.null(g$filename)) {
		if (length(zones)>1)
			g$filename <- paste(specList[[ID]]$samplecode, min(zones), max(zones), sep='_')
		else
			g$filename <- paste(specList[[ID]]$samplecode, zones, sep='_')
	}

	pkfit <- PROFILE$fitting
	pkfit <- pkfit[pkfit$zone %in% zones, ]
	ppm_range <<- c( min(pkfit$ppm1), max(pkfit$ppm2))
	if (g$verbose) cat('ppm range =',min(ppm_range),'-',max(ppm_range),"\n")

	fig <- view_spectra(ID, plotmodel=g$plotmodel, plotTrueSpec=g$plotTrueSpec,
				plotresidus=g$plotresidus, plotzones=g$plotzones,
				tags=g$tags, lw=g$lw, showlegend=g$showlegend,
				legendhoriz=g$legendhoriz, showgrid=g$showgrid, 
				title=g$title, colspecs=g$colspecs, colcpmds=g$colcpmds,
				opacity=g$opacity, verbose=g$verbose)
    plotFinish(fig, g)
})

internalClass$set("public", "plotCmpds", function(ID, cpmd, params=NULL)
{
	g <- getPlotParams(params)

	if (is.null(g$filename)) {
		if (length(cpmd)>1)
			g$filename <- paste(specList[[ID]]$samplecode,
						paste(sapply(cpmd, function(s) { substring(s,1,1) }), collapse='_'), sep='_')
		else
			g$filename <- paste(specList[[ID]]$samplecode, cpmd, sep='_')
	}

	fig <- plot_spectra(ID, cpmd, plotmodel=g$plotmodel, plotTrueSpec=g$plotTrueSpec, 
						plotresidus=g$plotresidus, plotzones=g$plotzones,
						tags=g$tags, lw=g$lw, showlegend=g$showlegend, legendhoriz=g$legendhoriz,
						showgrid=g$showgrid, title=g$title, colspecs=g$colspecs, colcpmds=g$colcpmds,
						verbose=g$verbose)
    plotFinish(fig, g)
})

