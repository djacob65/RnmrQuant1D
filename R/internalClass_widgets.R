#=====================================================================
# Plots / DT
#=====================================================================

internalClass$set("public", "saveWidgetFix", function(widget,file, ...)
{
	## A wrapper to saveWidget which compensates for arguable BUG in
	## saveWidget which requires `file` to be in current working
	## directory.
	wd<-getwd()
	on.exit(setwd(wd))
	outDir<-dirname(file)
	file<-basename(file)
	setwd(outDir);
	DT::saveWidget(widget,file=file, ...)
	setwd(wd);
})


internalClass$set("public", "displayWidget", function(widget, tmpdir='tmp', width='auto', height=400)
{
	type <- OUTTYPE
	if (!is.null(widget) && type %in% c("png", "jpeg", "svg")) {
		wd<-getwd()
		on.exit(setwd(wd))
		IMGname <- floor(runif(1, 0, 10^12))
		if (is.null(tmpdir))
			tmpdir <- paste0("tmp", floor(runif(1, 0, 10^12)))
		suppressWarnings(dir.create(tmpdir))
		setwd(tmpdir)
		tryCatch({
			sink("error.log")
			suppressWarnings(plotly::orca(widget, file=paste0(IMGname, ".",type), width=width, height=height))
			sink()
		}, error=function(cond) { sink() } )
		IMGfile <- paste0(IMGname, ".",type)
		if (!file.exists(IMGfile)) IMGfile <- paste0(IMGname, "_1.",type)
		if (file.exists(IMGfile)) {
			IRdisplay::display_html(paste0('<style type="text/css">div.output_',type,', img { max-width: 100%; height: ',height,'; }</style>'))
			if (type == "png") IRdisplay::display_png(file=IMGfile)
			if (type == "jpeg") IRdisplay::display_jpeg(file=IMGfile)
			if (type == "svg") IRdisplay::display_svg(file=IMGfile)
		}
		setwd(wd);
		#unlink(tmpdir, recursive=TRUE, force = TRUE)
	}
	if (!is.null(widget) && type %in% c("html")) {
		# see https://plotly.com/r/configuration-options/
		plotly::as_widget(widget) |> plotly::config(displaylogo = FALSE, displayModeBar = TRUE, scrollZoom = TRUE,
					toImageButtonOptions = list(format= 'svg'),
					modeBarButtonsToRemove = c('zoom2d','pan2d', 'select2d', 'lasso2d', 'zoomIn2d', 'zoomOut2d',
											'hoverClosestCartesian', 'hoverCompareCartesian','resetScale'))
	}
})


internalClass$set("public", "displayTable", function(M, nbdec=2, tmpdir='tmp')
{
	type <- OUTTYPE
	options(DT.options = list(pageLength = nrow(M)))
	optstyle <- DT::JS(
		"function(settings, json){$(this.api().table().header()).css({'font-size':'12px', 'background-color':'#c2d1f0', 'color':'#000'});}"
	)
	df=as.data.frame(M)
	V <- sapply(df, is.numeric)
	for (k in 1:length(V)) if (V[k]) df[names(V)[k]] <- round(df[names(V)[k]],nbdec)
	m <- DT::datatable(df, options = list(dom = 't', initComplete = optstyle)) |>
			DT::formatStyle( names(df), `font-size` = '12px') |>
			DT::formatStyle( 0, `font-size` = '12px')
	if (type == "html") {
		displayWidget(m, type)
	} else {
		IMGname <- floor(runif(1, 0, 10^12))
		IMGhtml <- file.path(tmpdir,paste0(IMGname, ".html"))
		saveWidgetFix(m, IMGhtml)
		IMGfile <- file.path(tmpdir,paste0(IMGname, ".",type))
		webshot2::webshot(IMGhtml, IMGfile)
		if (type == "png") IRdisplay::display_png(file=IMGfile)
		if (type == "jpeg") IRdisplay::display_jpeg(file=IMGfile)
		if (type == "svg") IRdisplay::display_svg(file=IMGfile)
	}
})
