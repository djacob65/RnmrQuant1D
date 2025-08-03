#=====================================================================
# Plotly / DT outputs
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
	tryCatch(DT::saveWidget(widget,file=file, ...), error=function(e){})
	setwd(wd);
})

# Note : installation of orca - see https://github.com/plotly/orca#installation
#        need to have the processx R package installed as well.
# For ubuntu :
# cd /tmp
# wget https://github.com/plotly/orca/releases/download/v1.3.1/orca-1.3.1.AppImage
# chmod +x orca-1.3.1.AppImage
# mv orca-1.3.1.AppImage /usr/bin/
# apt-get install -y xvfb
# echo '#!/bin/bash' > /usr/bin/orca
# echo 'xvfb-run -a /usr/bin/orca-1.3.1.AppImage "$@"' >> /usr/bin/orca
# chmod +x /usr/bin/orca
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
		IMGfile <- paste0(IMGname, ".",type)
		tryCatch({
			suppressWarnings(suppressMessages(
				plotly::orca(widget, file=IMGfile, width=width, height=height, verbose=TRUE, debug=TRUE)))
		}, error=function(cond) { cat("Failed:", paste0(cond, collapse="\n"), "\n");  } )
		if (!file.exists(IMGfile)) IMGfile <- paste0(IMGname, "_1.",type)
		repeat {
			if (! file.exists(IMGfile)) break
			# Jupyter : return the figure
			if (.Platform$GUI != "RStudio" & !interactive()) {
				IRdisplay::display_html(paste0('<style type="text/css">div.output_',type,', img { max-width: 100%; height: ',height,'; }</style>'))
				if (type == "png") IRdisplay::display_png(file=IMGfile)
				if (type == "jpeg") IRdisplay::display_jpeg(file=IMGfile)
				if (type == "svg") IRdisplay::display_svg(file=IMGfile)
				break
			}
			# RStudio or R GUI
			if (type=='svg' && (.Platform$GUI == "RStudio" || interactive())) {
				if ('svgtools' %in% installed.packages()) {
					suppressWarnings(suppressMessages(library(svgtools)))
					suppressWarnings(suppressMessages(svgtools::display_svg(readLines(IMGfile))))
				}
			}
			if (type %in% c("png", "jpeg") && (.Platform$GUI == "RStudio" || interactive())) {
				if ('magick' %in% installed.packages()) {
					suppressWarnings(suppressMessages(library(magick)))
					print(magick::image_read(IMGfile), info=FALSE)
				}
			}
			# RStudio, R GUI or Terminal : return the file name
			setwd(wd);
			return(IMGfile)
			break
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


internalClass$set("public", "displayTable", function(M, nbdec=2, tmpdir='tmp', container=NULL)
{
	type <- OUTTYPE
	options(DT.options = list(pageLength = nrow(M)))
	optstyle <- DT::JS(
		"function(settings, json){$(this.api().table().header()).css({'font-size':'12px', 'background-color':'#c2d1f0', 'color':'#000'});}"
	)
	df=as.data.frame(M)
    if (is.null(container))
        container <- htmltools::tags$table(DT::tableHeader(c('ID',colnames(df)), TRUE), class = 'display')
	V <- sapply(df, is.numeric)
	for (k in 1:length(V)) if (V[k]) df[names(V)[k]] <- round(df[names(V)[k]],nbdec)
	m <- DT::datatable(df, options = list(dom = 't', initComplete = optstyle), container=container) |>
			DT::formatStyle( names(df), `font-size` = '12px') |>
			DT::formatStyle( 0, `font-size` = '12px')
	if (type == "html") {
		displayWidget(m)
	} else {
		IMGname <- floor(runif(1, 0, 10^12))
		IMGhtml <- file.path(tmpdir,paste0(IMGname, ".html"))
		saveWidgetFix(m, IMGhtml)
		IMGfile <- file.path(tmpdir,paste0(IMGname, ".",type))
		tryCatch(webshot2::webshot(IMGhtml, IMGfile), error=function(e){})
		if (type == "png") IRdisplay::display_png(file=IMGfile)
		if (type == "jpeg") IRdisplay::display_jpeg(file=IMGfile)
		if (type == "svg") IRdisplay::display_svg(file=IMGfile)
	}
})
