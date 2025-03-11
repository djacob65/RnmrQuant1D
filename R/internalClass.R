internalClass <- R6Class("internalClass",
	portable = FALSE,
	cloneable = FALSE,

	# -- Public Fields --
	public = list(

	# Path
		RAWDIR       = "",
		QSDIR        = "",
		TMPDIR       = "tmp",
		RDATADIR     = "tmp",

	# Sample table
		SAMPLES      = data.frame(),

	# Profiles
		PROFILE       = list(),
		CALIBRATION  = data.frame(),

	# Processing parameters
		procParams   = list(),
		opars        = list(),
		SCALE_INT    = 1e4,
		allFilterset = c('daub8', 'symlet8', 'smooth1'),
		filterset    = c('daub8'),
		PPM_NOISE    = c(10.2,10.5), # ppm range corresponding to noise, and this, for all spectra
		TSPwidthMax  = 1,
		oblset       = 0,

	# Graphical parameters
		OUTTYPE      = "html",
		ppm_range    = c(0.5,10),

	# Instrument Field
		FIELD        = "unknown",

	# Instrument sequence
		SEQUENCE     = "unknown",

	# Sample Type
		TYPE         = "unknown",

	# Name of the dilution factor field
		FDILfield   = "F_dilution",

	# QC-QS type names
		QCtype       = "QC",
		QStype       = "QS",

	# Calibration compound names
		calib_cmpd_names  = c('TSP','TMSP','DSS'),

	# Results proc-process
		specList     = list(),
		res          = list(),
		quantpars    = list(),
		fP           = list(),
		expnolist    = "",

	#=====================================================================
	# Initialization
	#=====================================================================

		initialize = function() {
			options(stringsAsFactors=FALSE)
			options(warn=-1)
			options(warnings=-1)

			# Init Plotly & DT
			p <- plotly::plot_ly(x=1:100, y=rnorm(100))
			dt <- DT::datatable(data.frame(x=1:100, y=rnorm(100)))
			OUTTYPE <<- "html" # png

			# Preprocessing parameters
			procParams <<- Rnmr1D::Spec1rProcpar
			procParams$DEBUG <<- TRUE
			procParams$LOGFILE <<- ""
			procParams$VENDOR <<- 'bruker'
			procParams$INPUT_SIGNAL <<- 'fid'
			procParams$LB <<- 0.25
			procParams$OPTPHC0 <<- TRUE
			procParams$OPTPHC1 <<- FALSE
			procParams$ZEROFILLING <<- TRUE
			procParams$ZFFAC <<- 4
			procParams$TSP <<- TRUE
			procParams$ADJPZTSP <<- TRUE
			procParams$MVPZTSP <<- TRUE
			procParams$MVPZFAC <<- 10
			procParams$DHZPZRANGE <<- 323

			# Deconvolution Parameters
			opars <<- list(ratioPN=1, oneblk=1, pvoigt=1, oeta=1, oasym=1, asymmax=250, spcv=0.001, d2cv=0.05,
						lowPeaks=1, distPeaks=0.5, addPeaks=1, sndpass=0,  sigma_min=0.00025, addPsigmin=0.0005,
						R2limit=0.995)

			# Results proc-process
			res <<- list(allquantifs=NULL, peaklist=NULL, infos=NULL, proctype='')
			quantpars <<- list(cmpdlist=NULL, zones=NULL, tottime=0)
			specList <<- list()

		}
	)

	# -- Private Fields --
	# private = list(	)

# -- End Public Class --
)
