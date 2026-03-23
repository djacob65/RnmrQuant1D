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
		filters      = list(),
		filtersets   = list(),
		PPM_NOISE    = c(10.2,10.5), # ppm range corresponding to noise, and this, for all spectra
		TSPwidthMax  = 1,
		oblset       = 0,

	# Graphical parameters
		OUTTYPE      = "html",
		ppm_range    = c(0.5,10),

	# Instrument Field
		FIELD        = "unknown",

	# Instrument sequence
		SEQUENCE     = "noesy",

	# Sample Type
		TYPE         = "unknown",

	# Name of the dilution factor field
		FDILfield   = "F_dilution",

	# QC-QS type names
		QCtype       = "QC",
		QStype       = "QS",

	# Calibration compound names
		calib_cmpd_names  = c('TSP','TMSP','DSS'),

	# Default settings (P3) for metabolite identification rules/patterns
		patterns_defpars = list(),

	# Results proc-process
		specList     = list(),
		res          = list(),
		quantpars    = list(),
		fP           = list(),
		expnolist    = "",
 
	# Internally initialized sample list in RAWDIR directory
		RAWDIR_SLIST     = NULL,

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
			procParams$TSP <<- FALSE
			procParams$ADJPZTSP <<- TRUE
			procParams$MVPZTSP <<- TRUE
			procParams$MVPZFAC <<- 10
			procParams$DPHCPZTSP <<- 1.2
			procParams$DHZPZRANGE <<- 323

			# Filters
			filters1 <- list(main=c('daub8'), others=c('symlet8', 'smooth1'))
			filters2 <- list(main=c('smooth0'), others=c('smooth1'))
			filters3 <- list(main=c('smooth1'), others=c('smooth2', 'smooth3'))
			filtersets <<- list(filters1,filters2,filters3)
			filters  <<- filters1

			# Deconvolution Parameters
			opars <<- list(ratioPN=1, oneblk=1, pvoigt=1, oeta=1, oasym=1, asymmax=250, spcv=0.001, d2cv=0.05,
						lowPeaks=1, distPeaks=0.5, addPeaks=1, sndpass=0,  sigma_min=0.00025, addPsigmin=0.0005,
						R2limit=0.995)

			# Results proc-process
			res <<- list(allquantifs=NULL, peaklist=NULL, infos=NULL, proctype='')
			quantpars <<- list(cmpdlist=NULL, zones=NULL, tottime=0)
			specList <<- list()

			# Default settings (P3) for metabolite identification rules/patterns
			# See find_compounds()
			patterns_defpars <<- list(
				'b'  = c(0, 5),                  # nbpeaks, ratioPN
				's'  = c(1),                     # rank
				'd'  = c(0, 1.4, 0.3, 65),       # criterion, ratioA, dJ, dS
				't'  = c(0, 1.3, 0.35),          # criterion, ratioA, dJ
				'q'  = c(0, 3, 0.3),             # criterion, ratioA, dJ
				'dd' = c(0, 2.2),                # criterion, ratioA
				'm'  = c(0, 1.1, 0.3),           # criterion, ratioA, dJ
				'm2' = c(0, 1.1, 0.3),           # criterion, ratioA, dJ
				'r1' = c(1, 25),                 # rank, ratioPN
				'r2' = c(2, 4, 25),              # Jmin, Jmax, ratioPN
				'r3' = c(1.4),                   # ratioA
				'r4' = c(3),                     # dJ
				'r5' = c(0, 1.4, 0.3, 65),       # criterion, ratioA, dJ, dS
				'r6' = c(2, 2),                  # nbpeaks, dist
				'r7' = c(2, 1.2),                # ratioA, dist
				'r8' = c(2, 5, 5)                # J, nbpeaks, snrthres
			)
		}
	)

	# -- Private Fields --
	# private = list(	)

# -- End Public Class --
)


