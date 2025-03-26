#' RnmrQuant1D
#'
#' @importFrom R6 R6Class
#'
#' @author Daniel Jacob - (C) INRAE 2025
#'
#' @description
#' Package dedicated to 1D proton NMR quantification, including peak fitting and based on external calibration using standard spectra.
#'
#' @details
#' This package was initially developed as part of an ANR project on wine authenticity (\href{https://anr.fr/Project-ANR-21-CE21-0014}{ANR-21-CE21-0014}). However, it is generic enough to be used on other biological and/or food matrices. This involves the implementation of an analytical protocol allowing quantification from an external standard (see references).\cr\cr
#' This package provides all the methods needed: \cr 
#' \enumerate{
#'    \item to calculate different entities such as i) integration of patterns related to the targeted compounds of interest, ii) signal-to-noise ratio of motifs, iii) absolute quantifications of the targeted compounds,
#'    \item to calculate PULCON factors and verify using a quality control sample,
#'    \item  to exploit the results i) by producing a complete workbook gathering all the calculated data but also information needed for tracing, ii) by visualizing the spectra zones related to each targeted compound including the result of the deconvolution and the identification of patterns.
#' }
#' \cr
#' To work, this package needs a number of input files:\cr
#' \enumerate{
#'    \item a set of spectra concerning the samples. There may be for the same sample several repetitions acquired according to one or more sequences (zgpr, noesy)
#'    \item a set of spectra concerning quality control (QC) and quantification standards (QS). Similarly, there may be several repetitions acquired according to one or more sequences (zgpr, noesy)
#'    \item a file describing the samples by linking them to spectral data
#'    \item a file describing the QC and the QS samples in terms of compounds and their concentration, called "calibration profile"
#'    \item a file containing the parameters for deconvolution and quantification, zone by zone, called "quantification profile"
#' }
#' \cr
#' The user can refer to the \href{https://github.com/djacob65/RnmrQuant1D/wiki/}{online tutorial} for more details on the structure of the input files, and on the different steps to perform the quantifications, from data preparation to final results.\cr
#' \cr 
#' It was implemented based on the object-oriented programming using the R6 class.
#'
#'
#' @references
#'
#' Guillaume Leleu et al. (2024) Development of a standardized method for metabolite analysis by NMR to assess wine authenticity, IVES Conference Series, OIV, DOI:10.58233/IMGGSIME
#'
#' Teipel et al. (2020) Application of 1H Nuclear Magnetic Resonance Spectroscopy as Spirit Drinks Screener for Quality and Authenticity Control, Food 9, 1355; doi:10.3390/foods9101355
#'
#' @export
RnmrQuant1D <- R6Class("RnmrQuant1D",
	inherit = internalClass,
	portable = FALSE,
	cloneable = FALSE,

# -- Public Fields --
public = list(

#' @field RAWDIR Path of the sample spectra directory (required)
	RAWDIR       = "",

#' @field QSDIR Path of the QC-QS spectra directory (required)
	QSDIR        = "",

#' @field TMPDIR Output path for logs (required)
	TMPDIR       = "tmp",

#' @field RDATADIR Output path for RData files (required)
	RDATADIR     = "tmp",

#' @field SAMPLES The sample table (required)
	SAMPLES      = data.frame(),

#' @field PROFILE The quantification profile obtained by the \href{#method-RnmrQuant1D-readProfile}{\code{readProfile()}} method (required).
	PROFILE      = list(),

#' @field CALIBRATION The calibration profile as a data.frame (required)
	CALIBRATION  = data.frame(),

#' @field SEQUENCE The instrument sequence (required)
	SEQUENCE     = "unknown",

#' @field FIELD The instrument Field (recommended for tracing)
	FIELD        = "unknown",

#' @field TYPE The sample type (recommended for tracing)
	TYPE         = "unknown",

#' @field OUTTYPE specifies the output type for plotly graphs and data tables. The value can be "html", "png", "jpeg", or "svg" (default = "html")
	OUTTYPE      = "html",

#' @field FDILfield Name of the dilution factor column in the sample table (default = "F_dilution")
	FDILfield   = "F_dilution",

#' @field QCtype Identifier of the Quality Control (QC) in the calibration profile (1st column - default = "QC")
	QCtype       = "QC",

#' @field QStype Identifier of the Quantification Standards (QS) in the calibration profile (1st column - default = "QS")
	QStype       = "QS",

#' @description Class Constructor - Initializes internal parameters for preprocessing and deconvolution.
#' @param self the RnmrQuant1D instance
#' @return a RnmrQuant1D instance
	initialize = function() {
		super$initialize()
	},

#' @description
#' Checks if the sample table format is right and if sample names correspond to spectra under the RAWDIR directory
#' @param self The RnmrQuant1D instance 
#' @param verbose Display the success message if TRUE
#' @return nothing
	check_samples = function(verbose=FALSE) {
		super$check_samples(verbose)
	},

#' @description
#' Checks if the calibration profile is OK
#' @param self The RnmrQuant1D instance 
#' @param QC If specified, this is the name of the spectra directory for quality control (QC) (default is NULL)
#' @param QS If specified, this is the name of the spectra directory for qualification standards (QS) (default is NULL)
#' @param sequence If specified, this is the chosen sequence for testing. By default (NULL), the sequence in the RnmrQuant1D instance will be chosen.
#' @param verbose Display the success message if TRUE
#' @return nothing
	check_calibration = function(QC=NULL, QS=NULL, sequence=NULL, verbose=FALSE) {
		super$check_calibration(QC, QS, sequence, verbose)
	},

#' @description
#' Checks if the quantification profile is OK
#' @param self The RnmrQuant1D instance 
#' @param zones specifies the zone(s) to be checked
#' @param verbose Display the success message if TRUE
#' @return nothing
	check_profile = function(zones=NULL, verbose=FALSE) {
		super$check_profile(zones, verbose)
	},

#' @description
#' Checks if output directories (TMPDIR & RDATADIR) are OK
#' @param self The RnmrQuant1D instance 
#' @param verbose Display the success message if TRUE
#' @return nothing
	check_outdir = function(verbose=FALSE) {
		super$check_outdir(verbose)
	},

#' @description
#' Checks many things : directories, profiles (quantification & calibration), ...
#' @param self The RnmrQuant1D instance 
#' @param verbose Display the success message if TRUE
#' @return nothing
	check_all = function(verbose=FALSE) {
		super$check_all(verbose)
	},

#' @description
#' Generates the sample table for the corresponding sequence. 
#' @param self The RnmrQuant1D instance 
#' @param sequence If specified, take samples defined for this NMR sequence. Otherwise, takes samples for the one defined in the instance.
#' @param infos If TRUE, add some additional columns such as the TMSP width at half height and other acquisition parameters. 
#' @return the sample table
	get_samples_table = function(sequence=NULL, infos=FALSE) {
		super$get_samples_table(sequence, infos)
	},

#' @description
#' Reads the quantification profile. See the \href{https://github.com/djacob65/RnmrQuant1D/wiki/}{wiki page} for having more details.
#' @param self The RnmrQuant1D instance 
#' @param PROFILE the qualification profile filename
#' @return the qualification profile object
	readProfile = function(PROFILE) {
		super$readProfile(PROFILE)
	},

#' @description
#' Calculates the response factors according to the sample type ('QC' or 'QS'). 
#' @param self The RnmrQuant1D instance 
#' @param sampletype the sample type : 'QC' for Quality Control, or 'QS' for Quantification Standards.
#' @param samplename a valid sample name under the spectra directory (QSDIR)
#' @param thresfP defines the threshold of the CV (in percentage) of the response factor below which the spectrum (repetition) will not be taken into consideration. 
#' @param deconv defines whether the integration is calculated by deconvolution or by the Simpson approach.
#' @param verbose if TRUE, some messages are displayed
#' @return a response factor object, i.e a list with some attributes: fPUL, fP, fR, MC, INTG, fK
	get_response_factors = function(sampletype,  samplename, thresfP=5, deconv=TRUE, verbose=1) {
		super$get_response_factors(sampletype,  samplename, thresfP, deconv, verbose)
	},

#' @description
#' Returns the factor table according to the sample type ('QC' or 'QS'). 
#' @param self The RnmrQuant1D instance 
#' @param QS a response factor object returned by \href{#method-RnmrQuant1D-get_response_factors}{\code{get_response_factors()}}.
#' @return the factor table (matrix)
	get_factor_table = function(QS) {
		super$get_factor_table(QS)
	},

#' @description
#' Estimation of the compound concentrations (QC) based on the PULCON factor (QS)
#' @param self The RnmrQuant1D instance 
#' @param QC a response factor object for a QC type returned by \href{#method-RnmrQuant1D-get_response_factors}{\code{get_response_factors()}}.
#' @param QS a response factor object for a QS type returned by \href{#method-RnmrQuant1D-get_response_factors}{\code{get_response_factors()}}.
#' @param merge If TRUE, merge compounds that have the same names but with a different number at end, e.g. "citrate1", "citrate2" will be merged under "citrate".
#' @param verbose If TRUE, some messages are displayed
#' @return the estimation table (matrix)
	get_QC_estimation = function(QC, QS, merge=TRUE, verbose=1) {
		super$get_QC_estimation(QC, QS, merge, verbose)
	},

#' @description
#' Plots the estimation of the compound concentrations (QC) based on the PULCON factor (QS)
#' @param self The RnmrQuant1D instance 
#' @param QCest the estimation table returned by \href{#method-RnmrQuant1D-get_QC_estimation}{\code{get_QC_estimation()}}.
#' @return a \href{https://plotly.com/r/}{plotly} graph
	plot_QC_estimation = function(QCest) {
		super$plot_QC_estimation(QCest)
	},

#' @description
#' Performs the integration of all compounds defined by a pattern included in the specified zone(s).
#' @param self The RnmrQuant1D instance 
#' @param zones The selected zones for integration.
#' @param ncpu The number of processors (cores) used for parallel computing.
#' @param verbose if TRUE, some messages are displayed
#' @return The result will be stored in the RnmrQuant1D instance.
	proc_Integrals = function(zones, ncpu=2, verbose=1) {
		super$proc_Integrals(zones, ncpu, verbose)
	},

#' @description
#' Gets the integration matrix from the RnmrQuant1D instance, computed by \href{#method-RnmrQuant1D-proc_Integrals}{\code{proc_Integrals()}}.
#' @param self The RnmrQuant1D instance 
#' @return the integration matrix
	get_Matrix_Integrals = function() {
		super$get_Matrix_Integrals()
	},

#' @description
#' Gets the SNR matrix from the RnmrQuant1D instance, computed by \href{#method-RnmrQuant1D-proc_Integrals}{\code{proc_Integrals()}}.
#' @param self The RnmrQuant1D instance 
#' @return the SNR matrix
	get_Matrix_SNR = function() {
		super$get_Matrix_SNR()
	},

#' @description
#' Gets the CV matrix from the integration matrix if samples have repetitions.
#' @param self The RnmrQuant1D instance
#' @param MatInt The integration matrix. If NULL, it will be fetched from the the RnmrQuant1D instance. See \href{#method-RnmrQuant1D-get_Matrix_Integrals}{\code{get_Matrix_Integrals()}}.
#' @param nbrep the number of repetition. default=3.
#' @return the CV matrix
	get_Matrix_CV = function(MatInt=NULL, nbrep=3) {
		super$get_Matrix_CV(MatInt, nbrep)
	},

#' @description
#' Saves the result tables (samples, integration, SNR) computed by \href{#method-RnmrQuant1D-proc_Integrals}{\code{proc_Integrals()}} and other information for tracing (quantification profile, R environment, ...) into a workbook (\href{https://docs.fileformat.com/spreadsheet/xlsx/}{XLSX format}).
#' @param self The RnmrQuant1D instance
#' @param file the path to the output file
#' @param filelist a simple list containing two elements, just for tracing : 1) SAMPLEFILE, the path to the sample file, 2) PROFILE, the path to the quantification profile.
#' @return nothing
	save_Matrices = function(file, filelist=NULL) {
		super$save_Matrices(file, filelist)
	},

#' @description
#' Plots the ppm range defined by 'RnmrQuant1D$ppm_range', the latter being positioned by the zone(s) chosen when calculating the integrations by \href{#method-RnmrQuant1D-proc_Integrals}{\code{proc_Integrals()}}.
#' @param self The RnmrQuant1D instance
#' @param id the order number of the spectrum in the sample table.
#' @param plotmodel If TRUE, plot the model based on the deconvolution.
#' @param plotTrueSpec  If TRUE, plot the original spectrum without local baseline corrections.
#' @param plotzones If TRUE, add a semi-transparent rectangle superimposed on each fit zone.
#' @param tags add an arrow annotation on each peak resulting from the deconvolution. Among the possible values: 1) 'id' adds the number of the compound to which the peak belongs, 2) 'name' adds the name of the compound to which the peak belongs, 'auto' adds either the number if multiple zones, or the name of the compound, 3) 'none' adds no arrow annotation.
#' @param legendhoriz Put the legend at the bottom (TRUE) or at the left (FALSE) of the graph
#' @param verbose If TRUE, The peak list is displayed
#' @return a \href{https://plotly.com/r/}{plotly} graph
	view_spectra = function(id, plotmodel=TRUE, plotTrueSpec=TRUE, plotzones=TRUE, tags='none', legendhoriz=FALSE, verbose=FALSE) {
		super$view_spectra(id, plotmodel, plotTrueSpec, plotzones, tags, legendhoriz, verbose)
	},

#' @description
#' Computes the response factor for the QS sample.
#' @param self The RnmrQuant1D instance 
#' @param QSname a valid QS sample name under the spectra directory (QSDIR)
#' @param thresfP defines the threshold of the CV (in percentage) of the response factor below which the spectrum (repetition) will not be taken into consideration. 
#' @param deconv defines whether the integration is calculated by deconvolution or by the Simpson approach.
#' @param verbose if TRUE, some messages are displayed
#' @return he resulting object (fP) will be stored in the RnmrQuant1D instance.
	proc_fPULCON = function(QSname, thresfP=5, deconv=TRUE, verbose=1) {
		super$proc_fPULCON(QSname, thresfP, deconv, verbose)
	},

#' @description
#' Computes for each sample the concentration  of targeted compounds  based on the quantification profile.
#' @param self The RnmrQuant1D instance 
#' @param cmpdlist Sets the list of targeted compounds. It can be null and in this case, it will be the list of zones that will be taken into consideration.
#' @param zones In case the list of targeted compounds is null, then the targeted compounds will be defined by the patterns included in the indicated zones. If the list of zones is null, as well as that of the targeted compounds (cmpdlist), then takes into consideration all the zones defined in the quantification profile.
#' @param ncpu The number of processors (cores) used for parallel computing.
#' @param reset Indicates whether to delete all existing RData files or to proceed by accumulation. In this latter case, if some samples had already been processed then they will not be processed again.
#' @param CR Indicates whether to add a Carriage Return for a more presentable display (case with Jupyter notebooks)
#' @param verbose if TRUE, some messages are displayed
#' @return The results will be stored under the \code{"RDATADIR"} directory as RData files.
	proc_Quantification = function(cmpdlist=NULL, zones=NULL, ncpu=2, reset=FALSE, CR=FALSE, verbose=1) {
		super$proc_Quantification(cmpdlist, zones, ncpu, reset, CR, verbose)
	},

#' @description
#' Retrieves from RData files on disk generated by \href{#method-RnmrQuant1D-proc_Quantification}{\code{proc_Quantification()}}, the result tables within a list.
#' @param self The RnmrQuant1D instance 
#' @return the result tables within a list
	get_output_results = function() {
		super$get_output_results()
	},

#' @description
#' Retrieves from RData files generated on disk by \href{#method-RnmrQuant1D-proc_Quantification}{\code{proc_Quantification()}}, the spectra and the peaklist within lists in the RnmrQuant1D instance respectively named specList and res.
#' @param self The RnmrQuant1D instance 
#' @return nothing
	get_spectra_data = function() {
		super$get_spectra_data()
	},

#' @description
#' Saves the result tables (samples, integration, SNR, quantifications) computed by \href{#method-RnmrQuant1D-proc_Quantification}{\code{proc_Quantification()}} and other information for tracing (calibration profile, quantification profile, R environment, ...) into a workbook (\href{https://docs.fileformat.com/spreadsheet/xlsx/}{XLSX format}).
#' @param self The RnmrQuant1D instance
#' @param file the path to the output file
#' @param filelist a simple list containing three elements, just for tracing : 1) SAMPLEFILE, the path to the sample file, 2) PROFILE, the path to the quantification profile, 3) CALIBRATION, the path to the calibration profile
#' @return nothing
	save_Results = function(file, filelist=NULL) {
		super$save_Results(file, filelist)
	},

#' @description
#' Visualizes the compound 'compound' for spectrum 'id', where id is the order number of the spectrum in the sample table. 
#' @param id The order number of the spectrum in the sample table
#' @param compound The compound name taken from quantpars$cmpdlist in the RnmrQuant1D instance.
#' @param ... Arguments to be passed to the \href{#method-RnmrQuant1D-view_spectra}{\code{view_spectra()}} method.
#' @return a \href{https://plotly.com/r/}{plotly} graph
	plot_spectra = function(id, compound, ...) {
		super$plot_spectra(id, compound, ...)
	},

#' @description
#' Displays a widget (e.g. a plotly graph) - Depending on the OutType attribute in the RnmrQuant1D instance, returns a \href{https://plotly.com/r/}{plotly} graph (html) or displays the image (png,svg, ...)
#' @param widget a plotly graph
#' @param tmpdir the directory to temporarily save the image before displaying it (png,svg, ...)
#' @param width the width of the output image
#' @param height the height of the output image
#' @return depending on the OutType attribute in the RnmrQuant1D instance, returns a \href{https://plotly.com/r/}{plotly} graph (html) or display the image (png,svg, ...)
	displayWidget = function(widget, tmpdir='tmp', width='auto', height=400) {
		super$displayWidget(widget, tmpdir, width, height)
	},

#' @description
#' Saves a widget (e.g. a plotly graph) into a file - A wrapper to saveWidget which compensates for arguable BUG in htmlwidgets::saveWidget which requires `file` to be in current working directory.
#' @param widget a plotly graph
#' @param file the path of the output filename 
#' @param ... Arguments to be passed to htmlwidgets::saveWidget
#' @return nothing
	saveWidgetFix = function(widget,file, ...) {
		super$saveWidgetFix(widget,file, ...)
	},

#' @description
#' Beautifies a matrix using the DT package before displaying it
#' @param M the matrix to be embellished with the DT package
#' @param nbdec the number of decimal places for floating point numbers in the output table
#' @param tmpdir the directory to temporarily save the image before displaying it
#' @return Depending on the OutType attribute in the RnmrQuant1D instance, returns a \href{https://plotly.com/r/}{plotly} graph (html) or displays the image (png,svg, ...)
	displayTable = function(M, nbdec=2, tmpdir='tmp') {
		super$displayTable(M, nbdec, tmpdir)
	}

))


# --- Note ---
# To render the documentation in HTML format
# 1) devtools::document()
# 2) tools::Rd2HTML(file.path(getwd(),"man/RnmrQuant1D.Rd"), out=file.path(getwd(),"man/RnmrQuant1D.html"))
# ------------
