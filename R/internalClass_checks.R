# =====================================================================
# Some functions for checking directories / data required
# =====================================================================

internalClass$set("private", "check_sequence", function()
{
	# Check if the sequence is recognized
	if (! SEQUENCE %in% c('zg','zgpr', 'noesy'))
		stop_quietly(paste0("Error: ", SEQUENCE, "is not valid. Only 'zg','zgpr', 'noesy' are recognized."))
})

internalClass$set("public", "check_samples", function(verbose=FALSE)
{
	# Check if RAWDIR exist
	if (!dir.exists(RAWDIR)) 
		stop_quietly(paste0("Error: ", RAWDIR, " does not exist."))

	# Check if RAWDIR empty
	L1 <- get_list_samples(RAWDIR)
	if (length(L1)==0)
		stop_quietly(paste0("Error: ", RAWDIR, " is empty."))

	# Check if RAWDIR contains some bruker spectra
	M1 <- get_list_spectrum(RAWDIR, L1)
	if (is.null(M1))
		stop_quietly(paste0("Error: ", RAWDIR, " must contain some Bruker spectra."))

	# Check if the sequence is recognized
	check_sequence()

	# Check if RAWDIR contain some Bruker spectra acquired with the request sequence
	M1b <- M1[ M1[,3] == SEQUENCE, ,drop=F]
	if (nrow(M1b)==0)
		stop_quietly(paste0("Error: ", RAWDIR, " must contain some Bruker spectra acquired with a ",SEQUENCE," sequence"))
	if (verbose) cat(paste0("OK: ", RAWDIR, " contains some Bruker spectra acquired with a ",SEQUENCE," sequence\n"))

	# Check if som spectra in sample tables match a spectrum under RAWDIR
	M1c <- M1b[ M1b[,1] %in% SAMPLES[,1],, drop=F]
	if (nrow(M1c)==0)
		stop_quietly(paste0("Error: No spectrum declared in the sample table seems to match a spectrum under ", RAWDIR))

	# Check if there are more samples than spectra under RAWDIR
	if (nrow(SAMPLES)>nrow(M1c))
		stop_quietly(paste0("Error: there are more samples than spectra under ", RAWDIR))
	if (verbose) cat(paste0("OK: all spectra reported in the sample table appear to match a spectrum under ", RAWDIR, "\n"))

	# Check if the samples table has at least NCMIN columns and a 'F_dilition' columns
	NCMIN <- 4
	if (!ncol(SAMPLES)>=NCMIN)
		stop_quietly(paste0("Error: the sample table must have at least ",NCMIN," columns : Spectrum, Samplecode, ",FDILfield))

	# Check if the sample names in the sample table match the same expno as found previously
	V1 <- simplify2array(lapply(1:nrow(M1b), function(k){paste0(M1b[k,1],'-',M1b[k,2])}))
	V2 <- simplify2array(lapply(1:nrow(SAMPLES), function(k){paste0(SAMPLES[k,1],'-',SAMPLES[k,3])}))
	if (sum(V2 %in% V1) != nrow(SAMPLES))
		stop_quietly(paste0("Error: Some EXPNO numbers in the sample table  do not match the provided sequence"))

	# Check if the samples table has a 'F_dilition' columns
	if (! FDILfield %in% colnames(SAMPLES))
		stop_quietly(paste0("Error: the sample table must have a columns named '",FDILfield))
	if (verbose) cat(paste0("OK: the sample table has at least ",NCMIN," columns with one named '", FDILfield, "'\n"))
})

internalClass$set("public", "check_calibration", function(QC=NULL, QS=NULL, sequence=NULL, verbose=FALSE)
{
	if (is.null(sequence)) sequence <- SEQUENCE

	# Check if the sequence is recognized
	check_sequence()

 	# Check if QSDIR exist
	if (!dir.exists(QSDIR))
		stop_quietly(paste0("Error: ", QSDIR, " does not exist."))
	
	# Check if QSDIR empty
	L2 <- get_list_samples(QSDIR)
	if (length(L2)==0)
		stop_quietly(paste0("Error: ", QSDIR, " is empty."))
	
	# Check if QSDIR contains some bruker spectra
	M2 <- get_list_spectrum(QSDIR, L2)
	if (is.null(M2))
		stop_quietly(paste0("Error: ", QSDIR, " must contain some Bruker spectra."))

	# Check if RAWDIR contain some Bruker spectra acquired with the request sequence
	M2b <- M2[ M2[,3] == sequence, ,drop=F]
	if (nrow(M2b)==0)
		stop_quietly(paste0("Error: ", QSDIR, " must contain some Bruker spectra acquired with a ",sequence," sequence"))

	if (verbose) cat(paste0("OK: ", QSDIR, " contains some Bruker spectra acquired with a ",sequence," sequence\n"))

	# Check if the CALIBRATION table has a rigth column names
	if (sum(colnames(CALIBRATION) %in% c('Type','Compound','MW','NH','MC','PPM1','PPM2'))!=7)
		stop_quietly(paste0("Error: ", "the standard profile must contain 7 columns named : 'Type','Compound','MW','NH','MC','PPM1','PPM2'"))

	# Check if the Type column contains only 'QC' and 'QS'
	if (sum(CALIBRATION$Type %in% c(QCtype,QStype)) != nrow(CALIBRATION))
		stop_quietly(paste0("Error: ", "the Type column of the standard profile must contain only 'QC' and 'QS'"))

	if (verbose) cat(paste0("OK: the standard profile format seems correct\n"))

	# Check if QCname corresponds to some spectra under QSDIR
	if (!is.null(QC)) {
		M2c <- M2b[ M2b[,1]==QC, , drop=F]
		if (nrow(M2c)==0)
			stop_quietly(paste0("Error: No spectrum ",QC," seems to match a spectrum under ", QSDIR))
		if (verbose) cat(paste0("OK: some spectra under ",QSDIR," correspond to ",QC,"\n"))
	}

	# Check if QCname corresponds to some spectra under QSDIR
	if (!is.null(QS)) {
		M2c <- M2b[ M2b[,1]==QS, , drop=F]
		if (nrow(M2c)==0)
			stop_quietly(paste0("Error: No spectrum ",QS," seems to match a spectrum under ", QSDIR))
		if (verbose) cat(paste0("OK: some spectra under ",QSDIR," correspond to ",QS,"\n"))
	}

})

internalClass$set("public", "check_profile", function(zones=NULL, verbose=FALSE)
{
	# Note: We assume that the quantification profile can be read without error by the readProfile method

	# Check if PROFILE has the right sections : preprocess,  fiiting, quantif & compound
	if (is.null(PROFILE$preprocess))
		stop_quietly(paste0("Error: there no 'preprocess' section define in the quantification profile"))
	if (is.null(PROFILE$fitting))
		stop_quietly(paste0("Error: there no 'fitting' section define in the quantification profile"))
	if (is.null(PROFILE$quantif))
		stop_quietly(paste0("Error: there no 'quantif' section define in the quantification profile"))
	if (is.null(PROFILE$compound))
		stop_quietly(paste0("Error: there no 'compound' section define in the quantification profile"))

	# Check if all quantif zone match with a fitting zone
	L <- PROFILE$quantif$zone %in% PROFILE$fitting$zone
	if (length(L) != sum(L)) {
		bad <- paste(PROFILE$quantif$zone[which(L==FALSE)],sep=',')
		stop_quietly(paste0("Error: there some 'quantif' zones (",bad,") that do not match with a fitting zone in the quantification profile"))
	}

	# Check if the requested zones are defined in the profile
	if (!is.null(zones) && sum(zones %in% PROFILE$fitting$zone)<length(zones))
		stop_quietly(paste0("Error: The zone(s) '",paste0(zones,collapse=','),"'  do not correspond to those defined in the quantification profile"))

	# Check that no zone is included in another
	fit <- PROFILE$fitting
	L <- simplify2array(lapply(1:nrow(fit), function(k){sum(fit$ppm1[k]>fit$ppm1 & fit$ppm2[k]<fit$ppm2)}))
	if (sum(L)>0)
		stop_quietly(paste0("Error: The zone(s) '",paste0(which(L>0),collapse=','),"' is (are) included in another in the quantification profile"))

	# Check if all quantif compound match with a compound line
	L <-  PROFILE$quantif$compound %in% PROFILE$compound$name
	if (length(L) != sum(L)) {
		bad <- paste(PROFILE$quantif$compound[which(L==FALSE)],sep=',')
		stop_quietly(paste0("Error: there some 'quantif' compounds (",bad,") that do not match with a compound defined in the 'compound' section in the quantification profile"))
	}

	# Check if all ppm range in the fitting section are positive
	L <- which((PROFILE$fitting[2] - PROFILE$fitting[1])<0)
	if (length(L)) {
		bad <- paste(PROFILE$fitting$zone[L],sep=',')
		stop_quietly(paste0("Error: there some 'fitting' zones (",bad,") with wrong ppm range in the quantification 
		profile"))
	}

	# Check if all ppm range in the quantif section are positive
	if (!is.na(as.numeric(rq1d$FIELD))) {
		M <- get_quantif_ppmrange()
		L <- which((M[,2] - M[,1])<0)
		if (length(L)) {
			bad <- paste(PROFILE$quantif$compound[L],sep=',')
			stop_quietly(paste0("Error: there some 'quantif' compounds (",bad,") with wrong ppm range in the quantification profile"))
		}
	}

	# Check if all patterns/rules are recognized
	patterns <- c('b','s','d','t','q','dd','m','m2','r1','r2','r3','r4','r5','r6','r7','r8','r9','r10','r11')
	L<- unique(PROFILE$quantif$pattern) %in% patterns
	if (length(L) != sum(L)) {
		bad <- paste(unique(PROFILE$quantif$pattern)[which(L==FALSE)],sep=',')
		stop_quietly(paste0("Error: there some pattern in the 'quantif' lines (",bad,") that do not match with a recognized pattern in the quantification profile"))
	}

	if (verbose) cat(paste0("OK: the quantification profile format seems correct\n"))
})


internalClass$set("public", "check_outdir", function(verbose=FALSE)
{
	# Check if TMPDIR exist
	if (!dir.exists(TMPDIR))
		stop_quietly(paste0("Error: ", TMPDIR, " does not exist."))

	# Check if RDATADIR exist
	if (!dir.exists(RDATADIR))
		stop_quietly(paste0("Error: ", RDATADIR, " does not exist."))

	if (verbose) cat("OK: output directories are correctly defined\n")

})

internalClass$set("public", "check_all", function(verbose=FALSE)
{
	warn_cnt <- 0

	check_samples(verbose=verbose)
	check_calibration(verbose=verbose)
	check_profile(verbose=verbose)
	check_outdir(verbose=verbose)

	# Check if Field is defined
	if (is.na(as.numeric(FIELD))) {
		cat("Warning : the field 'FIELD' is not defined. Useful for tracing")
		warn_cnt <- warn_cnt + 1
	}

	# Check if Type is defined
	if (TYPE == "unknown") {
		cat("Warning : the field 'TYPE' is not defined. Useful for tracing")
		warn_cnt <- warn_cnt + 1
	}

	if (verbose) cat("Success: 0 errors, ",warn_cnt,"warnings\n")
})

