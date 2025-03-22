#=====================================================================
# Intermediate level functions for Processing / Deconvolution / Quantification
#=====================================================================

internalClass$set("private", "applyReadSpectrum", function(ACQDIR, verbose=1)
{
	# Retrieve preprocessing parameters using the profile file (PROFILE)
	procParams <<- get_procParams(PROFILE)

	# Read the NMR spectrum from the acquisition directory (ACQDIR)
	# - Uses the processing parameters (procParams)
	spec <- Rnmr1D::readSpectrum(ACQDIR, procParams, PPM_NOISE, NULL, SCALE_INT, verbose=verbose)

	# Return the processed spectrum
	spec
})
internalClass$set("private", "applyBLcorrection", function(spec, verbose=1)
{
	# Step 1: Zero out specific ppm zones defined in PROFILE$exclude_zones
	if( ! is.null(PROFILE$exclude_zones)) {
		if (verbose) cat("-------------------\n")
		for( k in 1:nrow(PROFILE$exclude_zones) ) {
			if (verbose) cat('Zeroing the ppm zone : [',paste(PROFILE$exclude_zones[k,],collapse=', '),"]\n")
			spec$int[getseq(spec,as.numeric(PROFILE$exclude_zones[k,]))] <- 0
		}
	}

	# Step 2: Zero out negative intensities in specific ppm ranges
	if( ! is.null(PROFILE$zeroneg)) {
		if (verbose) cat("-------------------\n")
		zz <- PROFILE$zeroneg
		for (k in 1:nrow(zz)) {
			if (verbose) cat('Zeroing negative intensities within the ppm zone : [',paste(c(zz$ppm1[k],zz$ppm2[k]),collapse=', '),"]\n")
			iseq <- getseq(spec, c(zz$ppm1[k],zz$ppm2[k]));
			V <- spec$int[iseq];
			V[V<0] <- 0
			spec$int[iseq] <- V
		}
	}

	# Step 3: Initialize the baseline correction array
	BL <- rep(0, length(spec$int))

	# Step 4: Apply baseline correction if PROFILE$baseline exists
	if( ! is.null(PROFILE$baseline)) {
		if (verbose) cat("-------------------\n")
		blpars <- PROFILE$baseline
		for (k in 1:nrow(blpars)) {
			if (verbose) cat('Baseline correction in the ppm zone : [',paste(c(blpars$ppm1[k],blpars$ppm2[k]),collapse=', '),
				'], lambda =',blpars$lambda[k],', order =',blpars$order[k],"\n")
			BLk <- getBLexternal(spec, blpars$lambda[k], blpars$order[k], c(blpars$ppm1[k],blpars$ppm2[k]), blpars$smooth[k])
			spec$int <- spec$int - BLk
			BL <- BL + BLk
		}
		# Apply an additional baseline correction offset based on the mean intensity in the noisy ppm range
		BL0 <- 0.25*mean(spec$int[getseq(spec,c(9.5,10))])
		spec$int <- spec$int - BL0
		BL <- BL + BL0
	}

	# Step 5: Store the corrected intensity and baseline in the spectrum object
	spec$intcorr <- spec$int
	spec$BL <- BL

	# Return the modified spectrum object
	spec
})

internalClass$set("private", "applyPeakFitting", function(spec, opars=NULL, zones=NULL, ncpu=2, verbose=1)
{
	# If no fitting parameters are provided, use the default ones stored in 'self$opars'
	if (is.null(opars)) opars <- self$opars

	# Initialize the variable to store the peak fitting result
	specfit <- NULL

	# Apply peak fitting depending on the number of available CPU
	if (ncpu>1) {
		# Use multi-threaded peak fitting for faster performance
		specfit <- applyPeakFitting2(spec, opars, zones, ncpu=ncpu, verbose=verbose)
	} else {
		# Use single-threaded peak fitting if only one CPU is available
		specfit <- applyPeakFitting1(spec, opars, zones, verbose=verbose)
	}
	specfit
})

internalClass$set("private", "applyPeakFitting1", function(spec, opars, zones=NULL, verbose=1)
{
	# If verbose mode is enabled, print the default parameters and peak fitting method
	if (verbose) {
		cat("-------------------\n")
		cat('Default Parameters: ratioPN =', opars$ratioPN,', asymetric =',opars$oasym,
		    ', lowPeaks =',opars$lowPeaks,', addPeaks =',opars$addPeaks,', sndpass =',opars$sndpass, "\n")
		cat('Peak fitting = Internal method',"\n")
		cat("-------------------\n")
	}

	# Extract peak fitting parameters from the profile
	pkfit <- PROFILE$fitting

	# If specific zones are provided, filter the peak fitting settings accordingly
	if (!is.null(zones)) pkfit <- pkfit[ pkfit$zone %in% zones, , drop=F]

	# If no valid zones exist, print an error message and return the spectrum unchanged
	if (nrow(pkfit) == 0) {
		cat("Error: There is no zone corresponding to those given as input")
		return(spec)
	}

	# Initialize variables
	Ymodel <- Y <- spec$int <- spec$intcorr  # Store intensity and corrected intensity values
	Peaks <- infos <- NULL  # To store peak fitting results and metadata
	BLSIG <- 10             # Baseline significance threshold
	ratioPN <- 5            # Ratio threshold for peak detection

	# Loop through each peak fitting range
	for (k in 1:nrow(pkfit)) {
		ppmrange <- c(pkfit[k,1], pkfit[k,2])  # Define PPM range
		if (verbose) cat("PPM range = [", pkfit[k,1],',',pkfit[k,2],"]","\n")

		# Extract and configure peak fitting parameters for this range
		opars.loc <- opars
		opars.loc$oasym <- ifelse(pkfit$asym[k] > 0, 1, 0)
		opars.loc$asymmax <- pkfit$asym[k]
		opars.loc$pvoigt <- ifelse(pkfit$etamin[k] < 1, 1, 0)
		opars.loc$etamin <- pkfit$etamin[k]
		opars.loc$addPeaks <- pkfit$addpeaks[k]
		opars.loc$qbl <- pkfit$qbl[k]

		# Select appropriate filter set based on the filter parameter
		filters <- filtersets[[min(length(filtersets), round(abs(pkfit$filters[k])))]]

		# Convert "obl" (optimized baseline correction i.e integrated into the model) values from string to integer
		if (grepl(',', pkfit$obl[k])) {
			obl <- strtoi(simplify2array(strsplit(pkfit$obl[k], ',')))
		} else {
			obl <- strtoi(pkfit$obl[k])
		}

		# If verbose mode is enabled, print the parameters and peak fitting method for this range
		if (verbose) {
			cat('Parameters: addPeaks =', opars.loc$addPeaks,', qbl =', opars.loc$qbl,', oblset =', pkfit$obl[k], "\n")
			cat('Parameters: asymmax =', opars.loc$asymmax,', etamin =', opars.loc$etamin, "\n")
			cat('filter =', paste(filters$main, collapse=','), "\n")
		}

		# If baseline correction (qbl) is enabled, perform baseline correction
		if (opars.loc$qbl != 0) {
			spec$int <- spec$intcorr
			if (opars.loc$qbl == 1) {
				BL <- qnmrbc(spec, ppmrange, BLSIG, 5) 
			} else { 
				BL <- getBLexternal(spec, opars.loc$qbl, 2, ppmrange) 
			}
			spec$int <- spec$int - BL
			Y <- Y - BL
		}

		# Perform peak fitting using the selected filters
		opars.save <- opars.loc  # Save parameters before running peak fitting
		t <- system.time({
			model <- Rnmr1D::LSDeconv(spec, ppmrange, opars.loc, filters$main, obl, verbose = verbose)
		})
		if (verbose) cat('elapsed time =', round(t[3], 2), ', Ended at ', format(Sys.time(), "%m/%d/%Y - %X"), "\n")

		# If the first fit doesn't meet the R² threshold, try alternative filters
		opars.loc <- opars.save
		if (is.null(model) || (!is.null(model) && model$R2 < opars.loc$R2limit)) {
			t <- system.time({
				opars.loc$peaks <- NULL
				others_filters <- filters$others[!filters$others %in% filters$main]
				model2 <- Rnmr1D::LSDeconv(spec, ppmrange, opars.loc, others_filters, obl, verbose = verbose)

				# Choose the best model (if model2 has a better R², use it)
				if ((is.null(model) && !is.null(model2)) || (!is.null(model) && !is.null(model2) && model2$R2 > model$R2)) {
					model <- model2
				}
			})
			if (verbose) cat('elapsed time =', round(t[3], 2), ', Ended at ', format(Sys.time(), "%m/%d/%Y - %X"), "\n")
		}

		# If no peaks were found, store default info and continue to the next range
		if (is.null(model) || nrow(model$peaks) == 0) {
			asym <- ifelse(opars.loc$oasym > 0, opars.loc$asymmax, 0)
			infos <- rbind(infos, c(pkfit[k,8], pkfit[k,1:2], 0, opars.loc$addPeaks, asym, model$params$obl, opars.loc$qbl, 0, 100, 100))
			next
		}

		# Store fitting results
		if (k > 1 && pkfit[k,1] >= pkfit[k-1,2]) {
			subsetPeaks <- model$peaks[model$peaks$ppm > pkfit[k-1,2], , drop=F]
			Peaks <- rbind(Peaks, subsetPeaks)
			Ymodel <- Ymodel + Rnmr1D::specModel(spec, c(pkfit[k,1], pkfit[k,2]), subsetPeaks)
			iseq <- getseq(spec, c(pkfit[k-1,2], pkfit[k,2]))
			Y[iseq] <- Y[iseq] - model$LB[iseq]
		} else {
			Peaks <- rbind(Peaks, model$peaks)
			Ymodel <- model$model
			Y <- Y - model$LB
		}

		# Store peak fitting metadata
		iseq <- model$iseq
		residus <- Y - Ymodel
		Imodel <- (Ymodel[iseq[1]] + Ymodel[iseq[length(iseq)]]) / 2 + sum(Ymodel[iseq])
		Ispec <- (Y[iseq[1]] + Y[iseq[length(iseq)]]) / 2 + sum(Y[iseq])
		Idiff <- round(100 * (Imodel - Ispec) / Ispec, 4)
		asym <- ifelse(opars.loc$oasym > 0, opars.loc$asymmax, 0)
		infos <- rbind(infos, c(spec$expno, pkfit[k,8], pkfit[k,1:2], model$nbpeak, opars.loc$addPeaks, asym, model$params$obl, opars.loc$qbl, round(model$R2, 4), Idiff))

		if (verbose) cat("-------------------\n")
	}

	# Cleanup to reduce the final object size
	spec$int <- spec$intcorr
	spec$fid <- NULL
	spec$intcorr <- NULL
	spec$BL <- NULL

	# Convert collected info to a structured dataframe
	if (!is.null(Peaks)) {
		rownames(Peaks) <- NULL
		infos <- data.frame(infos[,1:ncol(infos)])
		colnames(infos) <- c('expno', 'zone', 'ppm1', 'ppm2', 'nbpeaks', 'addpeaks', 'asym', 'obl', 'qbl', 'R2', 'Int. Diff. %')
	}

	# Construct the final fitting object
	ppmrange <- c(min(pkfit[,1]), max(pkfit[,2]))
	fit <- list(Y=Y, Ymodel=Ymodel, peaks=Peaks, infos=infos, ppmrange=ppmrange)
	spec$fit <- fit
	spec
})

# Similar to applyPeakFitting1, but introducing parallel processing to handle peak fitting in multiple zones simultaneously.
internalClass$set("private", "applyPeakFitting2", function(spec, opars, zones=NULL, ncpu=2, verbose=1)
{
	# Print default parameters if verbosity is enabled
	if (verbose) {
		cat("-------------------\n")
		cat('Default Parameters: ratioPN =', opars$ratioPN,', asymetric =',opars$oasym,', lowPeaks =',opars$lowPeaks,', addPeaks =',opars$addPeaks,', sndpass =',opars$sndpass, "\n")
		cat("NCPU =",ncpu,"\n")
		cat("-------------------\n")
	}

	# Retrieve peak fitting profile
	pkfit <- PROFILE$fitting

	# Filter peak fitting data based on specified zones
	if (!is.null(zones)) pkfit <- pkfit[ pkfit$zone %in% zones, , drop=F]

	# If no matching zones, return the original spectrum
	if (nrow(pkfit)==0) {
		cat("Error: There is no zone corresponding to those given as input")
		return(spec)
	}

	# Constants for baseline correction and peak fitting
	BLSIG <- 10             # Baseline significance threshold
	ratioPN <- 5            # Ratio threshold for peak detection

	# Function to merge peak fitting results from parallel computations
	combine_list <- function(LL1, LL2) {
		getList <- function(L) { list(id=L$id, spec=L$spec, LB=L$LB, Ymodel=L$Ymodel, peaks=L$peaks,  infos=L$infos) }
		mylist <- list()
		for ( i in 1:length(LL1) ) mylist[[LL1[[i]]$id]] <- getList(LL1[[i]])
		for ( i in 1:length(LL2) ) mylist[[LL2[[i]]$id]] <- getList(LL2[[i]])
		return(mylist)
	}

	# Initialize parallel processing cluster
	NCPU <- min(nrow(pkfit), ncpu)
	cl <- parallel::makeCluster(NCPU)
	doParallel::registerDoParallel(cl)
	Sys.sleep(1)

	# Process all peak fitting zones in parallel
	rq1d <- self
	out <- foreach::foreach(k=1:nrow(pkfit),
		.combine=combine_list,
		.packages=c('Rnmr1D','RnmrQuant1D'),
		.export=ls(globalenv())) %dopar%
	{
		priv <- rq1d$.__enclos_env__$private

	# Initialize variables for peak fitting
		Y <- spec$int <- spec$intcorr
		LB <- Ymodel <- rep(0,length(spec$int))
		Peaks <- infos <- NULL

	# Define the PPM range for peak fitting
		ppmrange <- c(pkfit[k,1],pkfit[k,2])
		if (verbose) cat("PPM range = [", pkfit[k,1],',',pkfit[k,2],"]","\n")

	# Extract and adjust peak fitting parameters
		opars.loc <- opars
		opars.loc$oasym <- ifelse(pkfit$asym[k]>0, 1, 0)
		opars.loc$asymmax <- pkfit$asym[k]
		opars.loc$pvoigt <- ifelse(pkfit$etamin[k]<1, 1, 0)
		opars.loc$etamin <- pkfit$etamin[k]
		opars.loc$addPeaks <- pkfit$addpeaks[k]
		opars.loc$qbl <- pkfit$qbl[k]

	# Extract filter set for peak fitting
		filters <- rq1d$filtersets[[min(length(rq1d$filtersets), round(abs(pkfit$filters[k])))]]

	# Parse obl (baseline correction integrated into the model)
		if (grepl(',',pkfit$obl[k])) {
			obl <- strtoi(simplify2array(strsplit(pkfit$obl[k],',')))
		} else {
			obl <- strtoi(pkfit$obl[k])
		}
		
	# Print fitting parameters for this zone if verbose
		if (verbose) cat('Parameters: addPeaks =',opars.loc$addPeaks,', qbl =',opars.loc$qbl,', oblset =',pkfit$obl[k], "\n")
		if (verbose) cat('Parameters: asymmax =',opars.loc$asymmax,', etamin =',opars.loc$etamin, "\n")
		if (verbose) cat('filter =',paste(filters$main,collapse=','),"\n")

	# Apply baseline correction (qbl) if necessary
		if (opars.loc$qbl) {
			if (opars.loc$qbl==1) { BL <- priv$qnmrbc(spec, ppmrange, BLSIG, 5) }
			else                  { BL <- priv$getBLexternal(spec, opars.loc$qbl, 2, ppmrange) }
			spec$int <- spec$int - BL
			Y <- Y - BL
			LB <- BL
		}

	# Perform peak fitting using Rnmr1D::LSDeconv function
		opars.save <- opars.loc
		t <- system.time({
			model <- Rnmr1D::LSDeconv(spec, ppmrange, opars.loc, filters$main, obl, verbose = verbose)
		})
		if (verbose) cat('elapsed time =', round(t[3],2),', Ended at ',format(Sys.time(), "%m/%d/%Y - %X"),"\n")

	# In case R2 was under R2limit, proceed peak fitting based on the alternative filters
		opars.loc <- opars.save
		if (is.null(model) || (!is.null(model) && model$R2<opars.loc$R2limit)) {
			t<-system.time({
				opars.loc$peaks <- NULL
				others_filters <- filters$others[! filters$others %in% filters$main]
				model2 <- Rnmr1D::LSDeconv(spec, ppmrange, opars.loc, others_filters, obl, verbose = verbose)
				if ( (is.null(model) && !is.null(model2)) || (!is.null(model) && !is.null(model2) && model2$R2>model$R2) )
						model <- model2
			})
			if (verbose) cat('elapsed time =', round(t[3],2),', Ended at ',format(Sys.time(), "%m/%d/%Y - %X"),"\n")
		}

		if (!is.null(model) && nrow(model$peaks)>0) {
			# Accumulating fitting results
			if (k>1 && pkfit[k,1]>=pkfit[k-1,2]) {
				Peaks <- model$peaks[ model$peaks$ppm>pkfit[k-1,2], ]
				Ymodel <- Rnmr1D::specModel(spec, c(pkfit[k,1], pkfit[k,2]), Peaks)
				iseq <- priv$getseq(spec, c(pkfit[k-1,2], pkfit[k,2]))
				LB[iseq] <- LB[iseq] + model$LB[iseq]
				Y[iseq] <- Y[iseq] - model$LB[iseq]
			} else {
				Peaks <- model$peaks
				Ymodel <- model$model
				LB <- LB + model$LB
				Y <- Y - model$LB
			}
			residus <- Y - Ymodel
			# Accumulating fitting infos
			iseq <- model$iseq
			Imodel <- (Ymodel[iseq[1]]+Ymodel[iseq[length(iseq)]])/2 + sum(Ymodel[iseq])
			Ispec <- (Y[iseq[1]]+Y[iseq[length(iseq)]])/2 + sum(Y[iseq])
			Idiff <- round(100*(Imodel-Ispec)/Ispec,4)
			asym <- ifelse(opars.loc$oasym>0, opars.loc$asymmax, 0)
			infos <- c(spec$expno, pkfit[k,8], pkfit[k,1:2], model$nbpeak, opars.loc$addPeaks, asym, model$params$obl, opars.loc$qbl, round(model$R2,4), Idiff, round(t[3],2))
		} else  {
			Peaks <- NULL
			asym <- ifelse(opars.loc$oasym>0, opars.loc$asymmax, 0)
			infos <- c(spec$expno, pkfit[k,8], pkfit[k,1:2], 0, opars.loc$addPeaks, asym, 0, opars.loc$qbl, 0, 100, round(t[3],2))
		}
		if (verbose) cat("-------------------\n")

	# Store fitting results in a list
		mylist <- list()
		mylist[[paste0('S',k)]] <- list( id=paste0('S',k), spec=spec, LB=LB, Ymodel=Ymodel, peaks=Peaks, infos=infos )
		return(mylist)
	}

	# Stop parallel cluster
	invisible(gc())
	parallel::stopCluster(cl)

	# Merge results from parallel processing
	Y <- spec$int <- spec$intcorr
	Ymodel <- rep(0,length(spec$int))
	Peaks <- infos <- NULL
	for(k in 1:length(out)) {
		Y <- Y - out[[k]]$LB
		Ymodel <- Ymodel + out[[k]]$Ymodel
		Peaks <- rbind(Peaks, out[[k]]$peaks)
		infos <- rbind(infos, out[[k]]$infos)
	}

	# Clean up the final spectrum object
	spec$fid <- NULL
	spec$intcorr <- NULL
	spec$BL <- NULL

	infos <- data.frame(infos[,1:ncol(infos)])
	colnames(infos) <- c('expno','zone','ppm1','ppm2','nbpeaks','addpeaks','asym','obl','qbl','R2','Int. Diff. %','Time')

	# Format the peak fitting information table
	if (!is.null(Peaks)) rownames(Peaks) <- NULL

	# Store the final fitted spectrum in the output object
	ppmrange <- c(min(pkfit[,1]), max(pkfit[,2]))
	fit <- list(Y=Y, Ymodel=Ymodel, peaks=Peaks, infos=infos, ppmrange=ppmrange)
	spec$fit <- fit
	spec
})

internalClass$set("private", "applyQuantification", function(spec, fullPattern=TRUE, verbose=1)
{
	# Extract the quantification zones
	ppmview <- spec$fit$ppmrange
	quantifs <- cbind( PROFILE$quantif, get_quantif_ppmrange(spec, PROFILE) )
	pkcmpd <- quantifs[ quantifs$P1>ppmview[1] & quantifs$P1<ppmview[2], ]

	idxint <- 7 # # Index for accessing peak intensity in the peak list

	# Identify unique compounds in the quantification zones
	cmpds <- unique(sort(pkcmpd$compound))

	# Extract the peak list from the spectrum fit
	peaks <- spec$fit$peaks
	rownames(peaks) <- 1:nrow(peaks)

	# Compute the noise level (Vnoise) for peak intensity thresholding
	Vnoise <- get_Vnoise(spec, PPM_NOISE)
	if (verbose) cat("\nVnoise =",Vnoise,"\nB      =",spec$B,"\n\n")

	# Initialize storage for quantification results
	quantification <- NULL
	peaklist <- NULL

	# Loop through each compound to quantify its concentration
	for (cmpd in cmpds) {
		ZQ <- pkcmpd[pkcmpd$compound==cmpd, , drop=FALSE]
		if (verbose) cat("Compound :",cmpd,"\n")
		ISUM <- 0
		ppmrange <- c(999,-999)
		PKZQ <- NULL

		# Loop through each peak zone related to the compound
		for (i in 1:nrow(ZQ)) {
			Pattern <- ZQ$pattern[i]
			ppm1 <- ZQ$ppm1[i]; ppm2 <- ZQ$ppm2[i]
			PK <- NULL
		# Extract peak search parameters, which may be single values or comma-separated lists
			if (grepl(',',ZQ$P2[i])) {
				P2 <- as.numeric(simplify2array(strsplit(ZQ$P2[i],',')))
			} else {
				P2 <- as.numeric(ZQ$P2[i])
			}
			if (grepl(',',ZQ$P3[i])) {
				param <- as.numeric(simplify2array(strsplit(ZQ$P3[i],',')))
			} else {
				param <- as.numeric(ZQ$P3[i])
			}

		# peaks within a range (block)
			if (Pattern == 'b') {
				if (verbose) cat("\t\tPPM pattern = ", Pattern, ", Range = [",ppm1,",",ppm2,"], Np=",ZQ$np[i],", Nb=",param,"\n")
				PK <- find_compounds(spec, peaks, list('C1'=c(Pattern, ppm1, ppm2, param)))
		# peaks defined by a singulet
			} else if (Pattern == 's') {
				if (verbose) cat("\t\tPPM pattern = ", Pattern, ZQ$P1[i],", ppm tol.=",ZQ$P2[i],", Rank=",param[1],"\n")
				PK <- find_compounds(spec, peaks, list('C1'=c(Pattern, ZQ$P1[i], P2, param)))
		# peaks defined by a doublet, triplet or doublet of doublet
			} else if (Pattern  %in% c('d','t','dd')) {
				if (length(param)>3) {
					if (verbose) cat("\t\tPPM pattern = ", Pattern, ZQ$P1[i],", J=",ZQ$P2[i],", params=",param[1],", ",param[2],",",param[3],",",param[4],"\n")
				} else if (length(param)>2) {
					if (verbose) cat("\t\tPPM pattern = ", Pattern, ZQ$P1[i],", J=",ZQ$P2[i],", params=",param[1],", ",param[2],",",param[3],"\n")
				} else if (length(param)>1) {
					if (verbose) cat("\t\tPPM pattern = ", Pattern, ZQ$P1[i],", J=",ZQ$P2[i],", params=",param[1],", ",param[2],"\n")
				} else {
					if (verbose) cat("\t\tPPM pattern = ", Pattern, ZQ$P1[i],", J=",ZQ$P2[i],", Criterion=",param[1],"\n")
				}
				PK <- find_compounds(spec, peaks, list('C1'=c(Pattern, ZQ$P1[i], P2, param)))
		# peaks defined by a rule rN
			} else if (grepl("^r[0-9]+", Pattern)) {
				if (length(param)>2) {
					if (verbose) cat("\t\tPPM Rule = ", Pattern, ", ppm1=",ZQ$P1[i],", ppm2 or J=",ZQ$P2[i],", params=",param[1],",",param[2],",",param[3],"\n")
				} else if (length(param)>1) {
					if (verbose) cat("\t\tPPM Rule = ", Pattern, ", ppm1=",ZQ$P1[i],", ppm2 or J=",ZQ$P2[i],", params=",param[1],",",param[2],"\n")
				} else {
					if (verbose) cat("\t\tPPM Rule = ", Pattern, ", ppm1=",ZQ$P1[i],", ppm2 or J=",ZQ$P2[i],", param=",param,"\n")
				}
				PK <- find_compounds(spec, peaks, list('C1'=c(Pattern, ZQ$P1[i], P2, param)))
		# peaks defined by another pattern
			} else {
				if (verbose) cat("\t\tPPM pattern = ", Pattern, ", ppm1=",ZQ$P1[i],", ppm2 or J=",ZQ$P2[i],", param=",param,"\n")
				PK <- find_compounds(spec, peaks, list('C1'=c(Pattern, ZQ$P1[i], P2, param)))
			}

		# If peaks are found, compute the integral sum
			if (!is.null(PK)) {
				if (verbose) cat("\t\tPeaks = ", paste(PK, collapse=","), "\n")
				PK <- as.numeric(PK$C1)
				if (!is.na(ISUM)) for (k in 1:length(PK))
					ISUM <- ISUM + ZQ$factor[i]*peaks[PK[k],idxint]
			} else {
				ISUM <- NA # Mark as missing if no peaks found
				if (verbose) cat("\t\tPattern not found\n")
			}

		# Update ppm range covered by the compound
			ppmrange <- c( min(ppmrange[1], ppm1), max(ppmrange[2], ppm2) )

		# Store identified peaks if found
			if (!is.null(PK)) {
				PKZQ <- c(PKZQ, PK)
				P1 <- peaks[ rownames(peaks)[PK], ]
				if (verbose) print(P1)
			} else {
		# If fullPattern mode is on, all peaks must be found, otherwise, mark the compound as not found
				if (fullPattern) break
			}
		}
		if (verbose) cat("\n\n")

	# Compute Signal-to-Noise Ratio (SNR) using the highest peak intensity
		SNR <- NA
		if (!is.null(PKZQ)) SNR <- round(max(peaks[ rownames(peaks)[PKZQ], ]$amp)/(2*Vnoise))

	# Store quantification results and peaklist
		quantification <- rbind(quantification,c(cmpd, ppmrange, ISUM, SNR))
		if (!is.null(PKZQ)) peaklist <- rbind(peaklist, c( cmpd, paste0(PKZQ, collapse=',') ))
		else                peaklist <- rbind(peaklist, c( cmpd, NA ))
	}

	# Convert specific columns to numeric for proper formatting
	numformat <- function(df, cols) { for (i in cols) df[,i] <- as.numeric(df[,i]); df }
	quantification <- numformat( data.frame(quantification), c(2:5) )
	colnames(quantification) <- c('Compound','From','To','Integral', 'SNR')
	rownames(quantification) <- 1:nrow(quantification)
	
	# Return the final quantification results and peak list
	list(quantification=quantification, peaklist=peaklist)
})

