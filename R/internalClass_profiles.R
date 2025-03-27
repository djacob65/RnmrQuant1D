#=====================================================================
# Read the different sections defining the profile, namely :
#   'preprocess' : parameters for preprocessing (LB, ZFFAC, PHC0)
#   'exclude'  : zones to be excluded for all processing steps,
#   'zeroneg'  : zones to zero before the baseline correction,
#   'baseline' : zones for baseline correction along with their airPLS parameters (lambda & order)
#   'fitting'  : zones for deconvolution (peak fitting) along with their parameters
#      * obl : polynomial order for a local baseline correction optimized at the same time that the peak fitting (default:0)
#      * qbl : indicates if a q-NMR baseline correction will be applied before deconvolution
#      * asym : allows asymetric peaks (default:0)
#      * filters : indicates what type of filters will be applied on spectra
#      * addpeaks : indicates if peaks will be added after a first deconvolution to fill the "holes"
#      * zone : defines a zone id for better selection
#   'quantif'  : zones to be quantified and assigned to a compound name
#      * compound : compound name
#      * pattern : pattern for the peak seach 
#      * P1, P2, P3 : parameters relatives to pattern
#      * np : number of proton
#      * factor : correction factor to be applied to the final concentration
#      * zone : defines a zone id corresponding to a fitting zone
#   'compound' : list of compounds along with their Molecular Weight (MW)
#=====================================================================

internalClass$set("public", "readProfile", function(PROFILE)
{
	# Check if the profile file exists
	if (!file.exists(PROFILE))
		 stop("You must specify an existing profile")

	# Try to read the profile file, stop execution if an error occurs
	CMDTEXT <- tryCatch({
		readLines(PROFILE)  # Read file contents
	},
	error=function(cond) {
		stop(cond)
	})

	# Filter out lines that start with tabs, commas, hash (#), or spaces
	CMD <- CMDTEXT[ grep( "^[^(\t,#,\n) ]", CMDTEXT ) ]

	# Initialize variables to store different types of profile data
	preprocess <- NULL
	exclude_zones <- NULL
	zeroneg <- NULL
	baseline <- NULL
	fitting <- NULL
	quantif <- NULL
	compound <- NULL

	# Process each command line in CMD
	while ( length(CMD) > 0 ) {
		cmdLine <- CMD[1]
		cmdPars <- unlist(strsplit(cmdLine[1],"\t"))  # Split line into parameters
		type <- cmdPars[1]  # Get command type

		# Process different command types
		repeat {
			if (type == 'preprocess' && length(cmdPars) > 3) {
				# Extract preprocessing parameters
				preprocess <- list(LB=as.numeric(cmdPars[2]), ZFFAC=as.numeric(cmdPars[3]), PHC0=as.numeric(cmdPars[4]))
				if (length(cmdPars) > 4) preprocess$MVPZTSP <- as.numeric(cmdPars[5])
				if (length(cmdPars) > 5) preprocess$DHZPZRANGE <- as.numeric(cmdPars[6])
				CMD <- CMD[-1]  # Remove processed line
				break
			}
			# deprecated : this section will not be supported in the futur
			if (type == 'exclude' && length(cmdPars) > 2) {
				# Store exclusion zones
				exclude_zones <- rbind(exclude_zones, cmdPars[2:3])
				CMD <- CMD[-1]
				break
			}
			# deprecated : this section will not be supported in the futur
			if (type == 'zeroneg' && length(cmdPars) > 2) {
				# Store zero-negative zones
				zeroneg <- rbind(zeroneg, cmdPars[2:3])
				CMD <- CMD[-1]
				break
			}
			# deprecated : this section will not be supported in the futur
			if (type == 'baseline' && length(cmdPars) > 5) {
				# Store baseline correction parameters
				baseline <- rbind(baseline, cmdPars[2:6])
				CMD <- CMD[-1]
				break
			}
			if (type == 'fitting') {
				if (length(cmdPars) < 10)
					stop_quietly("Error: the 'fitting' section must have 10 columns in the quantification profile")
				# Store peak fitting parameters
				fitting <- rbind(fitting, cmdPars[2:10])
				CMD <- CMD[-1]
				break
			}
			if (type == 'quantif') {
				if (length(cmdPars) < 9)
					stop_quietly("Error: the 'quantif' section must have 9 columns in the quantification profile")
				# Store quantification parameters
				quantif <- rbind(quantif, cmdPars[2:9])
				CMD <- CMD[-1]
				break
			}
			if (type == 'compound') {
				if (length(cmdPars) < 3)
					stop_quietly("Error: the 'compound' section must have 3 columns in the quantification profile")
				# Store compound information
				compound <- rbind(compound, cmdPars[2:3])
				CMD <- CMD[-1]
				break
			}
			CMD <- CMD[-1]
			break
		}
	}

	# Helper function to convert specified columns to numeric
	numformat <- function(df, cols) {
		for (i in cols) df[,i] <- as.numeric(df[,i]) 
		df 
	}

	# Convert and assign column names to processed data
	if (!is.null(exclude_zones)) {
		exclude_zones <- numformat(data.frame(exclude_zones, stringsAsFactors = FALSE), 1:2)
		colnames(exclude_zones) <- c('ppm1', 'ppm2')
	}

	if (!is.null(zeroneg)) {
		zeroneg <- numformat(data.frame(zeroneg, stringsAsFactors = FALSE), 1:2)
		colnames(zeroneg) <- c('ppm1', 'ppm2')
	}

	if (!is.null(baseline)) {
		baseline <- numformat(data.frame(baseline, stringsAsFactors = FALSE), 1:5)
		colnames(baseline) <- c('ppm1', 'ppm2', 'lambda', 'order', 'smooth')
	}

	if (!is.null(fitting)) {
		fitting <- numformat(data.frame(fitting, stringsAsFactors = FALSE), c(1:2,4:9))
		colnames(fitting) <- c('ppm1', 'ppm2', 'obl', 'qbl', 'asym', 'etamin', 'filters', 'addpeaks', 'zone')
	}

	if (!is.null(quantif)) {
		quantif <- numformat(data.frame(quantif, stringsAsFactors = FALSE), c(3,6:8))
		colnames(quantif) <- c('compound', 'pattern', 'P1', 'P2', 'P3', 'np', 'factor', 'zone')
	}

	if (!is.null(compound)) {
		compound <- numformat(data.frame(compound, stringsAsFactors = FALSE), 2)
		colnames(compound) <- c('name', 'mw')
	}

	# Return processed profile data as a list
	list(preprocess=preprocess, exclude_zones=exclude_zones, zeroneg=zeroneg, baseline=baseline,
		 fitting=fitting, quantif=quantif, compound=compound)
})

# Orders the zones according to the ppm scale then renumbers the zones in both fitting and quantification section
internalClass$set("public", "reorderProfile", function()
{
	fitting <- PROFILE$fitting
	quantif <- PROFILE$quantif
	fitting <- fitting[order(fitting$ppm1), ]
	quantif$zone <- simplify2array(lapply(quantif$zone, function(k){which(k==fitting$zone)}))
	fitting$zone <- 1:nrow(fitting)
	quantif <- quantif[order(quantif$zone), ]
	PROFILE$fitting <<- fitting
	PROFILE$quantif <<- quantif
})

# Save the current quantification profile in an external file
internalClass$set("public", "saveProfile", function(PROFILENAME)
{
	V <- cbind('preprocess', t(PROFILE$preprocess))
	colnames(V)[1] <- '#TYPE'
	write.table(V, PROFILENAME, append = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE, quote=FALSE)
	
	V <- "\n\n"
	write.table(V, PROFILENAME, append = TRUE, sep = "\t", dec = ".", row.names = FALSE, col.names = FALSE, quote=FALSE)
	V <- cbind(rep('fitting',nrow(PROFILE$fitting)), PROFILE$fitting)
	colnames(V)[1] <- '#TYPE'
	write.table(V, PROFILENAME, append = TRUE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE, quote=FALSE)

	V <- "\n\n"
	write.table(V, PROFILENAME, append = TRUE, sep = "\t", dec = ".", row.names = FALSE, col.names = FALSE, quote=FALSE)

	V <- cbind(rep('quantif',nrow(PROFILE$quantif)), PROFILE$quantif)
	colnames(V)[1] <- '#TYPE'
	write.table(V, PROFILENAME, append = TRUE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE, quote=FALSE)

	V <- "\n\n"
	write.table(V, PROFILENAME, append = TRUE, sep = "\t", dec = ".", row.names = FALSE, col.names = FALSE, quote=FALSE)

	V <- cbind(rep('compound',nrow(PROFILE$compound)), PROFILE$compound)
	colnames(V)[1] <- '#TYPE'
	write.table(V, PROFILENAME, append = TRUE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE, quote=FALSE)

})

# Get quantif parameters as a list for a given zone
internalClass$set("private", "get_quantif_cmpds", function(profil, zone)
{
	quantifs <- PROFILE$quantif[ PROFILE$quantif$zone == zone, ]
	cmpds <- list()
	if (nrow(quantifs)==0) return(cmpds)
	n <- 1
	for (x in 1:nrow(quantifs)) {
		if (grepl(',',quantifs[x, 4])) {
			p2 <- as.numeric(simplify2array(strsplit(quantifs[x, 4],',')))
		} else {
			p2 <- as.numeric(quantifs[x, 4])
		}
		if (grepl(',',quantifs[x, 5])) {
			p3 <- as.numeric(simplify2array(strsplit(quantifs[x, 5],',')))
		} else {
			p3 <- as.numeric(quantifs[x, 5])
		}
		cmpd <- quantifs[x, 1]
		if (is.null(cmpds[[cmpd]])) { n <- 1 } else { n <- n + 1; cmpd <- paste0(cmpd,n) }
		cmpds[[cmpd]] <- c(quantifs[x, 2], as.numeric(quantifs[x, 3]), p2, p3 )
	}
	cmpds
})

