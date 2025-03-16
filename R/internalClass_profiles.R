#=====================================================================
# Read the different sections defining the profile, namely :
# --- All profiles ---
#   'preprocess' : parameters for preprocessing (LB, ZFFAC, PHC0)
#   'exclude'  : zones to be excluded for all processing steps,
#   'zeroneg'  : zones to zero before the baseline correction,
#   'baseline' : zones for baseline correction along with their airPLS parameters (lambda & order)
#   'fitting'  : zones for deconvolution (peak fitting) along with their parameters
#      * obl : polynomial order for a local baseline correction optimized at the same time that the peak fitting (default:0)
#      * asym : allows asymetric peaks (default:0)
#      * zone : defines a zone id for better selection
#   'quantif'  : zones to be quantified and assigned to a compound name
#   'compound' : list of compounds along with their features (Mw, Np, Pattern)
#=====================================================================

internalClass$set("public", "readProfile", function(PROFILE)
{
	if (!file.exists(PROFILE))
		 stop("You must specified an existing profile")

	CMDTEXT <- tryCatch({
		# Read the profile file
		readLines(PROFILE)
	},
	error=function(cond) {
		stop(cond)
	})
	CMD <- CMDTEXT[ grep( "^[^(\t,#) ]", CMDTEXT ) ]

	preprocess <- NULL
	exclude_zones <- NULL
	zeroneg <- NULL
	baseline <- NULL
	fitting <- NULL
	quantif <- NULL
	compound <- NULL

	while ( length(CMD)>0 ) {
		cmdLine <- CMD[1]
		cmdPars <- unlist(strsplit(cmdLine[1],"\t"))
		type <- cmdPars[1]
		repeat {
			if (type == 'preprocess' && length(cmdPars)>3) {
				preprocess <- list(LB=as.numeric(cmdPars[2]), ZFFAC=as.numeric(cmdPars[3]), PHC0=as.numeric(cmdPars[4]))
				if (length(cmdPars)>4) preprocess$MVPZTSP <- as.numeric(cmdPars[5])
				if (length(cmdPars)>4) preprocess$DHZPZRANGE <- as.numeric(cmdPars[6])
				CMD <- CMD[-1]
				break
			}
			if (type == 'exclude' && length(cmdPars)>2) {
				exclude_zones <- rbind( exclude_zones, cmdPars[2:3])
				CMD <- CMD[-1]
				break
			}
			if (type == 'zeroneg' && length(cmdPars)>2) {
				zeroneg <- rbind( zeroneg, cmdPars[2:3])
				CMD <- CMD[-1]
				break
			}
			if (type == 'baseline' && length(cmdPars)>5) {
				baseline <- rbind( baseline, cmdPars[2:6])
				CMD <- CMD[-1]
				break
			}
			if (type == 'fitting' && length(cmdPars)>8) {
				fitting <- rbind( fitting, cmdPars[2:9] )
				CMD <- CMD[-1]
				break
			}
			if (type == 'quantif' && length(cmdPars)>8) {
				quantif <- rbind( quantif, cmdPars[2:9] )
				CMD <- CMD[-1]
				break
			}
			if (type == 'compound' && length(cmdPars)>2) {
				compound <- rbind( compound, cmdPars[2:3] )
				CMD <- CMD[-1]
				break
			}
			break
		}
	}

	numformat <- function(df, cols) { for (i in cols) df[,i] <- as.numeric(df[,i]); df }

	if (!is.null(exclude_zones)) {
		exclude_zones <- numformat( data.frame(exclude_zones, stringsAsFactors = FALSE), 1:2 )
		colnames(exclude_zones) <- c('ppm1','ppm2')
	}

	if (!is.null(zeroneg)) {
		zeroneg <- numformat( data.frame(zeroneg, stringsAsFactors = FALSE), 1:2 )
		colnames(zeroneg) <- c('ppm1','ppm2')
	}

	if (!is.null(baseline)) {
		baseline <- numformat( data.frame(baseline, stringsAsFactors = FALSE), 1:5 )
		colnames(baseline) <- c('ppm1','ppm2','lambda','order', 'smooth')
	}

	if (!is.null(fitting)) {
		fitting <- numformat( data.frame(fitting, stringsAsFactors = FALSE), c(1:2,4:8) )
		colnames(fitting) <- c('ppm1','ppm2','obl','qbl','asym','etamin','filters','zone')
	}

	if (!is.null(quantif)) {
		quantif <- numformat( data.frame(quantif, stringsAsFactors = FALSE), c(3,6:8) )
		colnames(quantif) <- c('compound', 'pattern', 'P1', 'P2', 'P3', 'np', 'factor', 'zone')
	}

	if (!is.null(compound)) {
		compound <- numformat( data.frame(compound, stringsAsFactors = FALSE), 2 )
		colnames(compound) <- c('name','mw')
	}

	list(preprocess=preprocess, exclude_zones=exclude_zones, zeroneg=zeroneg, baseline=baseline,
		 fitting=fitting, quantif=quantif, compound=compound)
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
