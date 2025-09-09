#=====================================================================
# Find patterns / compounds
#=====================================================================

# Get the ppm range for each quantification zone
internalClass$set("private", "get_quantif_ppmrange", function(spec=NULL, profil=NULL)
{
	if (is.null(spec) && is.na(as.numeric(FIELD)))
		stop_quietly("Error in 'get_quantif_ppmrange' : no 'spec' specified and 'FIELD' is not a numeric!")

	SFO1 <- ifelse( is.null(spec), as.numeric(FIELD), spec$acq$SFO1 )
	if (is.null(profil)) profil <- PROFILE

	# Calculate the ppm for each quantification zone
	V <- NULL
	pkcmpd <- profil$quantif
	for (k in 1:nrow(pkcmpd)) {
		P1 <- pkcmpd[k,]$P1
		pattern <- pkcmpd[k,]$pattern
		if (grepl(',',pkcmpd[k,]$P2)) {
			P2 <- as.numeric(simplify2array(strsplit(pkcmpd[k,]$P2,',')))
		} else {
			P2 <- as.numeric(pkcmpd[k,]$P2)
		}
		dppm <- NULL
		if (pattern == 'r1')	V <- rbind( V, c(P1, P2) )
		if (pattern == 'r2')	V <- rbind( V, c(P1, P2) )
		if (pattern == 'r3')	dppm <- 0.75*P2/SFO1
		if (pattern == 'r4')	dppm <- 0.75*P2/SFO1
		if (pattern == 'r5')	dppm <- 0.75*P2/SFO1
		if (pattern == 'r6')	dppm <- 0.75*P2/SFO1
		if (pattern == 'r7')	V <- rbind( V, c(P1, P2) )
		if (pattern == 'r8')	V <- rbind( V, c(P1, P2) )
		if (pattern == 'r9')	V <- rbind( V, c(P1, P2) )
		if (pattern == 'b')	 	V <- rbind( V, c(P1, P2) )
		if (pattern == 's')		dppm <- 1.5*P2
		if (pattern == 'd')		dppm <- 0.75*P2/SFO1
		if (pattern == 't')		dppm <- 1.35*P2/SFO1
		if (pattern == 'dd')	dppm <- (0.5*P2[1]+P2[2])/SFO1
		if (pattern == 'q')		dppm <- 2*P2/SFO1
		if (pattern == 'm')		dppm <- 2.5*P2/SFO1
		if (pattern == 'm2')	dppm <- 2.5*P2/SFO1
		if (is.null(dppm)) 	 next
		V <- rbind( V, c(P1-dppm, P1+dppm) )
	}
	V <- as.data.frame(V)
	colnames(V) <- c('ppm1', 'ppm2')
	V
})

# Get the ppm shift so that the maximum intensity in a ppm zone (pzone1, pzone2) corresponds to the p0 value
internalClass$set("private", "get_ppm_shift", function(spec, pzone1, pzone2, p0, type='s')
{
	iseq <- getseq(spec,c(pzone1,pzone2))
	Y <- spec$int[iseq]
	if (type=='d') {
		DMIN <- round(10*spec$size/65535)
		V <- order(Y, decreasing=T)
		k <- 2; while(abs(V[k]-V[k-1])<DMIN) k <- k+1
		i0 <- round(0.5*(V[1]+V[k]))
	} else {
		i0 <- which(Y == max(Y))
	}
	spec$ppm[iseq[1] + i0 - 1] - p0
})

# singulet (s) : central ppm, tolerance for ppm
internalClass$set("private", "find_pattern_s", function(spec, peaks, ppm0, dppm=0.005, rank=1)
{
	ret <- NULL
	P1 <- peaks[peaks$ppm>ppm0-dppm, ,drop=F]
	P2 <- P1[P1$ppm<ppm0+dppm, ,drop=F]
	if (!is.null(P2) && nrow(P2)>0) {
		P2 <- P2[order(P2$amp, decreasing=T), ]
		rk <- rank
		repeat {
			if (rk==1 || abs(P2$ppm[1]-P2$ppm[rk])>2*P2$sigma[1]) {
				ret <- rownames(P2)[min(rk,nrow(P2))]
				break
			}
			rk <- rk + 1
			if (rk>nrow(P2)) break
		}
	}
	ret
})

# doublet (d) : central ppm, J, sel, (ratio, dJ, dS)
internalClass$set("private", "find_pattern_d", function(spec, peaks, ppm0, J, sel=0, ratio=1.4, dJ=0.3, dS=65, full=FALSE)
{
	# Tolerance settings
	dA <- 0.25*ratio
	dppm0 <- 5/spec$acq$SFO1
	dppm <- dppm0 +0.5*J/spec$acq$SFO1

	# Find the doublet
	P2 <- peaks[peaks$ppm>(ppm0-dppm) & peaks$ppm<(ppm0+dppm), , drop=F]
	JL <- NULL
	repeat {
		if (nrow(P2)<2) break
		for (i in 1:(nrow(P2)-1))
			for (j in (i+1):nrow(P2)) {
				Jij <- (P2$ppm[j]-P2$ppm[i])*spec$acq$SFO1
				if (!is.na(Jij) && Jij>(J-dJ) && Jij<(J+dJ)) {
					Xm <- 0.5*(P2$ppm[i]+P2$ppm[j])
					Am <- 0.5*(P2$amp[i]+P2$amp[j])
					if (P2$amp[i]>P2$amp[j]) { Aratio <- P2$amp[i]/P2$amp[j] } else { Aratio <- P2$amp[j]/P2$amp[i] }
					dSigma <- 200*abs(P2$sigma[i]-P2$sigma[j])/(P2$sigma[i]+P2$sigma[j])
					# full criterion
					if (sel==0) {
						crit <- abs((Jij-J)/(Jij+J)) + abs((Xm-ppm0)/(Xm+ppm0))
					} else if (sel==1) {
						crit <- abs((Jij-J)/(Jij+J)) + 2*abs((Xm-ppm0)/(Xm+ppm0)) + 0.5*(Aratio-1)*(Aratio-1)
					} else {
						crit <- abs((Jij-J)/(Jij+J)) + 0.5*(Aratio-1)*(Aratio-1)
					}
					JL <- rbind(JL, c(rownames(P2)[i],rownames(P2)[j], round(Jij,2), round(Xm,4), round(Am,4), round(Aratio,4), round(dSigma,4), round(crit,4)) )
				}
			}
		if (is.null(JL)) break
		if (nrow(JL)>0) JL <- JL[abs(as.numeric(JL[,6]))<(ratio+dA), , drop=F] # criterion based on the ratio
		if (nrow(JL)>0) JL <- JL[abs(as.numeric(JL[,7]))<dS, , drop=F] # criterion based on dSigma
		if (nrow(JL)>0) JL <- JL[abs(as.numeric(JL[,4])-ppm0)<dppm0, , drop=F] # criterion based on ppm0
		if (nrow(JL)>1 && sel<2) JL <- JL[order(as.numeric(JL[,8])), ] # order by increasing criterion
		if (nrow(JL)>1 && sel>1) JL <- JL[order(as.numeric(JL[,5]), decreasing=T), ] # order by decreasing intensities
		if (nrow(JL)<1) { JL <- NULL; break }
		if (!full) {
			if (nrow(JL)>2) JL <- JL[1:2,, drop=F] # take the 2 doublets having the smallest criteria/intensities values
			# if intensity for the first doublet is not higher than 10% of the second one then ...
			if (nrow(JL)>1 && (min(as.numeric(JL[1:2,5]))/max(as.numeric(JL[1:2,5])))<0.1)
				JL <- JL[which(max(as.numeric(JL[,5]))==as.numeric(JL[,5])), , drop=F] # take the highest amplitude
			JL <- JL[1, 1:2] # take the first doublet
		}
		break
	}
	JL
})

# triplet (t) : central ppm, J
internalClass$set("private", "find_pattern_t", function(spec, peaks, ppm0, J, sel=0, Ramp2=1.3, dJ=0.35)
{
	# Tolerance settings
	Ramp <- 3.5
	dppm <- J/(2*spec$acq$SFO1)

	groups <- NULL
	g <- unique(rbind(find_pattern_d(spec, peaks, ppm0-dppm, J, sel, Ramp, dJ, full=T),
					  find_pattern_d(spec, peaks, ppm0+dppm, J, sel, Ramp, dJ, full=T)))
	repeat {
		if (! c('matrix') %in% class(g)) break
		if (is.null(g) || is.null(nrow(g)) || nrow(g)<2) break
		g <- g[order(as.numeric(g[,8])),]
		rownames(g) <- 1:nrow(g)
		g <- g[rownames(unique(g[,1:7])), , drop=F]
		if (! c('matrix') %in% class(g)) break
	# Doublets must be linked
		h <- NULL
		for (i in 1:nrow(g)) if (g[i,1] %in% g[,2] || g[i,2] %in% g[,1]) h <- rbind(h, g[i,])
		if (is.null(h) || nrow(h)<2) break
	# Apply separation depending on links
		L <- NULL; k <- 0
		h <- h[order(as.numeric(h[,1])),]
		for (i in 1:nrow(h)) {
			if (k==0) { k <- 1; L[[k]] = c(h[i,1], h[i,2]) }
			else {
				ok <- 0
				for (j in 1:length(L)) if (h[i,1] %in% L[[j]]) { L[[j]] <- unique(c(L[[j]], h[i,1], h[i,2]) ); ok <- 1 }
				if (ok==0) { k <- k+1; L[[k]] = c(h[i,1], h[i,2]) }
			}
		}
		if (is.null(L)) break
		L2 <- NULL; k <- 0
		for (j in 1:length(L)) {
			if (length(L[[j]])==3) { k<-k+1; L2[[k]] <- L[[j]] }
			if (length(L[[j]])==4) { k<-k+1; L2[[k]] <- L[[j]][-length(L[[j]])]; k<-k+1; L2[[k]] <- L[[j]][-1] }
		}
		if (is.null(L2)) break
	# Criterion on amplitude
		L <- NULL; k <- 0
		for (j in 1:length(L2)) {
			L3 <- sort(as.numeric(unique(L2[[j]])))
			v <- peaks[L3,]$amp
			if ( (2*v[2]/(v[1]+v[3]))>Ramp2 ) { k<-k+1; L[[k]] <- L3 }
		}
		if (is.null(L)) break
	# Criterion on J (J +/- dJ)
		L2 <- NULL; k <- 0
		for (i in 1:length(L)) {
			v <- L[[i]]
			J1 <- round(abs(peaks$ppm[v[2]]-peaks$ppm[v[1]])*spec$acq$SFO1,2)
			J2 <- round(abs(peaks$ppm[v[3]]-peaks$ppm[v[2]])*spec$acq$SFO1,2)
			if (J1>(J-dJ) && J1<(J+dJ) && J2>(J-dJ) && J2<(J+dJ)) { k<-k+1; L2[[k]] <- v }
		}
		if (is.null(L2)) break
		v <- NULL
		for (i in 1:length(L2)) v <-c(v, peaks$amp[ L2[[i]][2] ])
		groups <- L2[[which(v==max(v))]]
		break
	}
	groups
})

# quadruplet (q) : central ppm, J
internalClass$set("private", "find_pattern_q", function(spec, peaks, ppm0, J, Ramp=3)
{
	# Tolerance settings
	dJ <- 2.5*spec$dppm*spec$acq$SFO1
	facR <- 0.25
	dppm <- J/spec$acq$SFO1

	groups <- NULL
	sel <- 0
	g <- rbind( find_pattern_d(spec, peaks, ppm0, J, sel, Ramp, full=T),
				find_pattern_d(spec, peaks, ppm0-dppm, J, sel, Ramp, full=T),
				find_pattern_d(spec, peaks, ppm0+dppm, J, sel, Ramp, full=T) )
	g <- unique(g)
	repeat {
		if (! c('matrix') %in% class(g)) break
		if (is.null(g) || is.null(nrow(g)) || nrow(g)<2) break
	# Doublets must be linked
		h <- NULL
		for (i in 1:nrow(g)) if (g[i,1] %in% g[,2] || g[i,2] %in% g[,1]) h <- rbind(h, g[i,])
		if (is.null(h) || nrow(h)<2) break
	# Apply separation depending on links
		L <- NULL; k <- 0
		h <- h[order(as.numeric(h[,1])),]
		for (m in 1:(nrow(h)-1)) {
			for (i in m:nrow(h)) {
				if (i==m) { k <- k + 1; L[[k]] = c(h[i,1], h[i,2]) }
				else if (h[i,1] == L[[k]][length(L[[k]])]) { L[[k]] = c(L[[k]], h[i,2]) }
			}
		}
		G <- NULL; k <- 0
		for (j in 1:length(L))
			if (length(L[[j]])==4) {
				V <- as.numeric(L[[j]])
				if (mean(peaks$amp[V[order(V)][2:3]])>mean(peaks$amp[V[order(V)][c(1,4)]])) { k<-k+1; G[[k]] <- V }
			}
		if (is.null(G) && length(L)>2) {
			for (j in 1:length(L)) {
				if (length(L[[j]]) != 3) { G <- c(G,99999); next }
				V <- as.numeric(L[[j]]); V <- V[order(V)]
				G <- c(G, 0.5*(abs(abs(peaks$ppm[V[1]]-peaks$ppm[V[2]])*spec$acq$SFO1-J) + abs(abs(peaks$ppm[V[2]]-peaks$ppm[V[3]])*spec$acq$SFO1-J)))
			}
			if (!is.null(G)) L <- list(L[[which(G==min(G))]]) else L <- NULL
		} else if (!is.null(G)) {
			L <- G
		}
		# In case we missed a peak
		if(length(L)==1 && length(L[[1]])==3) {
			P1 <- peaks[L[[1]], ]
			P1 <- P1[order(P1$ppm), ]
			if (P1$amp[2]>1.5*P1$amp[3]) {
				P2 <- peaks[peaks$ppm<(P1$ppm[1]-0.8*J/spec$acq$SFO1) & peaks$ppm>(P1$ppm[1]-1.2*J/spec$acq$SFO1), , drop=F]
				if (nrow(P2)>1) P2 <- P2[which( abs((P1$ppm[1]-P2$ppm)*spec$acq$SFO1-J)==min(abs((P1$ppm[1]-P2$ppm)*spec$acq$SFO1-J))), , drop=F]
			} else {
				P2 <- peaks[peaks$ppm>(P1$ppm[3]+0.8*J/spec$acq$SFO1) & peaks$ppm<(P1$ppm[3]+1.2*J/spec$acq$SFO1), , drop=F]
				if (nrow(P2)>1) P2 <- P2[which( abs((P1$ppm[3]-P2$ppm)*spec$acq$SFO1-J)==min(abs((P1$ppm[3]-P2$ppm)*spec$acq$SFO1-J))), , drop=F]
			}
			if (nrow(P2)==1) L[[1]] <- c(L[[1]], rownames(P2))
		}
		G <- NULL; k <- 0
		for (j in 1:length(L)) if (length(L[[j]])==4) { k<-k+1; G[[k]] <- L[[j]] }
		if (is.null(G)) break
	## In case there are 2 or more peak lists, take the one with the highest sum of amplitudes
		v <- sapply(1:length(G), function(k) sum(peaks[as.numeric(G[[k]]), ]$amp))
		L <- sort(as.numeric(G[[which(v==max(v))]]))
		if (length(L)<4) break
		if (length(L)>4) {
	# Criterion on J (J +/- dJ)
			M <- matrix(rep(0,length(L)*length(L)), nrow=length(L), byrow=T)
			colnames(M) <- rownames(M) <- L
			for (i in 1:length(L)) for (j in 1:length(L)) {
				if (i==j) next
				Jij <- round(abs(peaks$ppm[L[i]]-peaks$ppm[L[j]])*spec$acq$SFO1,2)
				if (Jij>(J-dJ) && Jij<(J+dJ) )  M[i,j] <- Jij
			}
			idx <- which(apply(M,1,sum)>0)
			if (length(idx)<4) break
			L <- L[ idx ]
			M <- M[idx,idx]
			idx <- which(apply(M,2,sum)>0)
			if (length(idx)<4) break
			groups <- L[ idx ]
		} else {
			groups <- L
		}
		# Criterion on amplitude (1, 3, 3, 1)
		if (length(L)==4 && mean(peaks$amp[L[order(L)][2:3]])<mean(peaks$amp[L[order(L)][c(1,4)]]))
			groups <- NULL

		break
	}
	groups
})

# multiplet (m1) : central ppm, J (quintet/pentet)
internalClass$set("private", "find_pattern_m", function(spec, peaks, ppm0, J)
{
	groups <- find_pattern_t(spec, peaks, ppm0, J, 0, 1.1)
	if (!is.null(groups)) {
		dp0 <- mean(peaks[groups, ]$sigma)
		p1 <- peaks[as.numeric(groups[1]), ]$ppm - J/spec$acq$SFO1
		P1 <- peaks[peaks$ppm>(p1-dp0) & peaks$ppm<(p1+dp0),,drop=F]
		if (nrow(P1)>0) groups <- c(rownames(P1), groups)
		p2 <- peaks[as.numeric(groups[length(groups)]), ]$ppm + J/spec$acq$SFO1
		P2 <- peaks[peaks$ppm>(p2-dp0) & peaks$ppm<(p2+dp0),,drop=F]
		if (nrow(P2)>0) groups <- c(groups, rownames(P2))
	}
	groups
})

# multiplet (m2) : central ppm, J (septet)
internalClass$set("private", "find_pattern_m2", function(spec, peaks, ppm0, J)
{
	groups <- find_pattern_m(spec, peaks, ppm0, J)
	if (!is.null(groups)) {
		dp0 <- median(peaks[groups, ]$sigma)
		p1 <- peaks[as.numeric(groups[1]), ]$ppm - J/spec$acq$SFO1
		P1 <- peaks[peaks$ppm>(p1-dp0) & peaks$ppm<(p1+dp0),,drop=F]
		if (nrow(P1)>0) groups <- c(rownames(P1), groups)
		p2 <- peaks[as.numeric(groups[length(groups)]), ]$ppm + J/spec$acq$SFO1
		P2 <- peaks[peaks$ppm>(p2-dp0) & peaks$ppm<(p2+dp0),,drop=F]
		if (nrow(P2)>0) groups <- c(groups, rownames(P2))
	}
	groups
})

# doublet of doublet (t) : central ppm, J1 (between the doublets) , J2 (each doublet)
internalClass$set("private", "find_pattern_dd", function(spec, peaks, ppm0, J1, J2, crit=0, ratio=2.2)
{
	# Tolerance settings
	dJ2 <- 1.5*spec$dppm*spec$acq$SFO1
	dJ1 <- 2*dJ2
	dppm <- J1/(2*spec$acq$SFO1)
	Ramp <- ratio
	sel <- 0
	ret <- NULL
	g <- rbind( find_pattern_d(spec, peaks, ppm0-dppm, J2, sel, Ramp, full=T),
				find_pattern_d(spec, peaks, ppm0+dppm, J2, sel, Ramp, full=T),
	            find_pattern_d(spec, peaks, ppm0-dppm, J2-dJ2, sel, Ramp, full=T),
				find_pattern_d(spec, peaks, ppm0+dppm, J2-dJ2, sel, Ramp, full=T),
	            find_pattern_d(spec, peaks, ppm0-dppm, J2+dJ2, sel, Ramp, full=T),
				find_pattern_d(spec, peaks, ppm0+dppm, J2+dJ2, sel, Ramp, full=T) )
	g <- unique(g)
	repeat {
		if (! c('matrix') %in% class(g)) break
		if (is.null(g) || is.null(nrow(g)) || nrow(g)<2 ) break
		id <- ifelse(crit==0, 8, 1)
		g <- g[order(as.numeric(g[,id])), ,drop=F]
		rownames(g) <- 1:nrow(g)
		g <- g[rownames(unique(g[,1:7])), ,drop=F]
		if (! c('matrix') %in% class(g)) break
		if (!is.null(g) && nrow(g)>1) {
			n <- nrow(g)
			vmin <- J1
			dpmin <- 999
			im <- jm <- 0
			for (i in 1:(n-1))
				for (j in (i+1):n) {
					dppmij <- abs(ppm0 - 0.5*(as.numeric(g[i,4])+as.numeric(g[j,4])))
					Jij <- abs(as.numeric(g[i,4])-as.numeric(g[j,4]))*spec$acq$SFO1
					Am <- mean(as.numeric(g[i,5]), as.numeric(g[j,5]))
					v <- abs(J1 - Jij)*dppmij/Am
					if (Jij>(J1-dJ1) && Jij<(J1+dJ1) & v<vmin) { im <- i; jm <- j; vmin <- v; dpmin <- dppmij }
				}
			if (im==0 && jm==0) break
			ret <- as.character(sort(as.numeric(c( g[im, 1:2], g[jm, 1:2] ))))
		}
		break
	}
	ret
})

# find peaks in the ppm range
internalClass$set("private", "find_peaks_range", function(spec, peaks, ppm1, ppm2, nbpeaks=0, ratioPN=5)
{
	groups <- NULL
	repeat {
		P1 <- peaks[peaks$ppm>ppm1 & peaks$ppm<ppm2, , drop=F]
		if (is.null(P1) || nrow(P1)<1) break
		P2 <- Rnmr1D::peakFiltering(spec, P1, ratioPN)
		if (is.null(P2) || nrow(P2)<1) break
		rownames(P2) <- rownames(P1)[which( P1$pos %in% P2$pos)]
		if (nbpeaks>0)
			P2 <- P2[which(P2$amp %in% sort(P2$amp, decreasing=T)[1:nbpeaks]), ]
		if (is.null(P2) || (nbpeaks>0 && nrow(P2)<nbpeaks)) break
		if (is.null(P2) || nrow(P2)<1) break
		groups <- rownames(P2)
		break
	}
	groups
})

# Rule r1 : Among the peaks having a S/N greater than 25 (default), take the nth peak (rank) in increasing ppm order, otherwise in decreasing ppm order
internalClass$set("private", "find_peaks_rule_r1", function(spec, peaks, ppm1, ppm2, rank)
{
	groups <- NULL
	P1 <- peaks[peaks$ppm>ppm1 & peaks$ppm<ppm2, , drop=F]
	ratioPN <- 25
	repeat {
		if (rank==0 || is.null(P1) || nrow(P1)==0) break
		pk <- Rnmr1D::peakFiltering(spec,P1, ratioPN)        # take peaks only with S/N above ratioPN
		if (is.null(pk) || nrow(pk)<abs(rank)) break         # break if not enough peaks
		rownames(pk) <- which( peaks$pos %in% pk$pos)        # Get the right indexes
		if (rank<0) pk <- pk[order(pk$ppm, decreasing=T), ]  # order peaks by decreasing ppm if rank<0
		groups <- rownames(pk[abs(rank),])                   # Retrieve the right index  of the peak with the right rank
		break
	}
	groups
})

# Rule r2:  Among the peaks having a S/N greater than 65 (default), take the first peak which is distant at least Jmin Hz but not more than Jmax Hz from the first peak of the selected ppm interval. If Jmin is positive, consider the first peak from the smallest ppm, otherwise from the largest ppm. Jmax will have automatically the same sign that Jmin.
internalClass$set("private", "find_peaks_rule_r2", function(spec, peaks, ppm1, ppm2, Jmin, Jmax=4)
{
	groups <- NULL
	P1 <- peaks[peaks$ppm>ppm1 & peaks$ppm<ppm2, , drop=F]
	ratioPN <- 65
	repeat {
		if (is.null(P1) || nrow(P1)<2) break
		pk <- Rnmr1D::peakFiltering(spec, P1, ratioPN)        # take peaks only with S/N above ratioPN
		if (is.null(pk) || nrow(pk)<2) break                  # break if not enough peaks
		rownames(pk) <- which( peaks$pos %in% pk$pos)         # # Get the right indexes
		if (Jmin<0) pk <- pk[order(pk$ppm, decreasing=T), ]   # order peaks by decreasing ppm if Jmin<0
		n <- ifelse( pk$amp[1]>pk$amp[2], 2, min(3,nrow(pk)))
		for (k in n:nrow(pk)) {                               # foreach peak after the first one
			J <- abs(pk$ppm[k]-pk$ppm[n-1])*spec$acq$SFO1       # Distant in Hz between this peak and the first peak
			if (J>abs(Jmin) && J<abs(Jmax) ) {                # takes this peak if between Jmin and Jmax Hz from the first peak of the interval
				groups <- rownames(pk[k,])                    # retrieve the peak based on its index
				break
			}
		}
		break
	}
	groups
})

# Rule r3 : Among doublets that have a ratio between intensities lower or equal to ratio, take one which has the greatest intensity
internalClass$set("private", "find_peaks_rule_r3", function(spec, peaks, ppm0, J, ratio)
{
	groups <- NULL
	g <- unique( find_pattern_d (spec, peaks, ppm0, J, 0, ratio, full=T) )
	if (!is.null(g)) {
		if (nrow(g)>0)
			g <- g[ as.numeric(g[,4])>(ppm0-0.666*J) | as.numeric(g[,4])<(ppm0+0.666*J), , drop=F]
		if (nrow(g)>1)
			g <- g[ order(as.numeric(g[,5]), decreasing=T), ]
		if (nrow(g)>0 && as.numeric(g[1,6])<1.5*ratio)
			groups <- g[1, c(1:2)]
	}
	groups
})

# Rule r4 : Find the doublet based on the ppm0 and J parameters , then also take the two peaks of greatest intensities located between the two peaks of the doublet.
internalClass$set("private", "find_peaks_rule_r4", function(spec, peaks, ppm0, J, dJ, sel=0)
{
	groups <- NULL
	g <- unique( find_pattern_d (spec, peaks, ppm0, J, sel, ratio=2, full=F) )
	if (!is.null(g) && length(g)==2) {
		P1 <- peaks[ which(rownames(peaks) %in% as.numeric(g)), , drop=F]
		P2 <- peaks[peaks$ppm>min(P1$ppm) & peaks$ppm<max(P1$ppm), , drop=F]
		if (nrow(P2)==2)  {
			groups <- sort(c(g, rownames(P2[1:2, ])))
		} else {
			i <- which((P2$ppm - P1$ppm[1])*spec$acq$SFO1>dJ)[1]
			j <- rev(which((P1$ppm[2] - P2$ppm)*spec$acq$SFO1>dJ))[1]
			if (!is.na(i) && !is.na(j))
				groups <- sort(c(g, rownames(P2[i:j, ])))
		}
		dJ <- 0.5 # 1
		i <- which((P2$ppm - P1$ppm[1])*spec$acq$SFO1<dJ)[1]
		if (!is.na(i))
			groups <- c(groups, rownames(P2[i,]))
		i <- rev(which((P1$ppm[2] - P2$ppm)*spec$acq$SFO1<dJ))[1]
		if (!is.na(i))
			groups <- c(groups, rownames(P2[i,]))

		#if (nrow(P2)>=2) {
		#	P2 <- P2[ order(P2$int, decreasing=T), ]
		#	groups <- sort(c(g, rownames(P2[1:2, ])))
		#}
	}
	unique(groups)
})

# Rule r5 : Find the doublet according to the parameters ppm0 and J but with the central value ppm0 located in the interval [ppm1, ppm2]. Take the doublet with the highest intensity
internalClass$set("private", "find_peaks_rule_r5", function(spec, peaks, ppm0, J, ppm1, ppm2, sel=0)
{
	groups <- NULL
	dppm0 <- 3.5/spec$acq$SFO1
	dppm <- 0.5*J/spec$acq$SFO1
	ratioPN <- 5
	repeat {
		P1 <- Rnmr1D::cleanPeaks(spec, peaks, ratioPN, keeprows=TRUE)
		if (nrow(P1)<2) break
		g <- rbind( find_pattern_d (spec, P1, ppm0-dppm0, J, sel, ratio=2, full=T),
					find_pattern_d (spec, P1, ppm0,       J, sel, ratio=2, full=T),
					find_pattern_d (spec, P1, ppm0+dppm0, J, sel, ratio=2, full=T) )
		g <- unique(g)
		if (is.null(g) || nrow(g)<2) break
		if (nrow(g)>0) g <- g[ as.numeric(g[,4])>(ppm1-dppm) & as.numeric(g[,4])<(ppm2+dppm), , drop=F]
		if (nrow(g)>1) g <- g[ order(as.numeric(g[,5]), decreasing=T),  ]
		if (nrow(g)>0) groups <- g[1, c(1:2)]
		break
	}
	groups
})

# Rule r6 : Find the doublet according to the parameters ppm0 and J but with the J value in the interval J +/- dJ.
internalClass$set("private", "find_peaks_rule_r6", function(spec, peaks, ppm0, J, dJ, sel=0)
{
	groups <- NULL
	g <- unique( find_pattern_d (spec, peaks, ppm0, J, sel, ratio=2, full=T) )
	if (!is.null(g) && nrow(g)>0) {
		if (nrow(g)>0) g <- g[ as.numeric(g[,3])>(J-dJ) & as.numeric(g[,3])<(J+dJ), , drop=F]
		if (nrow(g)>1) g <- g[ order(as.numeric(g[,8])),  ]
		if (nrow(g)>0) groups <- g[1, c(1:2)]
	}
	groups
})

# Rule r7 : Find a multiplet having N peaks (nbpeaks) separated by at most D Hz (dist) and at least 1 Hz. Peaks must be contiguous
internalClass$set("private", "find_peaks_rule_r7", function(spec, peaks, ppm1, ppm2, nbpeaks=2, dist=2)
{
	groups <- NULL
	repeat {
		P1 <- peaks[peaks$ppm>ppm1 & peaks$ppm<ppm2, , drop=F]
		if (is.null(P1) || nrow(P1)<nbpeaks) break
		V <- NULL;
		for (k in 2:nrow(P1)) { # Select contiguous peaks satisfying the distance criteria
			D <- (P1$ppm[k]-P1$ppm[k-1])*spec$acq$SFO1
			if (D<dist) V <- unique(c(V, k-1, k));
			if (D<1 && length(V)>=nbpeaks) break
			if (D<1) V <- NULL
		}
		P1 <- P1[unique(V), ]
		V <- round(100*P1$amp/median(P1$amp))
		P1 <- P1[which(V>50 & V<200), ] # Eliminates peaks with extreme intensities
		if (is.null(P1) || nrow(P1)<nbpeaks) break
		V <- NULL;
		for (k in 2:nrow(P1)) # Select contiguous peaks satisfying the distance criteria
			if (((P1$ppm[k]-P1$ppm[k-1])*spec$acq$SFO1)<dist) V <- c(unique(V), k-1, k)
		P1 <- P1[unique(V), ]
		if (is.null(P1) || nrow(P1)<nbpeaks) break
		groups <- rownames(P1[1:nbpeaks, ]) # Peaks must be contiguous
		break
	}
	groups
})


# Rule r8 : Take the highest intensity peak in the selected area. Then if the some peaks are close to the first one in intensity (ratio<2) and in distance (dist Hz) then also take these peaks.
internalClass$set("private", "find_peaks_rule_r8", function(spec, peaks, ppm1, ppm2, ratio=2, dist=1.2)
{
	groups <- NULL
	repeat {
		P1 <-peaks[ peaks$ppm>ppm1 & peaks$ppm<ppm2, , drop=F]
		if (is.null(P1) || nrow(P1)<1) break
		P1 <- P1[order(as.numeric(P1$amp), decreasing=T), , drop=F]
		groups <- rownames(P1[1,])
		if (nrow(P1)==1) break
		p1 <- P1[1,]
		P1 <- P1[-1, ]
		P1 <- P1[ P1$ppm<(p1$ppm+dist/spec$acq$SFO1) & P1$ppm>(p1$ppm-dist/spec$acq$SFO1), , drop=F]
		if (nrow(P1)>0)
			for (k in 1:nrow(P1))
				if ((p1$amp/P1$amp[k])<ratio) groups <- c(groups, rownames(P1[k, ]))
		break
	}
	groups
})

# r9: find 'nbpeaks' consecutive peaks in the ppm range (e.g shikimic acid)
internalClass$set("private", "find_peaks_rule_r9", function(spec, peaks, ppm1, ppm2, J=2, nbpeaks=5, snrthres=5)
{
	Dmax <- 1.02*(nbpeaks-1)*J
	Dmin <- 0.9*(nbpeaks-1)*J
	facJ <- 1.05
	ratioPN <- 1

	distHz <- function(gn, n1, n2) { (peaks$ppm[gn[n2]] - peaks$ppm[gn[n1]])*spec$acq$SFO1 }
	is_snr <- function(n) { (peaks$amp[n]/spec$Noise)<snrthres }

	groups <- NULL
	repeat {
		# Get peaks within the ppm range
		P1 <- peaks[peaks$ppm>ppm1 & peaks$ppm<ppm2, , drop=F]
		if (is.null(P1) || nrow(P1)<(nbpeaks-1)) break
		P2 <- Rnmr1D::peakFiltering(spec, P1, ratioPN)
		if (is.null(P2) || nrow(P2)<(nbpeaks-1)) break
		gn <- as.numeric(rownames(P1)[which( P1$pos %in% P2$pos)])

		# Get peaks with high intensities
		gn <- gn[order(peaks$amp[gn], decreasing=T)][1:(nbpeaks-2)]
		gn <- gn[order(gn)]

		# Remove high peaks on both side
		Amed <- median(peaks$amp[gn])
		while (peaks$amp[gn[1]]>2*Amed) gn <- gn[2:length(gn)]
		while (peaks$amp[gn[length(gn)]]>2*Amed) gn <- gn[1:(length(gn)-1)]

		# Add peaks on both side if peaks are not too high
		k1 <- 0; while(is_snr(min(gn)-k1)) k1 <- k1+1
		k2 <- 0; while(is_snr(max(gn)+k2)) k2 <- k2+1
		aL <- 0.75*median(peaks$amp[gn])
		if (k1>0 && peaks$amp[min(gn)-k1]>aL) k1 <- k1-1
		if (k2>0 && peaks$amp[max(gn)+k2]>aL) k2 <- k2-1
		gn <- (min(gn)-k1):(max(gn)+k2)

		# Add peaks on both side if distHz<Dmin and peaks are not too high
		while (distHz(gn,1,length(gn))<Dmin) {
			k1 <- 1; while(is_snr(min(gn)-k1)) k1 <- k1+1
			k2 <- 1; while(is_snr(max(gn)+k2)) k2 <- k2+1
			J1 <- distHz(c(min(gn)-k1,min(gn)),1,2)
			J2 <- distHz(c(max(gn),max(gn)+k2),1,2)
			d1 <- distHz(c(min(gn)-k1, gn),1,length(gn)+1); a1 <- peaks$amp[min(gn)-k1] 
			d2 <- distHz(c(gn, max(gn)+k2),1,length(gn)+1); a2 <- peaks$amp[max(gn)+k2]
			aL <- 0.95*max(peaks$amp[gn])
			ret <- 0
			if (d1<Dmax && J1<facJ*J && (J1>0.7*J || J1<0.4*J) && a1<aL && a1<1.25*peaks$amp[min(gn)]) { gn <- c(min(gn)-k1, gn); ret <- 1 }
			if (distHz(gn,1,length(gn))>Dmin) break
			if (d2<Dmax && J2<facJ*J && (J2>0.7*J || J2<0.4*J) && a2<aL && a2<1.25*peaks$amp[max(gn)]) { gn <- c(gn, max(gn)+k2); ret <- 1 }
			if (ret==0) break
		}

		# Remove peaks on both side if distHz>Dmax
		while (distHz(gn,1,length(gn))>Dmax) {
			g1 <- c((min(gn)+1):max(gn));  g2 <-c(min(gn):(max(gn)-1)) 
			d1 <- distHz(g1,1,length(g1)); d2 <- distHz(g2,1,length(g2));
			a1 <- peaks$amp[min(gn)+1];    a2 <- peaks$amp[max(gn)-1]
			J1 <- distHz(c(min(gn),min(gn)+1),1,2)
			J2 <- distHz(c(max(gn)-1,max(gn)),1,2)
			if (J1>facJ*J) { gn <- g1;  next }
			if (J2>facJ$J) { gn <- g2;  next }
			if (d1>Dmin && (d1>d2 || a1>a2)) { gn <- g1;  next }
			if (d2>Dmin && (d2>d1 || a2>a1)) { gn <- g2;  next }
			break
		}

		# Remove peaks below the snr on both side
		# so that a peak that satisfies the Dmax constraint can be included next
		if (is_snr(min(gn))) gn <- c((min(gn)+1):max(gn))
		if (is_snr(max(gn))) gn <- c(min(gn):(max(gn)-1)) 

		# In case we missed a peak
		nloop <- 1
		repeat {
			k1 <- k2 <- 1
			if (distHz(gn,1,length(gn))<Dmax) {
				aL <- 0.75*median(peaks$amp[gn])
				a1 <- peaks$amp[min(gn)-k1];
				k1 <- 1; while(is_snr(min(gn)-k1)) k1 <- k1+1
				g1 <- c((min(gn)-k1):max(gn));
				d1 <- distHz(g1,1,length(g1));
				if (d1<Dmax && !is_snr(min(gn)-k1) && a1<aL) { gn <- c((min(gn)-k1):max(gn)) }
			}
			if (distHz(gn,1,length(gn))<Dmax) {
				aL <- 0.75*median(peaks$amp[gn])
				a2 <- peaks$amp[max(gn)+k2];
				k2 <- 1; while(is_snr(max(gn)+k2)) k2 <- k2+1
				g2 <-c(min(gn):(max(gn)+k2)); 
				d2 <- distHz(g2,1,length(g2));
				if (d2<Dmax && !is_snr(max(gn)+k2) && a2<aL) { gn <- c(min(gn):(max(gn)+k2)) }
			}
			if (distHz(gn,1,length(gn))<Dmax) {
				if (k1>1 && distHz(c(min(gn)-1, gn),1,length(gn)+1)<Dmax) { gn <- c(min(gn)-1, gn) }
				if (k2>1 && distHz(c(gn, max(gn)+1),1,length(gn)+1)<Dmax) { gn <- c(gn, max(gn)+1) }
			}
			nloop <- nloop + 1
			if (nloop>2) break
		}

		gn <- (min(gn):max(gn))
		groups <- as.character(gn)
		break
	}
	groups
})

# Merge all cmpds with the same name but having a different index as postfix
internalClass$set("private", "merge_compounds", function(groups)
{
	g2 <- groups[!grepl("*[2-9]", names(groups))]
	for (ng in names(groups)[grepl("*[2-9]", names(groups))])
		g2[[gsub("[2-9]$", "", ng)]] <- c(g2[[gsub("[2-9]$", "", ng)]], groups[[ng]])
	g2
})

# Find a list of compounds based on their pattern
# Ex: compound <- list(
#         'C1'=c('d',3.365,6.7),        # doublet (d), central ppm, J
#         'C2'=c('s',3.35,0.05),        # singulet (s), central ppm, tolerance for ppm
#         'C3'=c('dd',2.89,4.5,16.3),   # doublet of doublet (dd), central ppm, J1, J2
#         'C4'=c('b',1.125,1.132)       # peaks within a range (b), ppm1, ppm2
#     )
# calibration <- c( pzone1, pzone2, p0, type) with type : 1 for d, 0 for s
internalClass$set("private", "find_compounds", function(spec, peaks, compounds, calibration=0)
{
	groups <- NULL
	if (!is.null(compounds) && length(compounds)>0)
	{
		# calibration before if required
		if (length(calibration)>1) {
			V <- as.numeric(calibration[1:3])
			dppm <- get_ppm_shift(spec, V[1], V[2], V[3], type=ifelse(calibration[4]==1,'d','s'))
			peaks$ppm <- peaks$ppm - dppm
		}
		# for each compound, search for peaks matching the requested pattern
		for (k in 1:length(compounds)) {
			cmpd <- names(compounds)[k]
			pattern <- compounds[[k]][1]
			params <- as.numeric(compounds[[k]][-1])
			if (pattern == 's') {
				rank <- ifelse( length(params)>2, params[3], 1)
				groups[[cmpd]] <- tryCatch({ find_pattern_s(spec, peaks, params[1], params[2], rank) },
									error = function(e) { return(NULL) })
				next
			}
			if (pattern == 'd') {
				crit <- ifelse( length(params)>2, params[3], 0)
				if( length(params)>5 ) {
					groups[[cmpd]] <- tryCatch({ find_pattern_d(spec, peaks, params[1], params[2], crit, params[4], params[5], params[6]) },
										error = function(e) { return(NULL) })
				} else if( length(params)>4 ) {
					groups[[cmpd]] <- tryCatch({ find_pattern_d(spec, peaks, params[1], params[2], crit, params[4], params[5]) },
										error = function(e) { return(NULL) })
				} else if( length(params)>3 ) {
					groups[[cmpd]] <- tryCatch({find_pattern_d(spec, peaks, params[1], params[2], crit, params[4]) },
										error = function(e) { return(NULL) })
				} else {
					groups[[cmpd]] <- tryCatch({ find_pattern_d(spec, peaks, params[1], params[2], crit) },
										error = function(e) { return(NULL) })
				}
				next
			}
			if (pattern == 't') {
				crit <- ifelse( length(params)>2, params[3], 0)
				if( length(params)>4 ) {
					groups[[cmpd]] <- tryCatch({ find_pattern_t(spec, peaks, params[1], params[2], crit, params[4], params[5]) },
										error = function(e) { return(NULL) })
				} else if( length(params)>3 ) {
					groups[[cmpd]] <- tryCatch({ find_pattern_t(spec, peaks, params[1], params[2], crit, params[4]) },
										error = function(e) { return(NULL) })
				} else {
					groups[[cmpd]] <- tryCatch({ find_pattern_t(spec, peaks, params[1], params[2], crit) },
										error = function(e) { return(NULL) })
				}
				next
			}
			if (pattern == 'dd') {
				crit <- ifelse( length(params)>3, params[4], 0)
				groups[[cmpd]] <- tryCatch({ find_pattern_dd(spec, peaks, params[1], params[2], params[3], crit) },
									error = function(e) { return(NULL) })
				next
			}
			if (pattern == 'q') {
				Ramp <- ifelse( length(params)>2 & params[3]!=0, params[3], 3)
				groups[[cmpd]] <- tryCatch({ find_pattern_q(spec, peaks, params[1], params[2], Ramp) },
									error = function(e) { return(NULL) })
				next
			}
			if (pattern == 'm') {
				groups[[cmpd]] <- tryCatch({ find_pattern_m(spec, peaks, params[1], params[2]) },
									error = function(e) { return(NULL) })
				next
			}
			if (pattern == 'm2') {
				groups[[cmpd]] <- tryCatch({ find_pattern_m2(spec, peaks, params[1], params[2])},
									error = function(e) { return(NULL) })
				next
			}
			if (pattern == 'b') {
				nbpeaks <- ifelse( length(params)>2, params[3], 0 )
				ratioPN <- ifelse( length(params)>3, params[4], 5 )
				groups[[cmpd]] <- tryCatch({ find_peaks_range(spec, peaks, params[1], params[2], nbpeaks, ratioPN) },
									error = function(e) { return(NULL) })
				next
			}
			if (pattern == 'r2') {
				groups[[cmpd]] <- tryCatch({ find_peaks_rule_r2(spec, peaks, params[1], params[2], params[3], params[4]) },
									error = function(e) { return(NULL) })
				next
			}
			if (pattern == 'r3') {
				groups[[cmpd]] <- tryCatch({ find_peaks_rule_r3(spec, peaks, params[1], params[2], params[3]) },
									error = function(e) { return(NULL) })
				next
			}
			if (pattern == 'r4') {
				groups[[cmpd]] <- tryCatch({ find_peaks_rule_r4(spec, peaks, params[1], params[2], params[3]) },
									error = function(e) { return(NULL) })
				next
			}
			if (pattern == 'r5') {
				groups[[cmpd]] <- tryCatch({ find_peaks_rule_r5(spec, peaks, params[1], params[2], params[3], params[4]) },
									error = function(e) { return(NULL) })
				next
			}
			if (pattern == 'r7') {
				groups[[cmpd]] <- tryCatch({ find_peaks_rule_r7(spec, peaks, params[1], params[2], params[3], params[4]) },
									error = function(e) { return(NULL) })
				next
			}
			if (pattern == 'r9') {
				J        <- ifelse( length(params)>2, params[3], 2 )
				nbpeaks  <- ifelse( length(params)>3, params[4], 5 )
				snrthres <- ifelse( length(params)>4, params[5], 5 )
				groups[[cmpd]] <- tryCatch({ find_peaks_rule_r9(spec, peaks, params[1], params[2], J, nbpeaks, snrthres) }, 
									error = function(e) { return(NULL) })
				next
			}
			next
		}
	}
	if (!is.null(groups))
		groups <- merge_compounds(groups)
	groups
})

# Find the list of compounds in the request zones
# output : a list of compound => (list of peaks)
internalClass$set("private", "find_compounds_by_zone", function(spec, peaks, zones)
{
	groups <- NULL
	if (!is.null(zones) && is.numeric(zones) && zones %in% PROFILE$quantif$zone)
	{
		groups <- list()
		for(k in 1:length(zones)) {
			P <- PROFILE$quantif[PROFILE$quantif$zone==zones[k],,drop=F]
			for(i in 1:nrow(P)) {
				if (grepl(',',P$P3[i])) {
					P3 <- as.numeric(simplify2array(strsplit(P$P3[i],',')))
				} else {
					P3 <- as.numeric(P$P3[i])
				}
				P4 <- 0
				if (grepl(',',P$P4[i])) {
					type <- 's'
					V <- as.vector(simplify2array(strsplit(P$P4[i],',')))
					if (length(V)>3 && V[4] %in% c('s','d')) {
						type <- V[4]
						V[4] <- ifelse(type=='d', 1, 0)
					}
					P4 <- as.numeric(V)
				}
			    L <- find_compounds(spec, peaks, list('C1'=c(P$pattern[i], P$P1[i], P$P2[i], P3)), P4)
				if (!is.null(L))
					if (P$compound[i] %in% names(groups)) {
						groups[[P$compound[i]]] <- c(groups[[P$compound[i]]], L$C1)
					} else {
						groups[[P$compound[i]]] <- L$C1
					}
			}
		}
		if (length(groups)==0) groups <- NULL
	}
	groups
})
