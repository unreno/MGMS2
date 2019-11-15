
#' @export
species_list <- 	list(
	'Ab' = 'Acinetobacter baumannii',
	'Ec' = 'Enterobacter cloacae',
	'Ef' = 'Enterococus faecalis',
	'Kp' = 'Klebsiella pneumoniae',
	'Pa' = 'Pseudomonas aerugonusa',
	'Sa' = 'Staphylococcus aureus'
)

#' @export
species = names(species_list)





#' testfunction
#'
##' testfunction
##' @keywords cats
#' @export
#' @examples
#' testfunction()

testfunction <- function(){
	message("MGMS2 testfunction")
}

#' filtermass
#'
##' This function allows you to express your love of cats.
##' @param spectra	Add description
##' @param mass.range	Add description
##' @keywords cats
#' @export
##' @examples
##' filtermass()

filtermass <- function(spectra, mass.range){
	for (i in 1:length(spectra)){
		mass = spectra[[i]]@mass
		index=which(mass >= mass.range[1] & mass <= mass.range[2])
		spectra[[i]]@mass = spectra[[i]]@mass[index]
		spectra[[i]]@intensity = spectra[[i]]@intensity[index]
	}
	return(spectra)
}

#' gather_summary
#'
##' This function allows you to express your love of cats.
##' @param x	Add description
##' @keywords cats
#' @export
##' @examples
##' testfunction()

gather_summary <- function(x){
	species <- NULL
	for (i in 1:length(x)){
		species[i] <- as.character(x[[i]]$species[1])
	}
	return(list(pool=x, species=species))
}

#' gather_summary_file
#'
##' This function allows you to express your love of cats.
##' @param directory	Add description
##' @keywords cats
#' @export
##' @examples
##' testfunction()

gather_summary_file <- function(directory){
	files <- list.files(path=directory, pattern=".csv", full.names = TRUE)
	pool <- list(); species <- NULL
	for (i in 1:length(files)){
		pool[[i]] <- read.csv(file=files[i])
		species[i] <- as.character(pool[[i]]$species[1])
	}
	return(list(pool=pool, species=species))
}


#' preprocessMS
#'
##' This function allows you to express your love of cats.
##' @param spectra	Add description
##' @param halfWindowSize	Add description
##' @param SNIP.iteration	Add description
##' @keywords cats
#' @export
##' @examples
##' testfunction()

preprocessMS <- function(spectra, halfWindowSize=20, SNIP.iteration=60){
	spectra <- transformIntensity(spectra, method="sqrt")
	spectra <- smoothIntensity(spectra, method="SavitzkyGolay", halfWindowSize=halfWindowSize)
	spectra <- removeBaseline(spectra, method="SNIP", iterations=SNIP.iteration)
	spectra <- calibrateIntensity(spectra,method="TIC") #intensity calibration/normalization
	return(spectra)
}

#' summarize_monospectra
#'
##' This function allows you to express your love of cats.
##' @param processed.obj	Add description
##' @param species	Add description
##' @param directory	Add description
##' @param minFrequency	Add description
##' @param align.tolerance	Add description
##' @param snr	Add description
##' @param halfWindowSize	Add description
##' @param top.N	Add description
##' @keywords cats
#' @export
##' @examples
##' testfunction()

summarize_monospectra <- function(processed.obj, species, directory, minFrequency=0.50, align.tolerance=0.0005, snr=3, halfWindowSize=20, top.N=50){
	spec.summary.comb <- list()
	for (i in processed.obj$strain.unique){
		index <- which(processed.obj$strain.no==processed.obj$strain.unique[i])
		spectra.interest <- list()
		for (j in 1:length(index)){
			spectra.interest[[j]] <- processed.obj$spectra[[index[j]]]
		}
		spec.summary=summary_mono(spectra.interest, minFrequency, align.tolerance, snr, halfWindowSize, 	top.N)
		strain = processed.obj$strain.name[processed.obj$strain.no==i][1]
		spec.summary.comb[[i]] = cbind(spec.summary, species, strain)
		write.csv(spec.summary.comb[[i]], file=paste(directory, "/", species, "_",strain , ".csv", sep=""), row.names = FALSE)
	}
	return(spec.summary.comb)
}


#' process_monospectra
#'
##' This function allows you to express your love of cats.
##' @param file	Add description
##' @param mass.range	Add description
##' @param halfWindowSize	Add description
##' @param SNIP.iteration	Add description
##' @keywords cats
#' @export
##' @examples
##' testfunction()

process_monospectra <- function(file, mass.range=c(1000,2200), halfWindowSize=20, SNIP.iteration=60){
	file.interest <- read.csv(file=file, "\t", header=TRUE)

	#	These files all contained the absolute path.
	#	This won't work for most users.
	#	I trimmed the path from the content and will prepend the path of this file.
	dir=dirname(file)
	spectra <- importMzXml(paste(dir,file.interest$file.name,sep='/'))
	#spectra <- importMzXml(paste(file.interest$file.name))

	spectra <- filtermass(spectra, mass.range)
	spectra <- preprocessMS(spectra, halfWindowSize=20, SNIP.iteration=60)
	strain.no <- file.interest$strain.no
	strain.unique <- unique(strain.no)
	strain.name <- file.interest$strain
	return(list(spectra=spectra, strain.no=strain.no, strain.unique=strain.unique, strain.name=strain.name))
}


#' summary_mono
#'
##' This function allows you to express your love of cats.
##' @param spectra.interest	Add description
##' @param minFrequency	Add description
##' @param align.tolerance	Add description
##' @param snr	Add description
##' @param halfWindowSize	Add description
##' @param top.N	Add description
##' @keywords cats
#' @export
##' @examples
##' testfunction()

summary_mono <- function(spectra.interest, minFrequency=0.50, align.tolerance=0.0005, snr=3, halfWindowSize=20, top.N=50){
	ref.peak <- detectPeaks(spectra.interest, method="MAD", halfWindowSize = halfWindowSize, SNR=snr)
	reference.peaks <- referencePeaks(ref.peak, method="relaxed",
		minFrequency = minFrequency, tolerance=align.tolerance)
	spectra.interest2 <- alignSpectra(spectra.interest, reference=reference.peaks,
		halfWindowSize = halfWindowSize,
		SNR=snr, tolerance=align.tolerance, warpingMethod = "lowess")
	peaks2 <- detectPeaks(spectra.interest2, SNR=snr, halfWindowSize=halfWindowSize,
		method="MAD")
	peaks3 <- binPeaks(peaks2, tolerance=align.tolerance)
	peaks4 <- filterPeaks(peaks3, minFrequency=minFrequency)
	featureMatrix <- intensityMatrix(peaks4)

	mean.log.int = sd.log.int = missing.rate = NULL
	for (k in 1:ncol(featureMatrix)){
		featureMatrix[,k] = ifelse(featureMatrix[,k]==0, NA, log(featureMatrix[,k]))
		y <- featureMatrix[,k]
		mean.log.int[k] = mean(y, na.rm=TRUE)
		sd.log.int[k] = sd(y, na.rm=TRUE)
		missing.rate[k] = mean(is.na(y))
	}
	summary.info <- data.frame(
		mz = as.numeric(colnames(featureMatrix)),
		mean.log.int = mean.log.int,
		sd.log.int = sd.log.int,
		missing.rate = missing.rate)
	th <- sort(summary.info$mean.log.int, decreasing=TRUE)[top.N][1]
	if (!is.na(th)){
		summary.info <- summary.info[summary.info$mean.log.int >= th, ]
	}
	return(summary.info)
}

#' read_summary_file
#'
##' This function allows you to express your love of cats.
##' @param files	Add description
##' @keywords cats
#' @export
##' @examples
##' testfunction()

read_summary_file <- function(files){
	pool <- list()
	for (i in 1:length(files)){
		pool[[i]] <- read.csv(file=files[i])
		tmp <-	strsplit(files[i], "/")[[1]]
		tmp <- strsplit(tmp[length(tmp)], "[.]")[[1]][1]
		tmp <- strsplit(tmp, "_")[[1]]
		index=which(species==tmp[1])
		if (length(index)>0){
			pool[[i]]$species <- tmp[1]
			pool[[i]]$strain <- tmp[2]
		}
	}
	return(pool)
}

#' simulate_ind_spec_single
#'
##' This function allows you to express your love of cats.
##' @param interest	Add description
##' @param mz.tol	Add description
##' @param species	Add description
##' @param strain	Add description
##' @keywords cats
#' @export
##' @examples
##' testfunction()

simulate_ind_spec_single <- function(interest, mz.tol, species, strain){
	spec <- NULL
	for (i in 1:nrow(interest)){
		current <- interest[i,]
		log.int <- rnorm(1, current$mean.log.int, current$sd.log.int)
		mz <- rnorm(1, current$mz, mz.tol/qnorm(0.9995)) #Change
		tmp <- data.frame(mz=mz, log.int=log.int, species=species, strain=strain, missing.rate=current$missing.rate)
		if (length(spec)<1){
			spec <- tmp
		}else{
			spec <- rbind(spec, tmp)
		}
	}
	spec$normalized.int <- exp(spec$log.int)
	return(spec)
}

#' build_bin_sim_spec
#'
##' This function allows you to express your love of cats.
##' @param bin.mass	Add description
##' @param target	Add description
##' @param bin.size	Add description
##' @keywords cats
#' @export
##' @examples
##' testfunction()

build_bin_sim_spec <- function(bin.mass, target, bin.size=1){
	bin.sim = chosen.mz = mid.mass = array(c(0), dim=length(bin.mass))
	for (j in 1:length(bin.mass)){
		current.min.mass = bin.mass[j]
		current.max.mass = current.min.mass+bin.size
		mid.mass[j] = (current.min.mass+current.max.mass)/2
		index = which((target$mz >= current.min.mass) & (target$mz < current.max.mass))
		if (length(index)>0){
			index2 = which.max(target$normalized.int[index])
			bin.sim[j] = max(target$normalized.int[index])
			chosen.mz[j] = target$mz[index][index2]
		}
	}
	binned.target = as.data.frame(cbind(mid.mass, chosen.mz, bin.sim))
	return(binned.target)
}

#' simulate_poly_spectra
#'
##' This function allows you to express your love of cats.
##' @param sim.template	Add description
##' @param mixture.ratio	Add description
##' @param spectrum.name	Add description
##' @param mixture.missing.prob.peak	Add description
##' @param noise.peak.ratio	Add description
##' @param snr.basepeak	Add description
##' @param noise.cv	Add description
##' @param mz.range	Add description
##' @keywords cats
#' @export
##' @examples
##' testfunction()

simulate_poly_spectra <- function(sim.template, mixture.ratio, spectrum.name='Spectrum',
		mixture.missing.prob.peak = 0.05, noise.peak.ratio = 0.05, snr.basepeak = 500 , noise.cv = 0.25, mz.range=c(1000,2200)){

	#mixture ratio
	mixture.ratio.list <- NULL
	for (i in 1:nrow(sim.template)){
		mixture.ratio.list[i]=unlist(mixture.ratio[as.character(sim.template$species[i])])
	}
	sim.template$normalized.int = sim.template$normalized.int *	mixture.ratio.list
	sim.template = sim.template[sim.template$normalized.int>0,]
	sim.template = sim.template[sim.template$mz >= mz.range[1] & sim.template$mz <= mz.range[2], ]

	#mixture missing prob
	missing.probability <- mixture.missing.prob.peak+sim.template$missing.rate
	missing.probability <- ifelse(missing.probability<0, 0, ifelse(missing.probability>1, 1, missing.probability))
	missing <- rbinom(nrow(sim.template), 1, missing.probability)
	if (length(missing)>0){
		sim.template <- sim.template[missing==0,]
	}else{
		stop("too many missing.\n")
	}
	
	#add noise peaks
	if (noise.peak.ratio>0){
		noise.mean.int <- 1/snr.basepeak
		noise.n <- round(nrow(sim.template)*noise.peak.ratio,0)
		noise.int <- exp(rnorm(noise.n, log(noise.mean.int), log(noise.cv+1)))
		noise.mz <- runif(noise.n, mz.range[1], mz.range[2])
		noise.temp <- data.frame(mz=noise.mz, log.int=NA, species=c("Noise"), strain=c("Noise"), normalized.int=noise.int, missing.rate=0)
		sim.template = rbind(sim.template, noise.temp)
	}
	
	#normalize by the highest peak
	sim.template$normalized.int = sim.template$normalized.int / max(sim.template$normalized.int)
	

	#species.interest= names(mixture.ratio)[mixture.ratio>0]
	#	Keep all species
	species.interest= names(mixture.ratio)


#change the palette to have purple instead of yellow
col.pal <- palette()
col.pal[7] <- "purple"
palette(col.pal)


	for (p in species.interest){
		if (p==species.interest[1]){
			plot(sim.template$mz[sim.template$species==p],
				sim.template$normalized.int[sim.template$species==p],
				ylim=c(0,1.1), xlim=mz.range, col=2, type="h", main=spectrum.name,
				cex.axis=1.5, cex.lab=1.5,
				xlab="m/z", ylab="Relative intensity")
			legend("topleft",
				c(species.interest, "Noise"),
				lty=rep(1, length(species.interest)+1),
				cex=1.5,
				col=c((1+1):(length(species.interest)+1),1)
			)
		}else{
			points(sim.template$mz[sim.template$species==p],
				sim.template$normalized.int[sim.template$species==p], col=which(species.interest==p)+1, type="h")
		}
		if (length(sim.template$species=="Noise")>0){
			points(sim.template$mz[sim.template$species=="Noise"],
			sim.template$normalized.int[sim.template$species=="Noise"], col='black', type="h")
		}
	}
	return(sim.template=sim.template[,c(1,3,4,6)])
}


#' simulate_many_poly_spectra
#'
##' This function allows you to express your love of cats.
##' @param mono.info	Add description
##' @param nsim	Add description
##' @param file	Add description
##' @param mixture.ratio	Add description
##' @param mixture.missing.prob.peak	Add description
##' @param noise.peak.ratio	Add description
##' @param snr.basepeak	Add description
##' @param noise.cv	Add description
##' @param mz.range	Add description
##' @param mz.tol	Add description
##' @keywords cats
#' @export
##' @examples
##' testfunction()

simulate_many_poly_spectra <- function(mono.info, nsim=10000, file='MGMS2_insilico_spectra.pdf', mixture.ratio,
		mixture.missing.prob.peak = 0.05, noise.peak.ratio = 0.05,
		snr.basepeak = 500 , noise.cv = 0.25, mz.range=c(1000,2200), mz.tol=0.5) {
	sim.spectra.collection <- NULL
	pdf(file)
	for (i in 1:nsim){
		sim.template <- create_insilico_mixture_template(mono.info, mz.tol=0.5) #simulated polyspectra (users can modify this)
		sim.spectra.collection[[i]]=simulate_poly_spectra(sim.template, mixture.ratio=mixture.ratio,
			spectrum.name=paste('Spectrum', i, sep="_"),
			mixture.missing.prob.peak =mixture.missing.prob.peak,
			noise.peak.ratio = noise.peak.ratio,
			snr.basepeak = snr.basepeak , noise.cv =noise.cv, mz.range=mz.range)
	}
	dev.off()
	return(sim.spectra.collection)
}


#' create_insilico_mixture_template
#'
##' This function allows you to express your love of cats.
##' @param mono.info	Add description
##' @param mz.tol	Add description
##' @keywords cats
#' @export
##' @examples
##' create_insilico_mixture_template()

create_insilico_mixture_template <- function(mono.info, mz.tol=0.5){
	spec.mixture <- NULL
	unique.species <- unique(mono.info$species)
	for (i in 1:length(unique.species)){
		index=which(mono.info$species==unique.species[i])
		if (length(index)>0){
			if (length(index)==1){j=index}else{j=sample(index)[1]}
			ind.species <- simulate_ind_spec_single(mono.info$pool[[j]], mz.tol=mz.tol, species=mono.info$species[j], mono.info$pool[[j]]$strain[1])
			spec.mixture = rbind(spec.mixture, ind.species)
		}
	}
	return(spec.mixture)
}

