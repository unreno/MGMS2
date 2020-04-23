
#' @import utils
#' @importFrom grDevices dev.off palette pdf
#' @importFrom graphics legend
#' @importFrom stats qnorm rbinom rnorm runif sd
#' @import MALDIquant
#' @import MALDIquantForeign
NULL


utils::globalVariables(c("species"))

#' filtermass
#'
#' Internal function. This function removes peaks with their mass values (m/z values) outside a given mass range.
#' This function is used in \code{\link{process_monospectra}}.
#' @param spectra Mass Spectra (A MALDIquant MassSpectrum (S4) object). An output of \code{\link[MALDIquantForeign]{importMzXml}}.
#' @param mass.range Mass (m/z) range (a vector). For exmaple, c(1000,2200).
#' @return A list of filtered mass spectra (MALDIquant MassSpectrum (S4) objects) which contains mass, intensity, and metaData.

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
#' This function combines outputs from \code{\link{summarize_monospectra}}.
#' @param x A list of multiple monomicrobial mass spectra information from \code{\link{summarize_monospectra}}.
#' @return A list of combined summaries (data frames) of mass spectra from \code{\link{summarize_monospectra}} and the corresponding species (a vector).
#' @examples
#' spectra.processed.A <- process_monospectra(
#'    file=system.file("extdata", "listA.txt", package="MGMS2"),
#'    mass.range=c(1000,2200))
#' spectra.processed.B <- process_monospectra(
#'    file=system.file("extdata", "listB.txt", package="MGMS2"),
#'    mass.range=c(1000,2200))
#' spectra.processed.C <- process_monospectra(
#'    file=system.file("extdata", "listC.txt", package="MGMS2"),
#'    mass.range=c(1000,2200))
#' spectra.mono.summary.A <- summarize_monospectra(
#'    processed.obj=spectra.processed.A,
#'    species='A', directory=tempdir())
#' spectra.mono.summary.B <- summarize_monospectra(
#'    processed.obj=spectra.processed.B,
#'    species='B', directory=tempdir())
#' spectra.mono.summary.C <- summarize_monospectra(
#'    processed.obj=spectra.processed.C,
#'    species='C', directory=tempdir())
#' mono.info=gather_summary(c(spectra.mono.summary.A, spectra.mono.summary.B, spectra.mono.summary.C))
#' @export

gather_summary <- function(x){
	species <- NULL
	for (i in 1:length(x)){
		species[i] <- as.character(x[[i]]$species[1])
	}
	return(list(pool=x, species=species))
}

#' gather_summary_file
#'
#' This function combines output files from \code{\link{summarize_monospectra}}.
#' @param directory A directory that contains summary files from \code{\link{summarize_monospectra}}.
#' @return A list of combined summaries of mass spectra (data frames) from \code{\link{summarize_monospectra}} and the corresponding species (a vector).
#' @examples
#' spectra.processed.A <- process_monospectra(
#'    file=system.file("extdata", "listA.txt", package="MGMS2"),
#'    mass.range=c(1000,2200))
#' spectra.processed.B <- process_monospectra(
#'    file=system.file("extdata", "listB.txt", package="MGMS2"),
#'    mass.range=c(1000,2200))
#' spectra.processed.C <- process_monospectra(
#'    file=system.file("extdata", "listC.txt", package="MGMS2"),
#'    mass.range=c(1000,2200))
#' spectra.mono.summary.A <- summarize_monospectra(
#'    processed.obj=spectra.processed.A,
#'    species='A', directory=tempdir())
#' spectra.mono.summary.B <- summarize_monospectra(
#'    processed.obj=spectra.processed.B,
#'    species='B', directory=tempdir())
#' spectra.mono.summary.C <- summarize_monospectra(
#'    processed.obj=spectra.processed.C,
#'    species='C', directory=tempdir())
#' summary <- gather_summary_file(directory=tempdir())
#' @export

gather_summary_file <- function(directory){
	files <- list.files(path=directory, pattern=".csv", full.names = TRUE)
	pool <- list()
	species <- NULL
	for (i in 1:length(files)){
		pool[[i]] <- read.csv(file=files[i])
		species[i] <- as.character(pool[[i]]$species[1])
	}
	return(list(pool=pool, species=species))
}


#' preprocessMS
#'
#' Internal function. This function preprocesses spectra by transforming/smoothing intensity, removing baseline, and calibrating intensities.
#' @param spectra Spectra. A MALDIquant object. An output of either \code{\link[MALDIquantForeign]{importMzXml}} or \code{\link{filtermass}}.
#' @param halfWindowSize halfWindowSize The highest peaks in the given window (+/-halfWindowSize) will be recognized as peaks. (Default: 20). See \code{\link[MALDIquant]{detectPeaks}} for details.
#' @param SNIP.iteration SNIP.iteration An iteration used to remove the baseline of an spectrum. (Default: 60). See \code{\link[MALDIquant]{removeBaseline}} for details.
#' @return The processed mass spectra. A list of MALDIquant MassSpectrum objects (S4 objects).

preprocessMS <- function(spectra, halfWindowSize=20, SNIP.iteration=60){
	spectra <- transformIntensity(spectra, method="sqrt")
	spectra <- smoothIntensity(spectra, method="SavitzkyGolay", halfWindowSize=halfWindowSize)
	spectra <- removeBaseline(spectra, method="SNIP", iterations=SNIP.iteration)
	spectra <- calibrateIntensity(spectra,method="TIC") #intensity calibration/normalization
	return(spectra)
}

#' summarize_monospectra
#'
#' This function summarizes monomicrobial spectra and writes summary in the specified directory.
#' @param processed.obj A list from \code{\link{process_monospectra}} which contains peaks information for each strain.
#' @param species Species name.
#' @param directory Directory. (By default, no summary file will be generated.)
#' @param minFrequency Percentage value. A minimum occurrence proportion required for building a reference peaks. All peaks with their occurence proportion less than minFrequency will be moved. (Default: 0.50). See \code{\link[MALDIquant]{filterPeaks}} and \code{\link[MALDIquant]{referencePeaks}} for details.
#' @param align.tolerance Mass tolerance. Must be multiplied by 10^-6 for ppm. (Default: 0.0005).
#' @param snr Signal-to-noise ratio. (Default: 3).
#' @param halfWindowSize The highest peaks in the given window (+/-halfWindowSize) will be recognized as peaks. (Default: 20). See \code{\link[MALDIquant]{detectPeaks}} for details.
#' @param top.N The top N peaks will be chosen for the analysis. An integer value. (Default: 50).
#' @return A data frame that contains the peaks informations: m/z, mean log intensity, standard deviation of log intensity, missing rate of peaks. In addition, it also contains species and strain information.
#' @examples
#' spectra.processed.A <- process_monospectra(
#'    file=system.file("extdata", "listA.txt", package="MGMS2"),
#'    mass.range=c(1000,2200))
#' spectra.mono.summary.A <- summarize_monospectra(
#'    processed.obj=spectra.processed.A, species='A',
#'    directory=tempdir())
#' @export

summarize_monospectra <- function(processed.obj, species, directory=NULL, minFrequency=0.50, align.tolerance=0.0005, snr=3, halfWindowSize=20, top.N=50){
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
		if (length(directory)>0){
		write.csv(spec.summary.comb[[i]], file=paste(directory, "/", species, "_",strain , ".csv", sep=""), row.names = FALSE)
		}
	}
	return(spec.summary.comb)
}


#' process_monospectra
#'
#' This function processes multiple mzXML files which are listed in the file that an user specifies.
#' @param file A file name. This file is a tab-delimited file which contains the following columns: file names, strain.no, and strain. See below for details.
#' @param mass.range The m/z range that users want to consider for the analysis. (Default: c(1000,2200)).
#' @param halfWindowSize A half window size used for the smoothing the intensity values. (Default: 20). See \code{\link[MALDIquant]{smoothIntensity}} for details.
#' @param SNIP.iteration An iteration used to remove the baseline of an spectrum. (Default: 60). See \code{\link[MALDIquant]{removeBaseline}} for details.
#' @return A list of processed monobacterial mass spectra (S4 objects, MALDIquant MassSpectrum objects), and their strain numbers (a vector), unique strains (a vector), and strain names (a vector).
#' @examples
#' spectra.processed.A <- process_monospectra(
#'    file=system.file("extdata", "listA.txt", package="MGMS2"),
#'    mass.range=c(1000,2200))
#' @export

process_monospectra <- function(file, mass.range=c(1000,2200), halfWindowSize=20, SNIP.iteration=60){
	file.interest <- read.csv(file=file, "\t", header=TRUE)
	dir=dirname(file)
	spectra <- importMzXml(paste(dir,file.interest$file.name,sep='/'))

	spectra <- filtermass(spectra, mass.range)
	spectra <- preprocessMS(spectra, halfWindowSize=20, SNIP.iteration=60)
	strain.no <- file.interest$strain.no
	strain.unique <- unique(strain.no)
	strain.name <- file.interest$strain
	return(list(spectra=spectra, strain.no=strain.no, strain.unique=strain.unique, strain.name=strain.name))
}


#' summary_mono
#'
#' Internal function. This function calculates summary statistics for peaks afterling aligning spectra of interest.
#' @param spectra.interest A list which contains peaks information for a strain of interest.
#' @param minFrequency Percentage value. A minimum occurrence proportion required for building a reference peaks. All peaks with their occurence proportion less than minFrequency will be moved. (Default: 0.50). See \code{\link[MALDIquant]{filterPeaks}} and \code{\link[MALDIquant]{referencePeaks}} for details.
#' @param align.tolerance Mass tolerance. Must be multiplied by 10^-6 for ppm. (Default: 0.0005).
#' @param snr Signal-to-noise ratio. (Default: 3).
#' @param halfWindowSize The highest peaks in the given window (+/-halfWindowSize) will be recognized as peaks. (Default: 20). See \code{\link[MALDIquant]{detectPeaks}} for details.
#' @param top.N The top N peaks will be chosen for the analysis. An integer value. (Default: 50).
#' @return Summary information (Data frame) of spectra of interest.

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


#' simulate_ind_spec_single
#'
#' Internal function. The function simulates m/z and intensity values using given summary statistics.
#' @param interest Summary statistics of spectra.
#' @param mz.tol The tolerance of m/z. This is used to generate m/z values of peaks.
#' @param species Species.
#' @param strain Strain name.
#' @return A data frame that contains m/z, (normalized) intensity values, missing rates of peaks, species name, and strain name.

simulate_ind_spec_single <- function(interest, mz.tol, species, strain){
	spec <- NULL
	for (i in 1:nrow(interest)){
		current <- interest[i,]
		log.int <- rnorm(1, current$mean.log.int, current$sd.log.int)
		mz <- rnorm(1, current$mz, mz.tol/qnorm(0.9995))
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

#' simulate_poly_spectra
#'
#' This function takes simulated m/z and intensities of peaks from \code{\link{create_insilico_mixture_template}} and modifies them based on given parameters.
#' @param sim.template A data frame which contains m/z, log intensitiy, normalized intensity values and missing rates of peaks. There are also species and strain information. An object of \code{\link{create_insilico_mixture_template}}.
#' @param mixture.ratio A list of bacterial mixture ratios for given bacterial species in sim.template.
#' @param spectrum.name A character. An user can define the spectrum name. (Default: 'Spectrum').
#' @param mixture.missing.prob.peak A real value. The missing probability caused by mixing multiple bacteria species. (Default: 0.05)
#' @param noise.peak.ratio A ratio between the numbers of noise and signal peaks. (Default: 0.05)
#' @param snr.basepeak A (base peak) signal to noise ratio. (Default: 500)
#' @param noise.cv A coefficient of variation of noise peaks. (Default: 0.25)
#' @param mz.range A range of m/z values. (Default: c(1000,2200))
#' @return A data frame that contains m/z values of peaks, normalized intensities of peaks, species names, and strain names. A modified version of \code{sim.template}.
#' @examples
#' spectra.processed.A <- process_monospectra(
#'    file=system.file("extdata", "listA.txt", package="MGMS2"),
#'    mass.range=c(1000,2200))
#' spectra.processed.B <- process_monospectra(
#'    file=system.file("extdata", "listB.txt", package="MGMS2"),
#'    mass.range=c(1000,2200))
#' spectra.processed.C <- process_monospectra(
#'    file=system.file("extdata", "listC.txt", package="MGMS2"),
#'    mass.range=c(1000,2200))
#' spectra.mono.summary.A <- summarize_monospectra(
#'    processed.obj=spectra.processed.A,
#'    species='A', directory=tempdir())
#' spectra.mono.summary.B <- summarize_monospectra(
#'    processed.obj=spectra.processed.B,
#'    species='B', directory=tempdir())
#' spectra.mono.summary.C <- summarize_monospectra(
#'    processed.obj=spectra.processed.C,
#'    species='C', directory=tempdir())
#' mono.info=gather_summary(c(spectra.mono.summary.A, spectra.mono.summary.B, spectra.mono.summary.C))
#' mixture.ratio <- list()
#' mixture.ratio['A']=1
#' mixture.ratio['B']=0.5
#' mixture.ratio['C']=0
#' sim.template <- create_insilico_mixture_template(mono.info)
#' insilico.spectrum <- simulate_poly_spectra(sim.template, mixture.ratio)
#' @export

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
	return(sim.template=sim.template[,c(1,6,3,4)])
}


#' simulate_many_poly_spectra
#'
#' The function creates simulated mass spectra in pdf file and returns simulated mass spectra (m/z and intensity values of peaks).
#' @param mono.info A list output of \code{\link{gather_summary}} or {\code{\link{gather_summary_file}}}.
#' @param nsim The number of simulated spectra. (Default: 10000)
#' @param file An output file name. (By default, file=NULL. No pdf file will be generated.)
#' @param mixture.ratio A list of bacterial mixture ratios for given bacterial species in sim.template.
#' @param mixture.missing.prob.peak A real value. The missing probability caused by mixing multiple bacteria species. (Default: 0.05)
#' @param noise.peak.ratio A ratio between the numbers of noise and signal peaks. (Default: 0.05)
#' @param snr.basepeak A (base peak) signal to noise ratio. (Default: 5000)
#' @param noise.cv A coefficient of variation of noise peaks. (Default: 0.25)
#' @param mz.range A range of m/z values. (Default: c(1000,2200))
#' @param mz.tol m/z tolerance. (Default: 0.5)
#' @return A list of data frames. A list of simulated mass spectra (data frames) that contains m/z values of peaks, normalized intensities of peaks, species names, and strain names. This function also creates pdf files which contain simulated spectra.
#' @examples
#' spectra.processed.A <- process_monospectra(
#'    file=system.file("extdata", "listA.txt", package="MGMS2"),
#'    mass.range=c(1000,2200))
#' spectra.processed.B <- process_monospectra(
#'    file=system.file("extdata", "listB.txt", package="MGMS2"),
#'    mass.range=c(1000,2200))
#' spectra.processed.C <- process_monospectra(
#'    file=system.file("extdata", "listC.txt", package="MGMS2"),
#'    mass.range=c(1000,2200))
#' spectra.mono.summary.A <- summarize_monospectra(
#'    processed.obj=spectra.processed.A,
#'    species='A', directory=tempdir())
#' spectra.mono.summary.B <- summarize_monospectra(
#'    processed.obj=spectra.processed.B,
#'    species='B', directory=tempdir())
#' spectra.mono.summary.C <- summarize_monospectra(
#'    processed.obj=spectra.processed.C,
#'    species='C', directory=tempdir())
#' mono.info=gather_summary(c(spectra.mono.summary.A, spectra.mono.summary.B, spectra.mono.summary.C))
#' mixture.ratio <- list()
#' mixture.ratio['A']=1
#' mixture.ratio['B']=0.5
#' mixture.ratio['C']=0
#' insilico.spectra <- simulate_many_poly_spectra(mono.info, mixture.ratio=mixture.ratio, nsim=10)
#' @export

simulate_many_poly_spectra <- function(mono.info, nsim=10000, file=NULL, mixture.ratio,
		mixture.missing.prob.peak = 0.05, noise.peak.ratio = 0.05,
		snr.basepeak = 500 , noise.cv = 0.25, mz.range=c(1000,2200), mz.tol=0.5) {
	sim.spectra.collection <- NULL
	if (length(file)>0){pdf(file)}
	for (i in 1:nsim){
		sim.template <- create_insilico_mixture_template(mono.info, mz.tol=0.5) #simulated polyspectra (users can modify this)
		sim.spectra.collection[[i]]=simulate_poly_spectra(sim.template, mixture.ratio=mixture.ratio,
			spectrum.name=paste('Spectrum', i, sep="_"),
			mixture.missing.prob.peak =mixture.missing.prob.peak,
			noise.peak.ratio = noise.peak.ratio,
			snr.basepeak = snr.basepeak , noise.cv =noise.cv, mz.range=mz.range)
	}
	if (length(file)>0){dev.off()}
	return(sim.spectra.collection)
}


#' create_insilico_mixture_template
#'
#' This function generates an intial template for simulated mass spectra.
#' @param mono.info An output of \code{\link{gather_summary}}.
#' @param mz.tol A m/z tolerance in Da. (Default: 0.5)
#' @return A data frame which contains simulated m/z, log intensity, and normalized intensity values of peaks.
#' @examples
#' spectra.processed.A <- process_monospectra(
#'    file=system.file("extdata", "listA.txt", package="MGMS2"),
#'    mass.range=c(1000,2200))
#' spectra.processed.B <- process_monospectra(
#'    file=system.file("extdata", "listB.txt", package="MGMS2"),
#'    mass.range=c(1000,2200))
#' spectra.processed.C <- process_monospectra(
#'    file=system.file("extdata", "listC.txt", package="MGMS2"),
#'    mass.range=c(1000,2200))
#' spectra.mono.summary.A <- summarize_monospectra(
#'    processed.obj=spectra.processed.A,
#'    species='A', directory=tempdir())
#' spectra.mono.summary.B <- summarize_monospectra(
#'    processed.obj=spectra.processed.B,
#'    species='B', directory=tempdir())
#' spectra.mono.summary.C <- summarize_monospectra(
#'    processed.obj=spectra.processed.C,
#'    species='C', directory=tempdir())
#' mono.info=gather_summary(c(spectra.mono.summary.A, spectra.mono.summary.B, spectra.mono.summary.C))
#' template <- create_insilico_mixture_template(mono.info)
#' @export

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

#' characterize_peak
#'
#' This function characterizes peaks by species/strain in a simulated spectrum after taking the highest peak or merging peaks in each bin.
#' @param spec A data frame that contains m/z values of peaks, normalized intensities of peaks, species names, and strain names. Either an output of \code{\link{simulate_poly_spectra}} or one elements of a list output from \code{\link{simulate_many_poly_spectra}}. 
#' @param option An option on how to merge peaks. There are two options: 1) no merge, thus take the highest intensity peak in each bin after binning a spectrum by bin.size, or 2) take a sum of intensity within each bin after binning a spectrum by bin.size.
#' @return A data frame that contains m/z values of peaks (mz), intensities of peaks (int), species names (species), and strain names (strain). Species and strain columns may contain more than one species/strain if an option 2 is chosen.
#' @examples
#' spectra.processed.A <- process_monospectra(
#'    file=system.file("extdata", "listA.txt", package="MGMS2"),
#'    mass.range=c(1000,2200))
#' spectra.processed.B <- process_monospectra(
#'    file=system.file("extdata", "listB.txt", package="MGMS2"),
#'    mass.range=c(1000,2200))
#' spectra.processed.C <- process_monospectra(
#'    file=system.file("extdata", "listC.txt", package="MGMS2"),
#'    mass.range=c(1000,2200))
#' spectra.mono.summary.A <- summarize_monospectra(
#'    processed.obj=spectra.processed.A,
#'    species='A', directory=tempdir())
#' spectra.mono.summary.B <- summarize_monospectra(
#'    processed.obj=spectra.processed.B,
#'    species='B', directory=tempdir())
#' spectra.mono.summary.C <- summarize_monospectra(
#'    processed.obj=spectra.processed.C,
#'    species='C', directory=tempdir())
#' mono.info=gather_summary(c(spectra.mono.summary.A, spectra.mono.summary.B, spectra.mono.summary.C))
#' mixture.ratio <- list()
#' mixture.ratio['A']=1
#' mixture.ratio['B']=0.5
#' mixture.ratio['C']=0
#' sim.template <- create_insilico_mixture_template(mono.info)
#' insilico.spectrum <- simulate_poly_spectra(sim.template, mixture.ratio)
#' merged.spectrum <- characterize_peak(insilico.spectrum, option=2) 
#' @export
characterize_peak <- function(spec, option=1, bin.size=1, min.mz=1000, max.mz=2200){	
	bin.mass = seq(min.mz, max.mz, by=bin.size)
	bin.mass = bin.mass[-length(bin.mass)]
	chosen.mz = chosen.int = array(c(0), dim=length(bin.mass)) 
	chosen.species = chosen.strain = array(c(NA), dim=length(bin.mass))
	for (j in 1:length(bin.mass)){
		current.min.mass = bin.mass[j]
		current.max.mass = current.min.mass+bin.size
		index = which((spec$mz >= current.min.mass) & (spec$mz < current.max.mass))
		if (option==1){
			if (length(index)>0){
				index2 = which.max(spec$normalized.int[index])
				chosen.int[j] = max(spec$normalized.int[index])
				chosen.mz[j] = spec$mz[index][index2]
				chosen.species[j] = as.character(spec$species[index][index2])
				chosen.strain[j] = as.character(spec$strain[index][index2])
			}
		}else{
			if (length(index)>0){
				chosen.int[j] = sum(spec$normalized.int[index])
				chosen.mz[j] = mean(spec$mz[index])
				chosen.species[j] = paste(unique(as.character(spec$species[index])), collapse=';')
				chosen.strain[j] = paste(unique(as.character(spec$strain[index])), collapse=';')
			}		
		}
	}
	binned.spec = data.frame(mz=chosen.mz, 
							int=chosen.int,
							species=chosen.species,
							strain=chosen.strain)
	binned.spec = binned.spec[binned.spec$mz>0,]
	binned.spec = binned.spec[sort.int(binned.spec$mz, index.return=TRUE)$ix,]
	return(binned.spec)
}
	

