library(MALDIquant)
library(MALDIquantForeign)

#' A Testing Function
#'
#' This function allows you to test your love of functions.
#' @param
#' @keywords
#' @export
#' @examples
#' testfunction()

simulate.ind.spec.single <- function(interest, mz.tol, species, strain){
	spec <- NULL
	for (i in 1:nrow(interest)){
		current <- interest[i,] 
		log.int <- rnorm(1, current$mean.log.int, current$sd.log.int)
		mz <- rnorm(1, current$mz, mz.tol/qnorm(0.9995)) #0.999 to 0.9995
			tmp <- data.frame(mz=mz, log.int=log.int, species=species, strain=strain, missing.rate=current$missing.rate)
			if (length(spec)<1){
				spec <- tmp
			}else{
				spec <- rbind(spec, tmp)
			}        
		}
	spec$normalized.int <- exp(spec$log.int)
	spec$normalized.int <- spec$normalized.int/max(spec$normalized.int)
	return(spec)
}


