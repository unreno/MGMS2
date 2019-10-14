
#' A Testing Function
#'
#' This function allows you to test your love of functions.
#' @param
#' @keywords
#' @export
#' @examples
#' testfunction()

read.summary.file <- function(files, species){
	pool <- list()
	for (i in 1:length(files)){
		pool[[i]] <- read.csv(file=files[i])
		tmp <-  strsplit(files[i], "/")[[1]]
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
