
#	http://r-pkgs.had.co.nz/tests.html


suppressWarnings({
	
	spectra.processed.A <- process_monospectra(
		file=system.file("extdata", "listA.txt", package="MGMS2"),
		mass.range=c(1000,2200))
	spectra.processed.B <- process_monospectra(
		file=system.file("extdata", "listB.txt", package="MGMS2"),
		mass.range=c(1000,2200))
	spectra.processed.C <- process_monospectra(
		file=system.file("extdata", "listC.txt", package="MGMS2"),
		mass.range=c(1000,2200))
	
	spectra.mono.summary.A <- summarize_monospectra(
		processed.obj=spectra.processed.A,
		species='A', directory=tempdir())
	spectra.mono.summary.B <- summarize_monospectra(
		processed.obj=spectra.processed.B,
		species='B', directory=tempdir())
	spectra.mono.summary.C <- summarize_monospectra(
		processed.obj=spectra.processed.C,
		species='C', directory=tempdir())
	
	mono.info=gather_summary(c(spectra.mono.summary.A, spectra.mono.summary.B, spectra.mono.summary.C))
	
	summary=gather_summary_file(directory=tempdir())
	
	mixture.ratio <- list()
	mixture.ratio['A']=1
	mixture.ratio['B']=0.5
	mixture.ratio['C']=0
	
	sim.template <- create_insilico_mixture_template(mono.info)
	
	#insilico.spectrum <- simulate_poly_spectra(sim.template, mixture.ratio)
	
	#insilico.spectrums <- simulate_many_poly_spectra(mono.info, mixture.ratio=mixture.ratio, nsim=10 )
	
})




test_that("under_development",{

})





#library(stringr)
#test_that("str_length is number of characters", {
#  expect_equal(str_length("a"), 1)
#  expect_equal(str_length("ab"), 2)
#  expect_equal(str_length("abc"), 3)
#})


