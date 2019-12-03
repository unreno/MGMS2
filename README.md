# MGMS2

Library for the Interactive Poly-Bacterial in-Silico Glycolipid Spectrum Simulator


##	User


###	Installation


####	Production Version

```R
install.packages(c("MALDIquant","MALDIquantForeign"),type="source")
install.packages("MGMS2")
library("MGMS2")
```
or
```BASH
R -e 'install.packages("MGMS2")'
```

####	Development Version

```R
install.packages(c("MALDIquant","MALDIquantForeign"),type="source")
devtools::install_github("unreno/MGMS2")
library("MGMS2")
```
or
```BASH
R -e 'devtools::install_github("unreno/MGMS2")'
```




###	Usage

```R
library('MGMS2')
```


###	Removal

```R
remove.packages("MGMS2")
```








##	Developer

###	Initial Creation

```R
install.packages(c('devtools','roxygen2'))
library('devtools')
library('roxygen2')
create('MGMS2')
```

The "create()" function fails with "Error: Directory 'MGMS2' does not exist." after 
it creates `DESCRIPTION` and `NAMESPACE` files and an empty `R/` dir. Kinda pointless.
Could've just used a simple template for this.
Is that what its supposed to do?



###	Creation / Updation / Posting

Create/Edit `.R` files in `R/` folder with properly formated documentation above functions.

Then update `NAMESPACE` and `man/` from within R ...

```R
library('devtools')
document()
```

Or from the bash command line ...
```BASH
R -e 'library(devtools);document()'
```


Install from the bash command line ...
```BASH
R -e 'library(devtools);document();setwd("..");install("MGMS2",upgrade=FALSE)'
```






Or from R with some manual checking with sample data ...
```R
remove.packages("MGMS2")
library('devtools')
document()

setwd('..')
install('MGMS2',upgrade=FALSE)
library('MGMS2')
lsf.str("package:MGMS2")
"MGMS2" %in% rownames(installed.packages())



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

gather_summary_file(directory=tempdir())

mixture.ratio <- list()
mixture.ratio['A']=1
mixture.ratio['B']=0.5
mixture.ratio['C']=0

sim.template <- create_insilico_mixture_template(mono.info)

insilico.spectrum <- simulate_poly_spectra(sim.template, mixture.ratio)

insilico.spectrums <- simulate_many_poly_spectra(mono.info, mixture.ratio=mixture.ratio, nsim=10 )







remove.packages("MGMS2")
```





ADD AUTOMATED TESTING

ADD AUTOMATED VERSIONING







Then commit and push to github repo
```BASH
git commit -m "<LIST UPDATES>" -a
git push
```



https://cran.r-project.org/

https://xmpalantir.wu.ac.at/cransubmit/

Check for CRAN issues ...
```BASH
R -e 'library(devtools);document();setwd("..");install("MGMS2",upgrade=FALSE)'
cd ..
R CMD build MGMS2
R CMD check MGMS2_*.tar.gz
R CMD check --as-cran MGMS2_*.tar.gz
```

If all is well, try to submit to CRAN




##	Packaging References

https://hilaryparker.com/2014/04/29/writing-an-r-package-from-scratch/

https://kbroman.org/pkg_primer/

https://aberdeenstudygroup.github.io/studyGroup/lessons/SG-T4-Rpackages/makingRpackages1_Intro/

http://r-pkgs.had.co.nz

https://r-pkgs.org/release.html


