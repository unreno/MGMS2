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


Then update `NAMESPACE` and `man/`
Do I NEED "library('roxygen2')"? Works without?
```R
library('devtools')

document()
```

Or from the bash command line ...
```BASH
R -e 'library(devtools);document()'
```




Local, test installing
```R
library('devtools')
document()

setwd('..')
install('MGMS2',upgrade=FALSE)
library('MGMS2')
lsf.str("package:MGMS2")

remove.packages("MGMS2")
```

Or from the bash command line ...
```BASH
R -e 'library(devtools);document();setwd("..");install("MGMS2",upgrade=FALSE)'
```



ADD TESTING



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






##	Packaging References

https://hilaryparker.com/2014/04/29/writing-an-r-package-from-scratch/

https://kbroman.org/pkg_primer/

https://aberdeenstudygroup.github.io/studyGroup/lessons/SG-T4-Rpackages/makingRpackages1_Intro/

http://r-pkgs.had.co.nz

https://r-pkgs.org/release.html


