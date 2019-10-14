# specsim

Library for the Interactive Poly-bacterial In-silico Glycolipid Spectrum Simulator




##	Creation

https://github.com/kbroman/pkg_primer

https://kbroman.org/pkg_primer/pages/proper.html

https://kbroman.org/pkg_primer/pages/minimal.html

https://aberdeenstudygroup.github.io/studyGroup/lessons/SG-T4-Rpackages/makingRpackages1_Intro/

http://r-pkgs.had.co.nz/intro.html


```R
install.packages(c('devtools','roxygen2'))
library('devtools')
library('roxygen2')
create('specsim')
```

The "create()" function fails with "Error: Directory 'specsim' does not exist." after creating NAMESPACE.
It only creates `DESCRIPTION` and `NAMESPACE` files and an empty `R/` dir. Kinda pointless.
Could've just used a simple template for this.



##	Creation / Updation / Posting

Create/Edit `.R` files in `R/` folder with properly formated documentation above functions.


Then update `NAMESPACE` and `man/`
```R
library('devtools')
library('roxygen2')

document()
```

Then commit and push to github repo
```BASH

git commit -m "<LIST UPDATES>" -a

git push

```




##	Installation

```R

devtools::install_github("unreno/specsim")

library('specsim')
```


##	Usage


```R
library('specsim')


specsim::testfunction()

```





##	Removal

```R

remove.packages("specsim")


```


