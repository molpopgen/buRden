#Association test software from the Thornton Lab

##Obtaining the source

```
git clone https://github.com/molpopgen/buRden
```

##Dependencies

1. [R](http://r-project.org) v3.1.0 or greater
2. [Rcpp](http://cran.r-project.org/web/packages/Rcpp/index.html)  v0.11.1 or greater

The package has been tested with various version of [GCC](http://gcc.gnu.org) (all versions 4.6.4 or greater) and [clang](http://clang.llvm.org/) (version clang-503.0.38) on OS X Mavericks.

##Installation

The most straightforward approach is via R's [devtools](http://cran.r-project.org/web/packages/devtools/index.html) package:

~~~
library(devtools)
install_github("molpopgen/buRden")
~~~

The package may also be installed the "traditional" way:

```
R CMD INSTALL buRden
```

To install into a custom location:

```
R_LIBS=/path/to/special/place R CMD INSTALL buRden
```

Or, using the devtools approach, invoke R:

~~~
R_LIBS=/path/to/special/place R
~~~

And then, from the R shell:

~~~
library(devtools)
install_github("molpopgen/buRden")
~~~

##Documentation

All documentation is provided via R's help system.  You may build the pdf of the package manual using the following command:

```
R CMD Rd2pdf buRden
```

##Tests implemented:
1. [Madsen and Browning](http://www.plosgenetics.org/article/info%3Adoi%2F10.1371%2Fjournal.pgen.1000384) (2009)
2. [C-alpha](http://www.plosgenetics.org/article/info%3Adoi%2F10.1371%2Fjournal.pgen.1001322)
3. [ESM](http://www.plosgenetics.org/article/info%3Adoi%2F10.1371%2Fjournal.pgen.1003258)
4. A chi-squared test of the null model that the proportions of cases and controls carrying rare variants are the same.  This is equivalent to [Li and Leal's](http://www.ncbi.nlm.nih.gov/pubmed/18691683) "collapsing method".
