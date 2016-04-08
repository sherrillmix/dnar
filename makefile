VERSION:=$(shell grep Version: DESCRIPTION|sed 's/Version: //')
NAME:=$(shell grep Package: DESCRIPTION|sed 's/Package: //')
PACKAGEFILE:=$(NAME)_$(VERSION).tar.gz
#README.md
all: ../$(PACKAGEFILE)

.PHONY: all install

install:
	R -e 'devtools::install_github("sherrillmix/dnar")'

localInstall:
	R -e 'devtools::install()'

man: R/*.R
	R -e 'devtools::document()'
	touch man


data: data-raw/makeData.R
	R -e 'source("data-raw/makeData.R")'
	touch data

#README.md: README.Rmd
#	R -e 'knitr::opts_chunk$$set(fig.path="README_files/");knitr::knit("README.Rmd")'

#inst/doc: vignettes/*.Rnw R/*.R
	#make localInstall
	#R -e 'devtools::build_vignettes()'
	#touch inst/doc

#data: data-raw/makeAminoColors.R
	#Rscript data-raw/makeAminoColors.R
	#touch data

#inst/doc 
../$(PACKAGEFILE): man R/*.R DESCRIPTION tests/testthat/*.R data
	sed -i "s/^Date:.*$$/Date: `date +%Y-%m-%d`/" DESCRIPTION
	R -e 'devtools::check();devtools::build()'
