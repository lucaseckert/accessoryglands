## AG_statssem.html: AG_statssem.rmd corHMM.bib
##	Rscript -e 'rmarkdown::render("AG_statssem.rmd")'

AG_statssem.html: AG_statssem.qmd ../corHMM.bib
	quarto render $<

AG_csesem.html: AG_csesem.rmd corHMM.bib
	Rscript -e 'rmarkdown::render("AG_csesem.rmd")'




