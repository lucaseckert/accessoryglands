setup:
	Rscript -e "source('R/utils.R'); options('repos'='https://cloud.r-project.org'); install_pkgs()"

all:
	Rscript -e "library(targets); library(future); plan('multicore'); tar_make_future()"

## FIXME: avoid re-pushing?
push: ag_bayesdiag.html ag_supp.html
	rsync -uav $^ ms.mcmaster.ca:~/public_html/AG
#	scp $^ ms.mcmaster.ca:~/public_html/AG/

## FIXME: should be in targets file!
figures: flowfig.R R/plotmat.R figures.R
	R CMD BATCH --vanilla figures.R
	R CMD BATCH --vanilla flowfig.R

flowfig: flowfig.R R/plotmat.R figures.R
	R CMD BATCH --vanilla flowfig.R

ag_supp.pdf: ag_supp.rmd
	Rscript -e "rmarkdown::render('ag_supp.rmd', output_format = 'pdf_document')"

clean:
	rm -f *.aux *.log *.Rout *.tikz

## http://www.evolution.reading.ac.uk/BayesTraitsV4.1.2/Files/BayesTraitsV4.1.2-Linux.tar.gz

## From http://www.evolution.reading.ac.uk/BayesTraitsV4.1.2/BayesTraitsV4.1.2.html: "The threaded versions are no longer being supplied as binaries but can be built via the source code. New threaded methods are under development."
BTVER=4.1.2
get-bt:
	wget http://www.evolution.reading.ac.uk/BayesTraitsV$(BTVER)/Files/BayesTraitsV$(BTVER)-Linux.tar.gz
	- wget http://www.evolution.reading.ac.uk/BayesTraitsV$(BTVER)/Files/BayesTraitsV$(BTVER)-Linux-Threaded.tar.gz
	tar zxvf BayesTraitsV$(BTVER)-Linux.tar.gz
	- tar zxvf BayesTraitsV$(BTVER)-Linux-Threaded.tar.gz
	rm BayesTraitsV$(BTVER)*.tar.gz
## - means 'ignore errors'
