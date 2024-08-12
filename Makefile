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

get-bt:
	wget http://www.evolution.reading.ac.uk/BayesTraitsV4.0.0/Files/BayesTraitsV4.0.0-Linux.tar.gz
	wget http://www.evolution.reading.ac.uk/BayesTraitsV4.0.0/Files/BayesTraitsV4.0.0-Linux-Threaded.tar.gz
	tar zxvf BayesTraitsV4.0.0-Linux.tar.gz
	tar zxvf BayesTraitsV4.0.0-Linux-Threaded.tar.gz
	rm BayesTraitsV4.0.0*.tar.gz

