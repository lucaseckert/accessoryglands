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

