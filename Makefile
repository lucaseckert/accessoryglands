setup:
	Rscript -e "source('R/utils.R'); options('repos'='https://cloud.r-project.org'); install_pkgs()"

all:
	Rscript -e "library(targets); library(future); plan('multicore'); tar_make_future()"

## FIXME: avoid re-pushing?
push: ag_bayesdiag.html ag_supp.html
	rsync -uav $^ ms.mcmaster.ca:~/public_html/AG
#	scp $^ ms.mcmaster.ca:~/public_html/AG/
