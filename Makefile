setup:
	Rscript -e "source('R/pkg_install.R'); options('repos'='https://cloud.r-project.org'); install_pkgs()"

all:
	Rscript -e "library(targets); library(future); plan('multicore'); tar_make_future()"

## FIXME: avoid re-pushing?
push: ag_bayesdiag.html ag_supp.html
	rsync -uav $^ ms.mcmaster.ca:~/public_html/AG
#	scp $^ ms.mcmaster.ca:~/public_html/AG/

## FIXME: should be in targets file!
figures: mk_flowfig.R R/plotmat.R figures.R
	R CMD BATCH --vanilla figures.R
	R CMD BATCH --vanilla mk_flowfig.R

flowfig: flowfig.R R/plotmat.R figures.R
	R CMD BATCH --vanilla mk_flowfig.R

ag_supp.pdf: ag_supp.rmd
	Rscript -e "rmarkdown::render('ag_supp.rmd', output_format = 'pdf_document')"

clean:
	rm -f *.aux *.log *.Rout *.tikz
