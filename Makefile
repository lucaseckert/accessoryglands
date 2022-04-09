setup:
	Rscript -e "source('R/utils.R'); options('repos'='https://cloud.r-project.org'); install_pkgs()"

all:
	Rscript -e "library(targets); library(future); plan('multicore'); tar_make_future()"
