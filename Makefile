constrModel1.html: constrModel1.rmd cache/MK_3state_constr1_mcmc.rds
	echo "rmarkdown::render('constrModel1.rmd')" | R --slave

## cache/MK_3state_constr1_mcmc.rds: ...
##	
