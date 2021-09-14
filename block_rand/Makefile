SHELL = /bin/sh

.PHONY: all

reslst.rda: findopts.R utilityfns.R estimatorfns.R
	R --file=findopts.R

dd_sims.rda: dd_sims.R run_dd_sims.R \
	utilityfns.R \
	estimatorfns.R \
	reslst.rda
	R --file=run_dd_sims.R

g_all.png: summarize_dd_sims.R dd_sims.rda dd_parms.rda
	R --file=summarize_dd_sims.R

block_rand.html: block_rand.Rmd g_all.png
	Rscript -e "rmarkdown::render('block_rand.Rmd',output_format='all')"


