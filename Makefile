all: 
	make document

document:	
	docker-compose run --rm document

document-docker: 
	R -e "xfun::pkg_attach2('quarto')"
	R -e "quarto::quarto_render('manuscript.qmd')"
