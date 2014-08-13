# Report for job recommender system MAKEFILE
# Jimmy Lin <jimmylin@utexas.edu> 

#######################################################################
# DO NOT EDIT THIS MAKEFILE UNLESS YOU KNOW WHAT YOU ARE DOING. 
#######################################################################	

SRC_DIR = "$(shell pwd)"

DERIVATIVES = $(shell ls *.aux *.log *.out texput.log *.toc *.lot *.lof *.dvi *.bbl *.blg)

all: generate 
	
generate: main.tex main.bib
	pdflatex main 
	bibtex main
	pdflatex main 
	pdflatex main 

push:  clean
	git add main.tex main.pdf main.bib Makefile fig
	git commit -m "update on tex file"
	git push

clean:
	rm -rf $(DERIVATIVES)
