# Report for job recommender system MAKEFILE
# Jimmy Lin <jimmylin@utexas.edu> 

#######################################################################
# DO NOT EDIT THIS MAKEFILE UNLESS YOU KNOW WHAT YOU ARE DOING. 
#######################################################################	

SRC_DIR = "$(shell pwd)"

DERIVATIVES = $(shell ls *.aux *.log *.out texput.log *.toc *.lot *.lof)

all: clean generate clean
	
generate: main.tex
	pdflatex $(SRC_DIR)/main 
	biblatex $(SRC_DIR)/main
	pdflatex $(SRC_DIR)/main 
	pdflatex $(SRC_DIR)/main 

push:  clean
	git add main.tex main.pdf
	git commit -m "update on tex file"
	git push

clean:
	rm -rf $(DERIVATIVES)
