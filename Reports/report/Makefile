# Report for job recommender system MAKEFILE
# Jimmy Lin <jimmylin@utexas.edu> 

#######################################################################
# DO NOT EDIT THIS MAKEFILE UNLESS YOU KNOW WHAT YOU ARE DOING. 
#######################################################################	

SRC_DIR = "./src"
BIN_DIR = "./bin"
CUR_DIR = "'pdf'"

DERIVATIVES = $(shell ls *.aux *.log *.out texput.log *.toc *.lot *.lof)

all:
	pdflatex $(SRC_DIR)/main 
	biblatex $(SRC_DIR)/main
	pdflatex $(SRC_DIR)/main 
	pdflatex $(SRC_DIR)/main 

clean:
	rm -rf $(DERIVATIVES)
