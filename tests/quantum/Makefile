ALL=$(shell python3 script.py list)
PYTHON=python3
FLAGS=

all: $(ALL)
	@echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
	@echo "Finished all examples"
	@echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

case_studies:
	@$(PYTHON) script.py case -s 5 -s 6 -s 7 -s 8 -s 9 -s 10
	@echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
	@echo "Finished all case studies"
	@echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

compile:
	$(PYTHON) script.py compile

%:
	$(PYTHON) script.py $@ $(FLAGS)