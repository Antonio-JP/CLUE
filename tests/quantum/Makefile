ALL=table1.test table2.test
BASH=sh

all: clean $(ALL)
	@echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
	@echo "Finished all examples"
	@echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

install:
	@echo "Installing requirements..."
	@pip install --no-index --find-links ./packages clue
	@mkdir results

clean:
	@echo "Cleaning result files..."
	@find . -name "\[result\]*.csv" -type f -delete
	@echo "Cleaning log files..."
	@find . -name "*.log" -type f -delete

package:
	@echo "Cleaning the packages folder"
	@rm -rf ./packages/*.whl ./packages/*.tar.gz
	@pip cache purge
	@pip download -d packages -r requirements.txt
	@pip download -d packages wheel
	@echo "Creating the 'tar.gz' for CLUE"
	@cd ../../; python3 setup.py sdist; cp ./dist/CLUE-1.5.3+qiskit.tar.gz ./tests/quantum/packages 

zip: clean package
	@echo "Compressing the folder"
	@rm -f TACAS24_AE.zip
	@zip TACAS24_AE -r *.py *.sh circuits/*.qasm packages/*.whl packages/*.gz Makefile README.md LICENSE requirements.txt

pc_clean:
	@echo "Cleaning the Python compiled cache..."
	@rm -rf __pycache__

test: test.test
	@echo "Executed all tests for checking correctness"

early_light: clean early_light.test
	@echo "Executed Early Light Review Script"

%.test:
	$(BASH) $*.sh

.PHONY: clean pc_clean