.SUFFIXES: .rst

%.html: %.ipynb
	PYTHONPATH=../../.. jupyter nbconvert --ExecutePreprocessor.timeout=3600 --execute $<

%.md: %.ipynb %.html
	jupyter nbconvert --to markdown $<

%.rst: %.ipynb %.html
	jupyter nbconvert --to rst $<

RST=BGSmp.rst SelectedAlleleAges.rst
MD=BGSmp.md SelectedAlleleAges.md

all: $(RST) $(MD)

clean: 
	rm -f $(RST) $(MD)
	rm -f trajectories.nbconvert_files/*

