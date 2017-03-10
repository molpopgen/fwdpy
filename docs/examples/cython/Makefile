.SUFFIXES: .rst

%.html: %.ipynb
	PYTHONPATH=../../.. jupyter nbconvert --ExecutePreprocessor.timeout=3600 --execute $<

%.md: %.ipynb %.html
	jupyter nbconvert --to markdown $<

%.rst: %.ipynb %.html
	jupyter nbconvert --to rst $<

RST=WritingExtensions.rst
MD=WritingExtensions.md

	
all: $(RST) $(MD)

clean: 
	rm -f $(RST) $(MD)

