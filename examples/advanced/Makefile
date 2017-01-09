.SUFFIXES: .rst

%.rst: %.ipynb 
	$(eval NBF:=$(basename $@ .rst))
	PYTHONPATH=../.. jupyter nbconvert --ExecutePreprocessor.timeout=3600 --to notebook --execute $<
	$(eval XX=$(NBF))
	$(eval XX := ${XX:%=%.nbconvert.ipynb})
	jupyter nbconvert --to rst $(XX)
	$(eval XX=$(NBF))
	$(eval XX := ${XX:%=%.nbconvert.rst})
	mv -f $(XX) $@

all: BGSmp.rst

clean: 
	rm -f *.rst


