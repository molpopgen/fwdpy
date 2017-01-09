.SUFFIXES: .rst


%.nbconvert.rst: %.ipynb 
	echo $@
	echo $%
	echo $<
	$(eval NBF := `basename $@ .nbconvert.rst`)
	echo "NBF = " ${NBF}
	PYTHONPATH=../.. jupyter nbconvert --ExecutePreprocessor.timeout=3600 --to notebook --execute $<
	jupyter nbconvert --to rst ${NBF}.nbconvert.ipynb
	rm -f ${NBF}

all: BGSmp.nbconvert.rst

clean: rm -f *.nbconvert.rst


