.SUFFIXES: .rst

%.nbconvert.ipynb: %.ipynb
	PYTHONPATH=../.. jupyter nbconvert --ExecutePreprocessor.timeout=3600 --to notebook --execute $<

%.md: %.nbconvert.ipynb 
	jupyter nbconvert --to markdown $<
	$(eval NBF:=$(basename $@ .md))
	$(eval XX=$(NBF))
	$(eval XX := ${XX:%=%.nbconvert.md})
	mv -f $(XX) $@

%.rst: %.nbconvert.ipynb 
	jupyter nbconvert --to rst $<
	$(eval NBF:=$(basename $@ .rst))
	$(eval XX=$(NBF))
	$(eval XX := ${XX:%=%.nbconvert.rst})
	mv -f $(XX) $@

RST=BGS.rst copypops.rst FixationTimes1.rst temporalSampling.rst trajectories.rst views.rst windows.rst 
MD=BGS.md copypops.md FixationTimes1.md temporalSampling.md trajectories.md views.md windows.md 

all: $(RST) $(MD)

clean: 
	rm -f *.rst
	rm -f trajectories.nbconvert_files/*

