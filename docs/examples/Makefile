.SUFFIXES: .rst

%.html: %.ipynb
	PYTHONPATH=../.. jupyter nbconvert --ExecutePreprocessor.timeout=3600 --execute $<

%.md: %.ipynb %.html
	jupyter nbconvert --to markdown $<

%.rst: %.ipynb %.html
	jupyter nbconvert --to rst $<

RST=BGS.rst copypops.rst FixationTimes1.rst temporalSampling.rst trajectories.rst views.rst windows.rst 
MD=BGS.md copypops.md FixationTimes1.md temporalSampling.md trajectories.md views.md windows.md 

all: $(RST) $(MD)

clean: 
	rm -f $(RST) $(MD)
	rm -f trajectories.nbconvert_files/*

