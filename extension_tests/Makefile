all:
	PYTHONPATH=.. CPPFLAGS="-I../include" CXX=$(CXX) CC=$(CC) $(PYTHON) setup.py build_ext -i

clean:
	rm -f *.so *.cpp