CXX=c++
CXXFLAGS:=-std=c++11

DIR=check_deps

#all: $(DIR)/test_memlibs.o
all: $(DIR)/test_fwdpp.o $(DIR)/test_gsl.o
#	$(CXX) $(CXXFLAGS) -o  $(DIR)/test_memlibs  $(DIR)/test_memlibs.o $(LDFLAGS) $(TESTMEMLIBS)
	$(CXX) $(CXXFLAGS) -o  $(DIR)/test_fwdpp  $(DIR)/test_fwdpp.o $(LDFLAGS) -lgsl -lgslcblas
	$(CXX) $(CXXFLAGS) -o  $(DIR)/test_gsl  $(DIR)/test_gsl.o $(LDFLAGS) -lgsl -lgslcblas
clean:
#	rm -f  $(DIR)/test_memlibs  $(DIR)/test_memlibs.o
	rm -f  $(DIR)/test_fwdpp  $(DIR)/test_fwdpp.o
	rm -f  $(DIR)/test_gsl  $(DIR)/test_gsl.o
