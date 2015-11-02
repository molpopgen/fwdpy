CXX=c++
CXXFLAGS:=-std=c++11

DIR=check_deps

all: $(DIR)/test_tcmalloc.o $(DIR)/test_fwdpp.o $(DIR)/test_gsl.o
	$(CXX) $(CXXFLAGS) -o  $(DIR)/test_tcmalloc  $(DIR)/test_tcmalloc.o $(LDFLAGS) -ltcmalloc
	$(CXX) $(CXXFLAGS) -o  $(DIR)/test_fwdpp  $(DIR)/test_fwdpp.o $(LDFLAGS) 
	$(CXX) $(CXXFLAGS) -o  $(DIR)/test_gsl  $(DIR)/test_gsl.o $(LDFLAGS) -lgsl -lgslcblas
clean:
	rm -f  $(DIR)/test_tcmalloc  $(DIR)/test_tcmalloc.o
	rm -f  $(DIR)/test_fwdpp  $(DIR)/test_fwdpp.o
	rm -f  $(DIR)/test_gsl  $(DIR)/test_gsl.o
