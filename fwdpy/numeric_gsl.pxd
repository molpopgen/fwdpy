#Cython API for numeric routines based on the 
#GSL

from libcpp.utility cimport pair
from libcpp.vector cimport vector 
from cython_gsl cimport gsl_matrix,gsl_vector

#Regress v onto m via QR decomposition.
#Return value is the total sum of squares and then
#the vector of squared sums for each column in m.
#Note: m needs to be "ready to go"--no filtering,
#adding of intercept, etc., is done.
#Futher note: m is modified by the QR decomposition,
#so if you need it again later, copy it before calling this.
cdef pair[double,vector[double]] sum_of_squares(const gsl_vector * v,
                                                gsl_matrix * m) nogil
