#Expose zlib to Cython

from libc.stdint cimport int64_t

cdef extern from "zlib.h":
    ctypedef void *gzFile
    ctypedef int64_t z_off_t

    int gzclose(gzFile fp)
    gzFile gzopen(char *path, char *mode)
