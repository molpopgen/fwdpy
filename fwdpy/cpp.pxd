cdef extern from "<functional>" namespace "std" nogil:
    cdef cppclass hash[KEY]:
        size_t operator()(const KEY & k)
