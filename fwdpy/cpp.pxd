cdef extern from "<functional>" namespace "std" nogil:
    cdef cppclass hash[KEY]:
        size_t operator()(const KEY & k)

cdef extern from "<mutex>" namespace "std" nogil:
    cdef cppclass mutex:
        mutex()
