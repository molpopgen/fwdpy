cdef extern from "<memory>" namespace "std":
    cdef cppclass unique_ptr[T,DELETER]:
        T* get()
