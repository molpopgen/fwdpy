from libcpp.utility cimport pair

cdef extern from "<functional>" namespace "std" nogil:
    cdef cppclass hash[KEY]:
        size_t operator()(const KEY & k)

cdef extern from "<unordered_set>" namespace "std" nogil:
    cdef cppclass unordered_set[KEY,HASH=*,KEYEQUAL=*,ALLOCATOR=*]:
        #Iterators (it has no reverse iterator)
        cppclass iterator:
            KEY & operator*()
            iterator operator++()
            iterator operator--()
            bint operator==(iterator)
            bint operator!=(iterator)

        #need to work out local_iterator...

        cppclass const_iterator(iterator):
            pass

        iterator begin()
        const_iterator const_begin "begin"()
        iterator end()
        const_iterator const_end "end"()
        
        #Capacity:
        bint empty() const
        size_t size() const
        size_t max_size() const

        #Modifiers:
        void clear()
        pair[iterator,bint] insert(const KEY &)
        #pair[inerator,bint] insert(KEY &&)
        iterator insert(const_iterator hint, const KEY &)
        #iterator insert(const_iterator hint, KEY &&)
        void insert[INPUTITERATOR]( INPUTITERATOR first, INPUTITERATOR last)
        #void insert( initializer_list[KEY] )
        #emplace and emplace_hint skipped b/c Cython can't do that yet

        iterator erase(const_iterator pos)
        iterator erase(const_iterator first, const_iterator last)
        size_t erase(const KEY &)
        void swap( unordered_set[KEY,HASH,KEYEQUAL,ALLOCATOR] & )

        
        #Bucket interface
        size_t bucket_count() const
        size_t max_bucket_count() const
        void rehash(size_t)
        void reserve(size_t)
        
        #hash interface
        HASH hash_function() const
        KEYEQUAL key_eq() const
