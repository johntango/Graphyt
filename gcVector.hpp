#ifndef GC_VECTOR
#define GC_VECTOR

#ifdef __CUDACC__
#include "managed_allocator.hpp"
template<class T>
using managed_vector = std::vector<T, managed_allocator<T>>;
#else
template<class T>
using managed_vector = std::vector<T>;
#endif

template<class T>
class gcVector{
    managed_vector<T> managed_data;
    T *data;

    gcVector(){
        getPointer();
    }

    void getPointer(){
        data = managed_data.data();
    }

    CPUGPU T& operator[](size_t i)
    {
        return data[i];
    }

    void push_back(T item){
        managed_data.push();
        getPointer(); // May not be necesary, SAM TEST THIS!
    }
};


#endif
