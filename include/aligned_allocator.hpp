//////////////////////////////////////////////////////////////////////////
// This file is part of Tungsten
//
// Tungsten is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 2.1 of the License, or
// (at your option) any later version.
//
// Tungsten is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Tungsten. If not, see <http://www.gnu.org/licenses/>.
//
//////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2008 Microsoft Corporpation
// Author: Microsoft Corporpation
// Contributors: Robert McGibbon

// Allocator for aligned memory, for STL objects. This code is
// copied from the link below, with slight modifications to
// substitute in posix_memalign.
// http://blogs.msdn.com/b/vcblog/archive/2008/08/28/the-mallocator.aspx

#ifndef _ALIGNED_ALLOCATOR_H
#define _ALIGNED_ALLOCATOR_H
#include <stddef.h>  // Required for size_t and ptrdiff_t and NULL
#include <new>       // Required for placement new and std::bad_alloc
#include <stdexcept> // Required for std::length_error

// The following headers contain stuff that Mallocator uses.
#include <stdlib.h>  // For malloc() and free()
#include <iostream>  // For std::cout
#include <ostream>   // For std::endl
#include "config.h"  // use posix_memalign or _aligned_malloc

namespace Tungsten {

template <typename T, std::size_t Alignment>
class aligned_allocator {
public:

    // The following will be the same for virtually all allocators.
    typedef T * pointer;
    typedef const T * const_pointer;
    typedef T& reference;
    typedef const T& const_reference;
    typedef T value_type;
    typedef size_t size_type;
    typedef ptrdiff_t difference_type;

    T * address(T& r) const {
        return &r;
    }

    const T * address(const T& s) const {
        return &s;
    }

    size_t max_size() const {
        // The following has been carefully written to be independent of
        // the definition of size_t and to avoid signed/unsigned warnings.
        return (static_cast<size_t>(0) - static_cast<size_t>(1)) / sizeof(T);
    }


    // The following must be the same for all allocators.
    template <typename U> struct rebind {
      typedef aligned_allocator<U, Alignment> other;
    };

    bool operator!=(const aligned_allocator& other) const {
        return !(*this == other);
    }

    void construct(T * const p, const T& t) const {
        void * const pv = static_cast<void *>(p);

        new (pv) T(t);
    }

    void destroy(T * const p) const; // Defined below.


    // Returns true if and only if storage allocated from *this
    // can be deallocated from other, and vice versa.
    // Always returns true for stateless allocators.
    bool operator==(const aligned_allocator& other) const {
        return true;
    }


    // Default constructor, copy constructor, rebinding constructor, and destructor.
    // Empty for stateless allocators.
    aligned_allocator() { }

    aligned_allocator(const aligned_allocator&) { }

    template <typename U> aligned_allocator(const aligned_allocator<U, Alignment>&) { }

    ~aligned_allocator() { }


    // The following will be different for each allocator.
    T * allocate(const size_t n) const {
        // The return value of allocate(0) is unspecified.
        // aligned_allocator returns NULL in order to avoid depending
        // on malloc(0)'s implementation-defined behavior
        // (the implementation can define malloc(0) to return NULL,
        // in which case the bad_alloc check below would fire).
        // All allocators can return NULL in this case.
        if (n == 0) {
            return NULL;
        }

        // All allocators should contain an integer overflow check.
        // The Standardization Committee recommends that std::length_error
        // be thrown in the case of integer overflow.
        if (n > max_size()) {
            throw std::length_error("aligned_allocator<T>::allocate() - Integer overflow.");
        }

        // aligned_allocator wraps malloc().
        //void * const pv = malloc(n * sizeof(T));
        void* pv = NULL;
	int err = 0;
#ifdef HAVE_POSIX_MEMALIGN
        err = posix_memalign(&pv, Alignment, n*sizeof(T));
#elif HAVE_ALIGNED_MALLOC
	pv = _aligned_malloc(n*sizeof(T), Alignment);
#else
        #error "Unsupported platform. Neither posix_memalign nor _aligned_malloc found"
#endif

        if (err != 0 || pv == NULL) {
	  // Allocators should throw std::bad_alloc in the case of memory allocation failure.
            throw std::bad_alloc();
        }

        return static_cast<T *>(pv);
    }

    void deallocate(T * const p, const size_t n) const {
        // aligned_allocator prints a diagnostic message to demonstrate
        // what it's doing. Real allocators won't do this.
        // aligned_allocator wraps free().
        free(p);
    }


    // The following will be the same for all allocators that ignore hints.
    template <typename U> T * allocate(const size_t n, const U * /* const hint */) const {
        return allocate(n);
    }


    // Allocators are not required to be assignable, so
    // all allocators should have a private unimplemented
    // assignment operator. Note that this will trigger the
    // off-by-default (enabled under /Wall) warning C4626
    // "assignment operator could not be generated because a
    // base class assignment operator is inaccessible" within
    // the STL headers, but that warning is useless.
private:
    aligned_allocator& operator=(const aligned_allocator&);
};

// A compiler bug causes it to believe that p->~T() doesn't reference p.

#ifdef _MSC_VER
    #pragma warning(push)
    #pragma warning(disable: 4100) // unreferenced formal parameter
#endif

// The definition of destroy() must be the same for all allocators.
template <typename T, std::size_t Alignment> void aligned_allocator<T, Alignment>::destroy(T * const p) const {
    p->~T();
}

#ifdef _MSC_VER
 #pragma warning(pop)
#endif

}
#endif
