from libc.stdint cimport uint32_t, uint64_t, int32_t, int64_t

import cython
import numpy as np

cdef extern from "<algorithm>" namespace "std" nogil:
    Iter upper_bound[Iter, T](Iter, Iter, T)
    
ctypedef fused integral:
    int32_t
    int64_t

# highly non-intuitive function!
def permuterepeat(uint64_t[:] rep, int32_t[:] data, integral[:] order):
    """
    Effectively does the following but uses no additional memory:
        tmp = np.repeat(data, rep)
        order[:] = tmp[order]
        rep = np.cumsum(rep)
    """
    cdef size_t i
    for i in range(1, rep.shape[0]):
        rep[i] += rep[i - 1]
    for i in range(order.shape[0]):
        order[i] = data[upper_bound(&rep[0], &rep[0] + rep.shape[0], order[i]) - &rep[0]]

def permute_inplace(cython.floating[:] v, integral[:] order):
    """
    Sets `v` to `v[order]` in-place (almost; highest bit of `order` elements is exploited)
    Although O(N), substantially slower than doing `tmp = v[order]; v[:] = tmp` (which allocates).
    """
    cdef size_t n = order.shape[0]
    cdef int64_t i, j
    cdef cython.floating tmp
    for i in range(n):
        if order[i] < 0:
            continue
        tmp = v[i]
        j = i
        while order[j] != i:
            v[j] = v[order[j]]
            order[j] = ~order[j]
            j = ~order[j]
        v[j] = tmp
        order[j] = ~order[j]
    for i in range(n):
        order[i] = ~order[i]

def permute(values, order):
    """
    Tries to permute using an allocated temporary array; if that fails, uses `permute_inplace`
    """
    try:
        tmp = values[order]
        values[:] = tmp
    except MemoryError:
        permute_inplace(values, order)

# the following code is taken from PyTables library, which is BSD-licensed
#
# Created: May 18, 2006
# Author:  Francesc Alted - faltet@pytables.com

import numpy
cimport numpy as cnp

# Types, constants, functions, classes & other objects from everywhere
from numpy cimport import_array, ndarray, \
    npy_int8, npy_int16, npy_int32, npy_int64, \
    npy_uint8, npy_uint16, npy_uint32, npy_uint64, \
    npy_float32, npy_float64, \
    npy_float, npy_double, npy_longdouble

# These two types are defined in npy_common.h but not in cython's numpy.pxd
ctypedef unsigned char npy_bool
ctypedef npy_uint16 npy_float16

from libc.stdlib cimport malloc, free
from libc.string cimport memcpy, strncmp

import_array()

ctypedef fused number_type:
    npy_float32
    npy_float64

DEF PYA_QS_STACK = 100
DEF SMALL_QUICKSORT = 15

def keysort(ndarray array1, ndarray array2):
    """Sort array1 in-place. array2 is also sorted following the array1 order.

    array1 can be of any type, except complex or string.  array2 may be made of
    elements on any size.

    """
    cdef size_t size = cnp.PyArray_SIZE(array1)
    cdef size_t elsize1 = cnp.PyArray_ITEMSIZE(array1)
    cdef size_t elsize2 = cnp.PyArray_ITEMSIZE(array2)
    cdef int type_num = cnp.PyArray_TYPE(array1)

    # floating types
    if type_num == cnp.NPY_FLOAT32:
        _keysort[npy_float32](<npy_float32*>array1.data, array2.data, elsize2, size)
    elif type_num == cnp.NPY_FLOAT64:
        _keysort[npy_float64](<npy_float64*>array1.data, array2.data, elsize2, size)
    else:
        raise ValueError("Unknown array datatype")


cdef inline void swap_bytes(char *x, char *y, size_t n) nogil:
    if n == 8:
        (<npy_int64*>x)[0], (<npy_int64*>y)[0] = (<npy_int64*>y)[0], (<npy_int64*>x)[0]
    elif n == 4:
        (<npy_int32*>x)[0], (<npy_int32*>y)[0] = (<npy_int32*>y)[0], (<npy_int32*>x)[0]
    elif n == 2:
        (<npy_int16*>x)[0], (<npy_int16*>y)[0] = (<npy_int16*>y)[0], (<npy_int16*>x)[0]
    else:
        for i in range(n):
            x[i], y[i] = y[i], x[i]


cdef inline int less_than(number_type* a, number_type* b) nogil:
    return a[0] < b[0] or (b[0] != b[0] and a[0] == a[0])


@cython.cdivision(True)
cdef void _keysort(number_type* start1, char* start2, size_t elsize2, size_t n) nogil:
    cdef number_type *pl = start1
    cdef number_type *pr = start1 + (n - 1)

    cdef char *ipl = start2
    cdef char *ipr = start2 + (n - 1) * elsize2

    cdef number_type vp
    cdef char *ivp = <char *> malloc(elsize2)

    cdef number_type *stack[PYA_QS_STACK]
    cdef number_type **sptr = stack

    cdef char *istack[PYA_QS_STACK]
    cdef char **isptr = istack

    cdef size_t stack_index = 0

    cdef number_type *pm
    cdef number_type *pi
    cdef number_type *pj
    cdef number_type *pt
    cdef char *ipm
    cdef char *ipi
    cdef char *ipj
    cdef char *ipt

    while True:
        while pr - pl > SMALL_QUICKSORT:
            pm  = pl + ((pr - pl) >> 1)
            ipm  = ipl + ((ipr - ipl)/elsize2 >> 1)*elsize2

            if less_than(pm, pl):
                pm[0], pl[0] =  pl[0], pm[0]
                swap_bytes(ipm, ipl, elsize2)

            if less_than(pr, pm):
                pr[0], pm[0] =  pm[0], pr[0]
                swap_bytes(ipr, ipm, elsize2)

            if less_than(pm, pl):
                pm[0], pl[0] =  pl[0], pm[0]
                swap_bytes(ipm, ipl, elsize2)

            vp = pm[0]

            pi = pl
            ipi = ipl

            pj = pr - 1
            ipj = ipr - elsize2

            pm[0], pj[0] = pj[0], pm[0]
            swap_bytes(ipm, ipj, elsize2)

            while True:
                pi += 1
                ipi += elsize2
                while less_than(pi, &vp):
                    pi += 1
                    ipi += elsize2

                pj -= 1
                ipj -= elsize2
                while less_than(&vp, pj):
                    pj -= 1
                    ipj -= elsize2

                if pi >= pj:
                    break

                pi[0], pj[0] = pj[0], pi[0]
                swap_bytes(ipi, ipj, elsize2)

            pi[0], (pr-1)[0] = (pr-1)[0], pi[0]
            swap_bytes(ipi, ipr-elsize2, elsize2)

            # push largest partition on stack and proceed with the other
            if (pi - pl) < (pr - pi):
                sptr[0] = pi + 1
                sptr[1] = pr
                sptr += 2

                isptr[0] = ipi + elsize2
                isptr[1] = ipr
                isptr += 2

                pr = pi - 1
                ipr = ipi - elsize2
            else:
                sptr[0] = pl
                sptr[1] = pi - 1
                sptr += 2

                isptr[0] = ipl
                isptr[1] = ipi - elsize2
                isptr += 2

                pl = pi + 1
                ipl = ipi + elsize2

        pi = pl + 1
        ipi = ipl + elsize2
        while pi <= pr:
            vp = pi[0]
            memcpy(ivp, ipi, elsize2)

            pj = pi
            pt = pi - 1

            ipj = ipi
            ipt = ipi - elsize2

            while pj > pl and less_than(&vp, pt):
                pj[0] = pt[0]
                pj -= 1
                pt -= 1

                memcpy(ipj, ipt, elsize2)
                ipj -= elsize2
                ipt -= elsize2

            pj[0] = vp
            memcpy(ipj, ivp, elsize2)

            pi += 1
            ipi += elsize2

        if sptr == stack:
            break

        sptr -= 2
        pl = sptr[0]
        pr = sptr[1]

        isptr -= 2
        ipl = isptr[0]
        ipr = isptr[1]

    free(ivp)
