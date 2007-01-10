/***************************************************************************
 *     Copyright (C) 2004 by Bo Peng                                                                                 *
 *     bpeng@rice.edu                                                                                                                *
 *                                                                                                                                                 *
 *     $LastChangedDate$
 *     $Rev$
 *
 *     This program is free software; you can redistribute it and/or modify    *
 *     it under the terms of the GNU General Public License as published by    *
 *     the Free Software Foundation; either version 2 of the License, or         *
 *     (at your option) any later version.                                                                     *
 *                                                                                                                                                 *
 *     This program is distributed in the hope that it will be useful,             *
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of                *
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the                 *
 *     GNU General Public License for more details.                                                    *
 *                                                                                                                                                 *
 *     You should have received a copy of the GNU General Public License         *
 *     along with this program; if not, write to the                                                 *
 *     Free Software Foundation, Inc.,                                                                             *
 *     59 Temple Place - Suite 330, Boston, MA    02111-1307, USA.                         *
 ***************************************************************************/

/* Array object implementation */

/* An array is a uniform list -- all items have the same type.
     The item type is restricted to simple C types like int or float */
/* type bit has been added */

/* NOTE: this is modified from arraymodule.c from the standard python
    distribution. */

#include "Python.h"

#ifdef STDC_HEADERS
#include <stddef.h>
#else                                                                                         /* !STDC_HEADERS */
#ifndef DONT_HAVE_SYS_TYPES_H
#include <sys/types.h>                                                        /* For size_t */
#endif                                                                                        /* DONT_HAVE_SYS_TYPES_H */
#endif                                                                                        /* !STDC_HEADERS */

struct arrayobject;                                                             /* Forward */

/* All possible arraydescr values are defined in the vector "descriptors"
 * below.    That's defined later because the appropriate get and set
 * functions aren't visible yet.
 */
struct arraydescr
{
    int typecode;
    int itemsize;
    PyObject * (*getitem)(struct arrayobject *, int);
    int (*setitem)(struct arrayobject *, int, PyObject *);
};

typedef struct arrayobject
{
    PyObject_VAR_HEAD
    // pointer to the beginning of the item.
    struct iterator
    {
        char *ob_item;
        // this will be used by binary type only.
        GenoIterator ob_iter;
#ifdef SIMUMPI
		UINT ob_piece_size;
		UINT ob_piece_begin;
		UINT ob_piece_end;
#endif		
    } ob_iterator;
    // description of the type, the exact get and set item functions.
    struct arraydescr *ob_descr;
} arrayobject;

// redefinition of type...
// staticforward PyTypeObject Arraytype;

bool is_carrayobject(PyObject *op);

// #define is_carrayobject(op) ((op)->ob_type == &Arraytype)

/****************************************************************************
Get and Set functions for each type.
A Get function takes an arrayobject* and an integer index, returning the
array value at that index wrapped in an appropriate PyObject*.
A Set function takes an arrayobject, integer index, and PyObject*; sets
the array value at that index to the raw C data extracted from the PyObject*,
and returns 0 if successful, else nonzero on failure (PyObject* not of an
appropriate type or value).
Note that the basic Get and Set functions do NOT check that the index is
in bounds; that's the responsibility of the caller.
****************************************************************************/

// bit type
static PyObject *
a_getitem(arrayobject *ap, int i)
{
    return PyInt_FromLong( *(ap->ob_iterator.ob_iter+i) );
}


static int
a_setitem(arrayobject *ap, int i, PyObject *v)
{
    short x;
    /* PyArg_Parse's 'b' formatter is for an unsigned char, therefore
         must use the next size up that is signed ('h') and manually do
         the overflow checking */
    if (!PyArg_Parse(v, "h;array item must be integer", &x))
        return -1;
    // force the value to bool to avoid a warning
#ifdef BINARYALLELE
    *(ap->ob_iterator.ob_iter+i) = (x != 0);
#else
    *(ap->ob_iterator.ob_iter+i) = x;
#endif
    return 0;
}


static PyObject *
c_getitem(arrayobject *ap, int i)
{
    return PyString_FromStringAndSize(&((char *)ap->ob_iterator.ob_item)[i], 1);
}


static int
c_setitem(arrayobject *ap, int i, PyObject *v)
{
    char x;
    if (!PyArg_Parse(v, "c;array item must be char", &x))
        return -1;
    if (i >= 0)
        ((char *)ap->ob_iterator.ob_item)[i] = x;
    return 0;
}


static PyObject *
b_getitem(arrayobject *ap, int i)
{
    long x = ((char *)ap->ob_iterator.ob_item)[i];
    if (x >= 128)
        x -= 256;
    return PyInt_FromLong(x);
}


static int
b_setitem(arrayobject *ap, int i, PyObject *v)
{
    short x;
    /* PyArg_Parse's 'b' formatter is for an unsigned char, therefore
         must use the next size up that is signed ('h') and manually do
         the overflow checking */
    if (!PyArg_Parse(v, "h;array item must be integer", &x))
        return -1;
    else if (x < -128)
    {
        PyErr_SetString(PyExc_OverflowError,
            "signed char is less than minimum");
        return -1;
    }
    else if (x > 127)
    {
        PyErr_SetString(PyExc_OverflowError,
            "signed char is greater than maximum");
        return -1;
    }
    if (i >= 0)
        ((char *)ap->ob_iterator.ob_item)[i] = (char)x;
    return 0;
}


static PyObject *
BB_getitem(arrayobject *ap, int i)
{
    long x = ((unsigned char *)ap->ob_iterator.ob_item)[i];
    return PyInt_FromLong(x);
}


static int
BB_setitem(arrayobject *ap, int i, PyObject *v)
{
    unsigned char x;
    /* 'B' == unsigned char, maps to PyArg_Parse's 'b' formatter */
    if (!PyArg_Parse(v, "b;array item must be integer", &x))
        return -1;
    if (i >= 0)
        ((char *)ap->ob_iterator.ob_item)[i] = x;
    return 0;
}


static PyObject *
h_getitem(arrayobject *ap, int i)
{
    return PyInt_FromLong((long) ((short *)ap->ob_iterator.ob_item)[i]);
}


static int
h_setitem(arrayobject *ap, int i, PyObject *v)
{
    short x;
    /* 'h' == signed short, maps to PyArg_Parse's 'h' formatter */
    if (!PyArg_Parse(v, "h;array item must be integer", &x))
        return -1;
    if (i >= 0)
        ((short *)ap->ob_iterator.ob_item)[i] = x;
    return 0;
}


static PyObject *
HH_getitem(arrayobject *ap, int i)
{
    return PyInt_FromLong((long) ((unsigned short *)ap->ob_iterator.ob_item)[i]);
}


static int
HH_setitem(arrayobject *ap, int i, PyObject *v)
{
    int x;
    /* PyArg_Parse's 'h' formatter is for a signed short, therefore
         must use the next size up and manually do the overflow checking */
    if (!PyArg_Parse(v, "i;array item must be integer", &x))
        return -1;
    else if (x < 0)
    {
        PyErr_SetString(PyExc_OverflowError,
            "unsigned short is less than minimum");
        return -1;
    }
    else if (x > USHRT_MAX)
    {
        PyErr_SetString(PyExc_OverflowError,
            "unsigned short is greater than maximum");
        return -1;
    }
    if (i >= 0)
        ((short *)ap->ob_iterator.ob_item)[i] = (short)x;
    return 0;
}


static PyObject *
i_getitem(arrayobject *ap, int i)
{
    return PyInt_FromLong((long) ((int *)ap->ob_iterator.ob_item)[i]);
}


static int
i_setitem(arrayobject *ap, int i, PyObject *v)
{
    int x;
    /* 'i' == signed int, maps to PyArg_Parse's 'i' formatter */
    if (!PyArg_Parse(v, "i;array item must be integer", &x))
        return -1;
    if (i >= 0)
        ((int *)ap->ob_iterator.ob_item)[i] = x;
    return 0;
}


static PyObject *
II_getitem(arrayobject *ap, int i)
{
    return PyLong_FromUnsignedLong(
        (unsigned long) ((unsigned int *)ap->ob_iterator.ob_item)[i]);
}


static int
II_setitem(arrayobject *ap, int i, PyObject *v)
{
    unsigned long x;
    if (PyLong_Check(v))
    {
        x = PyLong_AsUnsignedLong(v);
        if (x == (unsigned long) -1 && PyErr_Occurred())
            return -1;
    }
    else
    {
        long y;
        if (!PyArg_Parse(v, "l;array item must be integer", &y))
            return -1;
        if (y < 0)
        {
            PyErr_SetString(PyExc_OverflowError,
                "unsigned int is less than minimum");
            return -1;
        }
        x = (unsigned long)y;

    }
    if (x > UINT_MAX)
    {
        PyErr_SetString(PyExc_OverflowError,
            "unsigned int is greater than maximum");
        return -1;
    }

    if (i >= 0)
        ((unsigned int *)ap->ob_iterator.ob_item)[i] = (unsigned int)x;
    return 0;
}


static PyObject *
l_getitem(arrayobject *ap, int i)
{
    return PyInt_FromLong(((long *)ap->ob_iterator.ob_item)[i]);
}


static int
l_setitem(arrayobject *ap, int i, PyObject *v)
{
    long x;
    if (!PyArg_Parse(v, "l;array item must be integer", &x))
        return -1;
    if (i >= 0)
        ((long *)ap->ob_iterator.ob_item)[i] = x;
    return 0;
}


static PyObject *
LL_getitem(arrayobject *ap, int i)
{
    return PyLong_FromUnsignedLong(((unsigned long *)ap->ob_iterator.ob_item)[i]);
}


static int
LL_setitem(arrayobject *ap, int i, PyObject *v)
{
    unsigned long x;
    if (PyLong_Check(v))
    {
        x = PyLong_AsUnsignedLong(v);
        if (x == (unsigned long) -1 && PyErr_Occurred())
            return -1;
    }
    else
    {
        long y;
        if (!PyArg_Parse(v, "l;array item must be integer", &y))
            return -1;
        if (y < 0)
        {
            PyErr_SetString(PyExc_OverflowError,
                "unsigned long is less than minimum");
            return -1;
        }
        x = (unsigned long)y;

    }
    if (x > ULONG_MAX)
    {
        PyErr_SetString(PyExc_OverflowError,
            "unsigned long is greater than maximum");
        return -1;
    }

    if (i >= 0)
        ((unsigned long *)ap->ob_iterator.ob_item)[i] = x;
    return 0;
}


static PyObject *
f_getitem(arrayobject *ap, int i)
{
    return PyFloat_FromDouble((double) ((float *)ap->ob_iterator.ob_item)[i]);
}


static int
f_setitem(arrayobject *ap, int i, PyObject *v)
{
    float x;
    if (!PyArg_Parse(v, "f;array item must be float", &x))
        return -1;
    if (i >= 0)
        ((float *)ap->ob_iterator.ob_item)[i] = x;
    return 0;
}


static PyObject *
d_getitem(arrayobject *ap, int i)
{
    return PyFloat_FromDouble(((double *)ap->ob_iterator.ob_item)[i]);
}


static int
d_setitem(arrayobject *ap, int i, PyObject *v)
{
    double x;
    if (!PyArg_Parse(v, "d;array item must be float", &x))
        return -1;
    if (i >= 0)
        ((double *)ap->ob_iterator.ob_item)[i] = x;
    return 0;
}


/* Description of types */
static struct arraydescr descriptors[] =
{
    {'a', 0, a_getitem, a_setitem},
    {'c', sizeof(char), c_getitem, c_setitem},
    {'b', sizeof(char), b_getitem, b_setitem},
    {'B', sizeof(char), BB_getitem, BB_setitem},
    {'h', sizeof(short), h_getitem, h_setitem},
    {'H', sizeof(short), HH_getitem, HH_setitem},
    {'i', sizeof(int), i_getitem, i_setitem},
    {'I', sizeof(int), II_getitem, II_setitem},
    {'l', sizeof(long), l_getitem, l_setitem},
    {'L', sizeof(long), LL_getitem, LL_setitem},
    {'f', sizeof(float), f_getitem, f_setitem},
    {'d', sizeof(double), d_getitem, d_setitem},
    {                                                                                             /* Sentinel */
        '\0', 0, 0, 0
    }
};

/****************************************************************************
Implementations of array object methods.
****************************************************************************/
// you can not create a object from python,
// error will occur
static PyObject *
carray_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    PyErr_SetString(PyExc_TypeError,
        "Can not create carray object from python.");
    return NULL;
}


// you can not init a object from python,
// error will occur
static PyObject *
carray_init(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    PyErr_SetString(PyExc_TypeError,
        "Can not create carray object from python.");
    return NULL;
}


// declaration only to avoid use of Arraytype
PyObject * newcarrayobject(char* ptr, char type, int size);
#ifdef SIMUMPI
PyObject * newcarrayiterobject(GenoIterator begin, GenoIterator end, 
	UINT s_size, UINT s_begin, UINT s_end);
#else
PyObject * newcarrayiterobject(GenoIterator begin, GenoIterator end);
#endif

static PyObject * getarrayitem(PyObject *op, int i)
{
    register arrayobject *ap;
    assert(is_carrayobject(op));
    ap = (arrayobject *)op;
    if (i < 0 || i >= ap->ob_size)
    {
        // use automatic increase of size?
        PyErr_SetString(PyExc_IndexError, "array index out of range");
        return NULL;
    }
    return (*ap->ob_descr->getitem)(ap, i);
}


static void
array_dealloc(arrayobject *op)
{
    PyObject_Del(op);
}


static PyObject *
array_richcompare(PyObject *v, PyObject *w, int op)
{
    // will really has this case?
    if (!is_carrayobject(v) && !is_carrayobject(w))
    {
        Py_INCREF(Py_NotImplemented);
        return Py_NotImplemented;
    }

    // both are array
    if( is_carrayobject(v) && is_carrayobject(w) )
    {
        arrayobject *va, *wa;
        PyObject *vi = NULL;
        PyObject *wi = NULL;
        int i, k;
        PyObject *res;

        va = (arrayobject *)v;
        wa = (arrayobject *)w;

        if (va->ob_size != wa->ob_size && (op == Py_EQ || op == Py_NE))
        {
            /* Shortcut: if the lengths differ, the arrays differ */
            if (op == Py_EQ)
                res = Py_False;
            else
                res = Py_True;
            Py_INCREF(res);
            return res;
        }

        /* Search for the first index where items are different */
        k = 1;
        for (i = 0; i < va->ob_size && i < wa->ob_size; i++)
        {
            vi = getarrayitem(v, i);
            wi = getarrayitem(w, i);
            if (vi == NULL || wi == NULL)
            {
                Py_XDECREF(vi);
                Py_XDECREF(wi);
                return NULL;
            }
            k = PyObject_RichCompareBool(vi, wi, Py_EQ);
            if (k == 0)
                break;                                                                        /* Keeping vi and wi alive! */
            Py_DECREF(vi);
            Py_DECREF(wi);
            if (k < 0)
                return NULL;
        }

        if (k)
        {
            /* No more items to compare -- compare sizes */
            int vs = va->ob_size;
            int ws = wa->ob_size;
            int cmp;
            switch (op)
            {
                case Py_LT: cmp = vs <    ws; break;
                case Py_LE: cmp = vs <= ws; break;
                case Py_EQ: cmp = vs == ws; break;
                case Py_NE: cmp = vs != ws; break;
                case Py_GT: cmp = vs >    ws; break;
                case Py_GE: cmp = vs >= ws; break;
                default: return NULL;                                         /* cannot happen */
            }
            if (cmp)
                res = Py_True;
            else
                res = Py_False;
            Py_INCREF(res);
            return res;
        }
        /* We have an item that differs.    First, shortcuts for EQ/NE */
        if (op == Py_EQ)
        {
            Py_INCREF(Py_False);
            res = Py_False;
        }
        else if (op == Py_NE)
        {
            Py_INCREF(Py_True);
            res = Py_True;
        }
        else
        {
            /* Compare the final item again using the proper operator */
            res = PyObject_RichCompare(vi, wi, op);
        }
        Py_DECREF(vi);
        Py_DECREF(wi);
        return res;
    }
    else
    {
        arrayobject *va;
        PyObject* wa, *res;
        bool dir;
        int vs, ws;                                                                     // direction

        // one of them is not array
        if( is_carrayobject(v) )
        {
            va = (arrayobject *)v;
            wa = w;
            dir = true;
        }
        else
        {
            va = (arrayobject *)w;
            wa = v;
            dir = false;
        }

        if( ! PySequence_Check(wa) )
        {
            // use automatic increase of size?
            PyErr_SetString(PyExc_IndexError, "only sequence can be compared");
            return NULL;
        }

        vs = va->ob_size;
        ws = PySequence_Size(wa);

        if (vs != ws && (op == Py_EQ || op == Py_NE))
        {
            /* Shortcut: if the lengths differ, the arrays differ */
            if (op == Py_EQ)
                res = Py_False;
            else
                res = Py_True;
            Py_INCREF(res);
            return res;
        }

        /* Search for the first index where items are different */
        PyObject* vi, *wi;
        int k = 1;
        for (int i = 0; i < vs && i < ws; i++)
        {
            vi = getarrayitem((PyObject*)(va), i);
            wi = PySequence_GetItem(wa, i);
            if (vi == NULL || wi == NULL)
            {
                Py_XDECREF(vi);
                Py_XDECREF(wi);
                return NULL;
            }
            k = PyObject_RichCompareBool(vi, wi, Py_EQ);
            if (k == 0)
                break;                                                                        /* Keeping vi and wi alive! */
            Py_DECREF(vi);
            Py_DECREF(wi);
            // -1 for error
            if (k < 0)
                return NULL;
        }

        if (k)                                                                                // if equal
        {
            /* No more items to compare -- compare sizes */
            int cmp;
            switch (op)
            {
                case Py_LT: cmp = vs <    ws; break;
                case Py_LE: cmp = vs <= ws; break;
                case Py_EQ: cmp = vs == ws; break;
                case Py_NE: cmp = vs != ws; break;
                case Py_GT: cmp = vs >    ws; break;
                case Py_GE: cmp = vs >= ws; break;
                default: return NULL;                                         /* cannot happen */
            }
            if ((cmp && dir) || (!cmp && !dir))
                res = Py_True;
            else
                res = Py_False;
            Py_INCREF(res);
            return res;
        }

        /* We have an item that differs.    First, shortcuts for EQ/NE */
        if (op == Py_EQ)
        {
            Py_INCREF(Py_False);
            res = Py_False;
        }
        else if (op == Py_NE)
        {
            Py_INCREF(Py_True);
            res = Py_True;
        }
        else
        {
            /* Compare the final item again using the proper operator */
            int r = PyObject_RichCompareBool(vi, wi, op);
            if( (r==0 && dir) || (r!=0 && !dir) )             // false
            {
                Py_INCREF(Py_False);
                res = Py_False;
            }
            else
            {
                Py_INCREF(Py_True);
                res = Py_True;
            }
        }
        Py_DECREF(vi);
        Py_DECREF(wi);
        return res;
    }
}


static int array_length(arrayobject *a)
{
    return a->ob_size;
}


static PyObject * array_concat(arrayobject *a, PyObject *bb)
{
    PyErr_SetString(PyExc_TypeError,
        "Can not concat carray object.");
    return NULL;
}


static PyObject * array_repeat(arrayobject *a, int n)
{
    PyErr_SetString(PyExc_TypeError,
        "Can not repeat carray object.");
    return NULL;
}


static PyObject * array_item(arrayobject *a, int i)
{
    if (i < 0 || i >= a->ob_size)
    {
        PyErr_SetString(PyExc_IndexError, "array index out of range");
        return NULL;
    }
    return getarrayitem((PyObject *)a, i);
}


static PyObject * array_slice(arrayobject *a, int ilow, int ihigh)
{
    arrayobject *np;
    if (ilow < 0)
        ilow = 0;
    else if (ilow > a->ob_size)
        ilow = a->ob_size;
    if (ihigh < 0)
        ihigh = 0;
    if (ihigh < ilow)
        ihigh = ilow;
    else if (ihigh > a->ob_size)
        ihigh = a->ob_size;
    if( a->ob_descr->typecode == 'a')
        np = (arrayobject *) newcarrayiterobject(a->ob_iterator.ob_iter + ilow,
            a->ob_iterator.ob_iter + ihigh);
    else
        np = (arrayobject *) newcarrayobject(a->ob_iterator.ob_item + ilow*a->ob_descr->itemsize,
            a->ob_descr->typecode, ihigh - ilow);
    if (np == NULL)
        return NULL;
    return (PyObject *)np;
}


static int array_ass_slice(arrayobject *a, int ilow, int ihigh, PyObject *v)
{
    if (v == NULL || a==(arrayobject*)v)
    {
        PyErr_BadArgument();
        return -1;
    }

    if (ilow < 0)
        ilow = 0;
    else if (ilow > a->ob_size)
        ilow = a->ob_size;
    if (ihigh < 0)
        ihigh = 0;
    if (ihigh < ilow)
        ihigh = ilow;
    else if (ihigh > a->ob_size)
        ihigh = a->ob_size;

    // use a single number to propagate v
    if( PyNumber_Check(v) )
    {
        for(int i=ilow; i<ihigh; ++i)
            (*a->ob_descr->setitem)(a, i, v);
        return 0;
    }
#define b ((arrayobject *)v)
    if(is_carrayobject(v))                                                    /* v is of array type */
    {
        int n = b->ob_size;
        if (b->ob_descr != a->ob_descr)
        {
            PyErr_BadArgument();
            return -1;
        }
        if( n != ihigh - ilow)
        {
            PyErr_SetString(PyExc_ValueError, "Can not extend or thrink slice");
            return -1;
        }
        if( a->ob_descr->typecode != 'a')
            memcpy(a->ob_iterator.ob_item + ilow * a->ob_descr->itemsize,
                b->ob_iterator.ob_item, (ihigh-ilow) * a->ob_descr->itemsize);
        else
        {
            for(int i=0; i<n; ++i)
                (*a->ob_descr->setitem)(a, i+ilow, (*b->ob_descr->getitem)(b,i) );
        }
        return 0;
    }
#undef b
    /* a general sequence */
    if( PySequence_Check(v) )
    {
        int n = PySequence_Size(v);
        if( n != ihigh - ilow)
        {
            PyErr_SetString(PyExc_ValueError, "Can not extend or thrink slice");
            return -1;
        }
        // iterator sequence
        for(int i=0; i<n; ++i)
        {
            PyObject* item = PySequence_GetItem(v, i);
            (*a->ob_descr->setitem)(a, i+ilow, item);
            Py_DECREF(item);
        }
        return 0;
    }
    PyErr_SetString(PyExc_ValueError, "Only number or list can be assigned");
    return -1;
}


static int array_ass_item(arrayobject *a, int i, PyObject *v)
{
    if (i < 0 || i >= a->ob_size)
    {
        PyErr_SetString(PyExc_IndexError,
            "array assignment index out of range");
        return -1;
    }
    if (v == NULL)
        return array_ass_slice(a, i, i+1, v);
    return (*a->ob_descr->setitem)(a, i, v);
}


/* not used */
/*
static int setarrayitem(PyObject *a, int i, PyObject *v)
{
    assert(is_carrayobject(a));
    return array_ass_item((arrayobject *)a, i, v);
}
*/

static PyObject * array_count(arrayobject *self, PyObject *args)
{
    int count = 0;
    int i;
    PyObject *v;

    if (!PyArg_ParseTuple(args, "O:count", &v))
        return NULL;
    for (i = 0; i < self->ob_size; i++)
    {
        PyObject *selfi = getarrayitem((PyObject *)self, i);
        int cmp = PyObject_RichCompareBool(selfi, v, Py_EQ);
        Py_DECREF(selfi);
        if (cmp > 0)
            count++;
        else if (cmp < 0)
            return NULL;
    }
    return PyInt_FromLong((long)count);
}


static char count_doc [] =
"count(x)\n\
\n\
Return number of occurences of x in the array.";

static PyObject * array_index(arrayobject *self, PyObject *args)
{
    int i;
    PyObject *v;

    if (!PyArg_ParseTuple(args, "O:index", &v))
        return NULL;
    for (i = 0; i < self->ob_size; i++)
    {
        PyObject *selfi = getarrayitem((PyObject *)self, i);
        int cmp = PyObject_RichCompareBool(selfi, v, Py_EQ);
        Py_DECREF(selfi);
        if (cmp > 0)
        {
            return PyInt_FromLong((long)i);
        }
        else if (cmp < 0)
            return NULL;
    }
    PyErr_SetString(PyExc_ValueError, "array.index(x): x not in list");
    return NULL;
}


static char index_doc [] =
"index(x)\n\
\n\
Return index of first occurence of x in the array.";

static PyObject * array_tolist(arrayobject *self, PyObject *args)
{
    PyObject *list = PyList_New(self->ob_size);
    int i;
    if (!PyArg_ParseTuple(args, ":tolist"))
        return NULL;
    if (list == NULL)
        return NULL;
    for (i = 0; i < self->ob_size; i++)
    {
        PyObject *v = getarrayitem((PyObject *)self, i);
        if (v == NULL)
        {
            Py_DECREF(list);
            return NULL;
        }
        PyList_SetItem(list, i, v);
    }
    return list;
}


static char tolist_doc [] =
"tolist() -> list\n\
\n\
Convert array to an ordinary list with the same items.";

PyMethodDef array_methods[] =
{
    {
        "count", (PyCFunction)array_count, METH_VARARGS,
        count_doc
    },
    {
        "index", (PyCFunction)array_index, METH_VARARGS,
        index_doc
    },
    {
        "tolist",    (PyCFunction)array_tolist,    METH_VARARGS,
        tolist_doc
    },
    {                                                                                             /* sentinel */
        NULL,        NULL
    }
};

static PyObject * array_getattr(arrayobject *a, char *name)
{
    if (strcmp(name, "typecode") == 0)
    {
        char tc = a->ob_descr->typecode;
        return PyString_FromStringAndSize(&tc, 1);
    }
    if (strcmp(name, "itemsize") == 0)
    {
        return PyInt_FromLong((long)a->ob_descr->itemsize);
    }
    if (strcmp(name, "__members__") == 0)
    {
        PyObject *list = PyList_New(2);
        if (list)
        {
            PyList_SetItem(list, 0,
                PyString_FromString("typecode"));
            PyList_SetItem(list, 1,
                PyString_FromString("itemsize"));
            if (PyErr_Occurred())
            {
                Py_DECREF(list);
                list = NULL;
            }
        }
        return list;
    }
    return Py_FindMethod(array_methods, (PyObject *)a, name);
}


static int array_print(arrayobject *a, FILE *fp, int flags)
{
    int ok = 0;
    int i, len;
    PyObject *v;
    len = a->ob_size;
    if (len == 0)
    {
        fprintf(fp, "[]");
        return ok;
    }
    fprintf(fp, "[");
    for (i = 0; i < len && ok == 0; i++)
    {
        if (i > 0)
            fprintf(fp, ", ");
        v = (a->ob_descr->getitem)(a, i);
        ok = PyObject_Print(v, fp, 0);
        Py_XDECREF(v);
    }
    fprintf(fp, "]");
    return ok;
}


static PyObject *
array_repr(arrayobject *a)
{
    char buf[256];
    PyObject *s, *t, *comma, *v;
    int i, len;
    len = a->ob_size;
    if (len == 0)
    {
        PyOS_snprintf(buf, sizeof(buf), "[]");
        return PyString_FromString(buf);
    }
    PyOS_snprintf(buf, sizeof(buf), "[");
    s = PyString_FromString(buf);
    comma = PyString_FromString(", ");
    for (i = 0; i < len && !PyErr_Occurred(); i++)
    {
        if (i > 0)
            PyString_Concat(&s, comma);
        v = (a->ob_descr->getitem)(a, i);
        t = PyObject_Repr(v);
        Py_XDECREF(v);
        PyString_ConcatAndDel(&s, t);
    }
    Py_XDECREF(comma);
    PyString_ConcatAndDel(&s, PyString_FromString("]"));
    return s;
}


static PySequenceMethods array_as_sequence =
{
    (inquiry)array_length,                                                    /*sq_length*/
    (binaryfunc)array_concat,                                             /*sq_concat*/
    (intargfunc)array_repeat,                                             /*sq_repeat*/
    (intargfunc)array_item,                                                 /*sq_item*/
    (intintargfunc)array_slice,                                         /*sq_slice*/
    (intobjargproc)array_ass_item,                                    /*sq_ass_item*/
    (intintobjargproc)array_ass_slice,                            /*sq_ass_slice*/
};

static char arraytype_doc [] =
"An array represents underlying memory of simuPOP structure \n\
so that you can edit the values in python. The type will behave \n\
very much like lists, except that you can change its size.\n\
\n\
Methods:\n\
\n\
count() -- return number of occurences of an object\n\
index() -- return index of first occurence of an object\n\
tolist() -- return the array converted to an ordinary list\n\
\n\
Variables:\n\
\n\
typecode -- the typecode character used to create the array\n\
itemsize -- the length in bytes of one array item\n\
        ";

PyTypeObject Arraytype =
{
    PyObject_HEAD_INIT(NULL)
    0,
    "simuPOP.carray",                                                             /* mudoule.type name */
    sizeof(arrayobject),
    0,
    (destructor)array_dealloc,                                            /* tp_dealloc */
    (printfunc)array_print,                                                 /* tp_print */
    (getattrfunc)array_getattr,                                         /* tp_getattr */
    0,                                                                                            /* tp_setattr */
    0,                                                                                            /* tp_compare */
    (reprfunc)array_repr,                                                     /* tp_repr */
    0,                                                                                            /* tp_as _number*/
    &array_as_sequence,                                                         /* tp_as _sequence*/
    0,                                                                                            /* tp_as _mapping*/
    0,                                                                                            /* tp_hash */
    0,                                                                                            /* tp_call */
    0,                                                                                            /* tp_str */
    0,                                                                                            /* tp_getattro */
    0,                                                                                            /* tp_setattro */
    0,                                                                                            /* tp_as_buffer*/
    Py_TPFLAGS_DEFAULT,                                                         /* tp_flags */
    arraytype_doc,                                                                    /* tp_doc */
    0,                                                                                            /* tp_traverse */
    0,                                                                                            /* tp_clear */
    array_richcompare,                                                            /* tp_richcompare */
    0,                                                                                            /* tp_weaklistoffset */
    0,                                                                                            /* tp_iter */
    0,                                                                                            /* tp_iternext */
    0,                                                                                            /* tp_methods */
    0,                                                                                            /* tp_members */
    0,                                                                                            /* tp_getset */
    0,                                                                                            /* tp_base */
    0,                                                                                            /* tp_dict */
    0,                                                                                            /* tp_descr_get */
    0,                                                                                            /* tp_descr_set */
    0,                                                                                            /* tp_dictoffset */
    (initproc)carray_init,                                                    /* tp_init */
    0,                                                                                            /* tp_alloc */
    carray_new,                                                                         /* tp_new */
};

// we do not import or export hings,
// carray is defined within simuPOP.
void initcarray(void)
{
    // this will be done in PyType_Ready() is your read this
    // from python reference manual.
    Arraytype.ob_type = &PyType_Type;
}


bool is_carrayobject(PyObject* op)
{
    return op->ob_type == &Arraytype;
}


int carray_length(PyObject*a)
{
    return ((arrayobject*)(a))->ob_size;
}


int carray_itemsize(PyObject*a)
{
    return ((arrayobject*)(a))->ob_descr->itemsize;
}


char carray_type(PyObject* a)
{
    return ((arrayobject*)(a))->ob_descr->typecode;
}


char * carray_data(PyObject*a)
{
    return ((arrayobject*)(a))->ob_iterator.ob_item;
}


PyObject * newcarrayobject(char* ptr, char type, int size)
{
    struct arraydescr * descr;

    if (size < 0)
    {
        PyErr_BadInternalCall();
        return NULL;
    }

    // skip the first one, which is for iterator
    for (descr = descriptors+1; descr->typecode != '\0'; descr++)
    {
        if (descr->typecode == type)
        {
            // create an object and copy data
            arrayobject *op;

            op = PyObject_New(arrayobject, &Arraytype);
            if (op == NULL)
            {
                PyObject_Del(op);
                return PyErr_NoMemory();
            }
            op->ob_size = size;
            op->ob_descr = descr;
            op->ob_iterator.ob_item = ptr;
            return (PyObject *) op;
        }
    }
    PyErr_SetString(PyExc_ValueError,
        "bad typecode (must be c, b, B, h, H, i, I, l, L, f or d)");
    return NULL;
}


#ifdef SIMUMPI
PyObject * newcarrayiterobject(GenoIterator begin, GenoIterator end, 
	UINT s_size, UINT s_begin, UINT s_end)
#else
PyObject * newcarrayiterobject(GenoIterator begin, GenoIterator end)
#endif
{
    // create an object and copy data
    arrayobject *op;

    op = PyObject_New(arrayobject, &Arraytype);
    if (op == NULL)
    {
        PyObject_Del(op);
        return PyErr_NoMemory();
    }
    //
    op->ob_size = end - begin;
#ifdef SIMUMPI	
	op->ob_piece_size = s_size;
	op->ob_piece_begin = s_begin;
	op->ob_piece_end = s_end;
#endif	
    op->ob_descr = descriptors;
    op->ob_iterator.ob_iter = begin;
    return (PyObject *) op;
}
