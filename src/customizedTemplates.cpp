/**
 *  $File: customizedTemplates.cpp $
 *  $LastChangedDate$
 *  $Rev$
 *
 *  This file is part of simuPOP, a forward-time population genetics
 *  simulation environment. Please visit http://simupop.sourceforge.net
 *  for details.
 *
 *  Copyright (C) 2004 - 2009 Bo Peng (bpeng@mdanderson.org)
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program. If not, see <http://www.gnu.org/licenses/>.
 */


/* this file defined carrary and defdict datatypes that are customzied from
   corresponding types defined in arraymodule.c and collections.c in the
   standard python distribution. */

#include "Python.h"
#include "structmember.h"

#ifdef STDC_HEADERS
#  include <stddef.h>
#else                                                                                           /* !STDC_HEADERS */
#  ifndef DONT_HAVE_SYS_TYPES_H
#    include <sys/types.h>                                                                      /* For size_t */
#  endif                                                                                        /* DONT_HAVE_SYS_TYPES_H */
#endif                                                                                          /* !STDC_HEADERS */

#if PY_VERSION_HEX < 0x03000000

#  if PY_VERSION_HEX < 0x02060000
#    define Py_SIZE(obj) (((PyVarObject *)(obj))->ob_size)
#  endif

/// CPPONLY
template <typename T>
struct arrayobject_template
{
	PyObject_VAR_HEAD
	// pointer to the beginning of the genotype
	T ob_iter;
};

/// CPPONLY
template <typename T>
bool is_carrayobject_template(PyObject * op);

/// CPPONLY
template <typename T>
PyObject *
getarrayitem_template(struct arrayobject_template<T> * op, Py_ssize_t i)
{
	register struct arrayobject_template<T> * ap;

	assert(is_carrayobject_template<T>(op));
	ap = (struct arrayobject_template<T> *)op;
	assert(i >= 0 && i < Py_SIZE(ap));
	return PyInt_FromLong(*(ap->ob_iter + i));
}


/// CPPONLY
template <>
PyObject *
getarrayitem_template<GenoIterator>(struct arrayobject_template<GenoIterator> * op, Py_ssize_t i)
{
	register struct arrayobject_template<GenoIterator> * ap;

	assert(is_carrayobject_template<GenoIterator>(op));
	ap = (struct arrayobject_template<GenoIterator> *)op;
	assert(i >= 0 && i < Py_SIZE(ap));
	// for read only purpose, make sure not to really insert a value using the
	// non-constant version of operator *.
	return PyInt_FromLong(DEREF_ALLELE(ap->ob_iter + i));
}


/// CPPONLY
template <typename T>
int
setarrayitem_template(struct arrayobject_template<T> * /* ap */, Py_ssize_t /* i */, PyObject * /* v */)
{
	return -1;
}


/// CPPONLY
template <>
int
setarrayitem_template<GenoIterator>(struct arrayobject_template<GenoIterator> * ap,
                                    Py_ssize_t i, PyObject * v)
{
	// right now, the longest allele is uint16_t, but we need to be careful.
	int x;

	/* PyArg_Parse's 'b' formatter is for an unsigned char, therefore
	     must use the next size up that is signed ('h') and manually do
	     the overflow checking */
	if (!PyArg_Parse(v, "i;array item must be integer", &x))
		return -1;
	// force the value to bool to avoid a warning
#  ifdef BINARYALLELE
	*(ap->ob_iter + i) = (x != 0);
#  else
	REF_ASSIGN_ALLELE(ap->ob_iter + i, TO_ALLELE(x));
#  endif
	return 0;
}


/// CPPONLY
template <>
int
setarrayitem_template<LineageIterator>(struct arrayobject_template<LineageIterator> * ap,
                                       Py_ssize_t i, PyObject * v)
{
	long x;

	/* PyArg_Parse's 'b' formatter is for an unsigned char, therefore
	     must use the next size up that is signed ('h') and manually do
	     the overflow checking */
	if (!PyArg_Parse(v, "l;array item must be long", &x))
		return -1;
	*(ap->ob_iter + i) = x;
	return 0;
}


/// CPPONLY
template <typename T>
PyObject *
carray_new_template(PyTypeObject * /* type */, PyObject * /* args */, PyObject * /* kwds */)
{
	PyErr_SetString(PyExc_TypeError,
		"Can not create carray object from python.");
	return NULL;
}


/// CPPONLY
template <typename T>
PyObject *
carray_init_template(PyTypeObject * /* type */, PyObject * /* args */, PyObject * /* kwds */)
{
	PyErr_SetString(PyExc_TypeError,
		"Can not create carray object from python.");
	return NULL;
}


/// CPPONLY
template <typename T>
PyObject * newcarrayobject_template(T begin, T end);

/// CPPONLY
template <typename T>
void
array_dealloc_template(arrayobject_template<T> * op)
{
	PyObject_Del(op);
}


/// CPPONLY
template <typename T>
PyObject *
array_richcompare_template(PyObject * v, PyObject * w, int op)
{
	// will really has this case?
	if (!is_carrayobject_template<T>(v) && !is_carrayobject_template<T>(w)) {
		Py_INCREF(Py_NotImplemented);
		return Py_NotImplemented;
	}

	// both are array
	if (is_carrayobject_template<T>(v) && is_carrayobject_template<T>(w) ) {
		struct arrayobject_template<T> * va, * wa;
		PyObject * vi = NULL;
		PyObject * wi = NULL;
		int i, k;
		PyObject * res;

		va = (struct arrayobject_template<T> *)v;
		wa = (struct arrayobject_template<T> *)w;

		if (Py_SIZE(va) != Py_SIZE(wa) && (op == Py_EQ || op == Py_NE)) {
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
		for (i = 0; i < Py_SIZE(va) && i < Py_SIZE(wa); i++) {
			vi = getarrayitem_template<T>((struct arrayobject_template<T> *)v, i);
			wi = getarrayitem_template<T>((struct arrayobject_template<T> *)w, i);
			if (vi == NULL || wi == NULL) {
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

		if (k) {
			/* No more items to compare -- compare sizes */
			Py_ssize_t vs = Py_SIZE(va);
			Py_ssize_t ws = Py_SIZE(wa);
			int cmp;
			switch (op) {
			case Py_LT: cmp = vs < ws; break;
			case Py_LE: cmp = vs <= ws; break;
			case Py_EQ: cmp = vs == ws; break;
			case Py_NE: cmp = vs != ws; break;
			case Py_GT: cmp = vs > ws; break;
			case Py_GE: cmp = vs >= ws; break;
			default: return NULL;                                             /* cannot happen */
			}
			if (cmp)
				res = Py_True;
			else
				res = Py_False;
			Py_INCREF(res);
			return res;
		}
		/* We have an item that differs.    First, shortcuts for EQ/NE */
		if (op == Py_EQ) {
			Py_INCREF(Py_False);
			res = Py_False;
		} else if (op == Py_NE) {
			Py_INCREF(Py_True);
			res = Py_True;
		} else {
			/* Compare the final item again using the proper operator */
			res = PyObject_RichCompare(vi, wi, op);
		}
		Py_DECREF(vi);
		Py_DECREF(wi);
		return res;
	} else {
		arrayobject_template<T> * va;
		PyObject * wa, * res;
		bool dir;
		Py_ssize_t vs, ws;                                                                     // direction

		// one of them is not array
		if (is_carrayobject_template<T>(v) ) {
			va = (struct arrayobject_template<T> *)v;
			wa = w;
			dir = true;
		} else {
			va = (struct arrayobject_template<T> *)w;
			wa = v;
			dir = false;
		}

		if (!PySequence_Check(wa) ) {
			// use automatic increase of size?
			PyErr_SetString(PyExc_IndexError, "only sequence can be compared");
			return NULL;
		}

		vs = Py_SIZE(va);
		ws = PySequence_Size(wa);

		if (vs != ws && (op == Py_EQ || op == Py_NE)) {
			/* Shortcut: if the lengths differ, the arrays differ */
			if (op == Py_EQ)
				res = Py_False;
			else
				res = Py_True;
			Py_INCREF(res);
			return res;
		}

		/* Search for the first index where items are different */
		PyObject * vi = NULL;
		PyObject * wi = NULL;
		int k = 1;
		for (int i = 0; i < vs && i < ws; i++) {
			vi = getarrayitem_template<T>(va, i);
			wi = PySequence_GetItem(wa, i);
			if (vi == NULL || wi == NULL) {
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

		if (k) {                                                                              // if equal
			/* No more items to compare -- compare sizes */
			int cmp;
			switch (op) {
			case Py_LT: cmp = vs < ws; break;
			case Py_LE: cmp = vs <= ws; break;
			case Py_EQ: cmp = vs == ws; break;
			case Py_NE: cmp = vs != ws; break;
			case Py_GT: cmp = vs > ws; break;
			case Py_GE: cmp = vs >= ws; break;
			default: return NULL;                                             /* cannot happen */
			}
			if ((cmp && dir) || (!cmp && !dir))
				res = Py_True;
			else
				res = Py_False;
			Py_INCREF(res);
			return res;
		}

		/* We have an item that differs.    First, shortcuts for EQ/NE */
		if (op == Py_EQ) {
			Py_INCREF(Py_False);
			res = Py_False;
		} else if (op == Py_NE) {
			Py_INCREF(Py_True);
			res = Py_True;
		} else {
			/* Compare the final item again using the proper operator */
			int r = PyObject_RichCompareBool(vi, wi, op);
			if ( (r == 0 && dir) || (r != 0 && !dir) ) {           // false
				Py_INCREF(Py_False);
				res = Py_False;
			} else {
				Py_INCREF(Py_True);
				res = Py_True;
			}
		}
		Py_DECREF(vi);
		Py_DECREF(wi);
		return res;
	}
}


/// CPPONLY
template<typename T>
Py_ssize_t array_length_template(struct arrayobject_template<T> * a)
{
	return Py_SIZE(a);
}


/// CPPONLY
template<typename T>
PyObject * array_concat_template(struct arrayobject_template<T> *, PyObject *)
{
	PyErr_SetString(PyExc_TypeError,
		"Can not concat carray object.");
	return NULL;
}


/// CPPONLY
template<typename T>
PyObject * array_repeat_template(struct arrayobject_template<T> *, Py_ssize_t)
{
	PyErr_SetString(PyExc_TypeError,
		"Can not repeat carray object.");
	return NULL;
}


/// CPPONLY
template<typename T>
PyObject * array_item_template(struct arrayobject_template<T> * a, Py_ssize_t i)
{
	if (i < 0 || i >= Py_SIZE(a)) {
		PyErr_SetString(PyExc_IndexError, "array index out of range");
		return NULL;
	}
	return getarrayitem_template<T>(a, i);
}


/// CPPONLY
template<typename T>
PyObject * array_slice_template(struct arrayobject_template<T> * a, Py_ssize_t ilow, Py_ssize_t ihigh)
{
	struct arrayobject_template<T> * np;

	if (ilow < 0)
		ilow = 0;
	else if (ilow > Py_SIZE(a))
		ilow = Py_SIZE(a);
	if (ihigh < 0)
		ihigh = 0;
	if (ihigh < ilow)
		ihigh = ilow;
	else if (ihigh > Py_SIZE(a))
		ihigh = Py_SIZE(a);
	np = (struct arrayobject_template<T> *)newcarrayobject_template<T>(a->ob_iter + ilow, a->ob_iter + ihigh);
	if (np == NULL)
		return NULL;
	return (PyObject *)np;
}


/// CPPONLY
template<typename T>
int array_ass_slice_template(struct arrayobject_template<T> * a, Py_ssize_t ilow, Py_ssize_t ihigh, PyObject * v)
{
	if (v == NULL || a == (struct arrayobject_template<T> *)v) {
		PyErr_BadArgument();
		return -1;
	}

	if (ilow < 0)
		ilow = 0;
	else if (ilow > Py_SIZE(a))
		ilow = Py_SIZE(a);
	if (ihigh < 0)
		ihigh = 0;
	if (ihigh < ilow)
		ihigh = ilow;
	else if (ihigh > Py_SIZE(a))
		ihigh = Py_SIZE(a);

	// use a single number to propagate v
	if (PyNumber_Check(v) ) {
		for (Py_ssize_t i = ilow; i < ihigh; ++i)
			setarrayitem_template<T>(a, i, v);
		return 0;
	}
#  define b ((struct arrayobject_template<T> *)v)
	if (is_carrayobject_template<T>(v)) {                                                  /* v is of array type */
		Py_ssize_t n = Py_SIZE(b);
		if (n != ihigh - ilow) {
			PyErr_SetString(PyExc_ValueError, "Can not extend or thrink slice");
			return -1;
		}
		for (Py_ssize_t i = 0; i < n; ++i)
			setarrayitem_template<T>(a, i + ilow, getarrayitem_template<T>(b, i) );
		return 0;
	}
#  undef b
	/* a general sequence */
	if (PySequence_Check(v) ) {
		Py_ssize_t n = PySequence_Size(v);
		if (n != ihigh - ilow) {
			PyErr_SetString(PyExc_ValueError, "Can not extend or thrink slice");
			return -1;
		}
		// iterator sequence
		for (Py_ssize_t i = 0; i < n; ++i) {
			PyObject * item = PySequence_GetItem(v, i);
			setarrayitem_template<T>(a, i + ilow, item);
			Py_DECREF(item);
		}
		return 0;
	}
	PyErr_SetString(PyExc_ValueError, "Only number or list can be assigned");
	return -1;
}


/// CPPONLY
template <typename T>
Py_ssize_t array_ass_item_template(struct arrayobject_template<T> * a, Py_ssize_t i, PyObject * v)
{
	if (i < 0 || i >= Py_SIZE(a)) {
		PyErr_SetString(PyExc_IndexError,
			"array assignment index out of range");
		return -1;
	}
	if (v == NULL)
		return array_ass_slice_template<T>(a, i, i + 1, v);
	return setarrayitem_template<T>(a, i, v);
}


/// CPPONLY
template <typename T>
PyObject * array_count_template(struct arrayobject_template<T> * self, PyObject * args)
{
	int count = 0;
	int i;
	PyObject * v;

	if (!PyArg_ParseTuple(args, "O:count", &v))
		return NULL;
	for (i = 0; i < Py_SIZE(self); i++) {
		PyObject * selfi = getarrayitem_template<T>(self, i);
		int cmp = PyObject_RichCompareBool(selfi, v, Py_EQ);
		Py_DECREF(selfi);
		if (cmp > 0)
			count++;
		else if (cmp < 0)
			return NULL;
	}
	return PyInt_FromLong((long)count);
}


/// CPPONLY
template <typename T>
PyObject * array_index_template(struct arrayobject_template<T> * self, PyObject * args)
{
	Py_ssize_t i;
	PyObject * v;
	Py_ssize_t start = 0, stop = Py_SIZE(self);

	if (!PyArg_ParseTuple(args, "O|O&O&:index", &v,
			_PyEval_SliceIndex, &start, _PyEval_SliceIndex, &stop))
		return NULL;
	if (start < 0) {
		start += Py_SIZE(self);
		if (start < 0)
			start = 0;
	}

	for (i = start; i < stop; i++) {
		PyObject * selfi = getarrayitem_template<T>(self, i);
		int cmp = PyObject_RichCompareBool(selfi, v, Py_EQ);
		Py_DECREF(selfi);
		if (cmp > 0) {
			return PyInt_FromLong((long)i);
		} else if (cmp < 0)
			return NULL;
	}
	PyErr_SetString(PyExc_ValueError, "array.index(x): x not in list");
	return NULL;
}


/// CPPONLY
template<typename T>
PyObject * array_tolist_template(struct arrayobject_template<T> * self, PyObject * args)
{
	PyObject * list = PyList_New(Py_SIZE(self));
	int i;

	if (!PyArg_ParseTuple(args, ":tolist"))
		return NULL;
	if (list == NULL)
		return NULL;
	for (i = 0; i < Py_SIZE(self); i++) {
		PyObject * v = getarrayitem_template<T>(self, i);
		if (v == NULL) {
			Py_DECREF(list);
			return NULL;
		}
		PyList_SetItem(list, i, v);
	}
	return list;
}


/// CPPONLY
template<typename T>
int array_print_template(struct arrayobject_template<T> * a, FILE * fp, int /* flags */)
{
	int ok = 0;
	Py_ssize_t i, len;
	PyObject * v;

	len = Py_SIZE(a);
	if (len == 0) {
		fprintf(fp, "[]");
		return ok;
	}
	fprintf(fp, "[");
	for (i = 0; i < len && ok == 0; i++) {
		if (i > 0)
			fprintf(fp, ", ");
		v = getarrayitem_template<T>(a, i);
		ok = PyObject_Print(v, fp, 0);
		Py_XDECREF(v);
	}
	fprintf(fp, "]");
	return ok;
}


/// CPPONLY
template<typename T>
PyObject *
array_repr_template(struct arrayobject_template<T> * a)
{
	char buf[256];
	PyObject * s, * t, * comma, * v;
	Py_ssize_t i, len;

	len = Py_SIZE(a);
	if (len == 0) {
		PyOS_snprintf(buf, sizeof(buf), "[]");
		return PyString_FromString(buf);
	}
	PyOS_snprintf(buf, sizeof(buf), "[");
	s = PyString_FromString(buf);
	comma = PyString_FromString(", ");
	for (i = 0; i < len && !PyErr_Occurred(); i++) {
		if (i > 0)
			PyString_Concat(&s, comma);
		v = getarrayitem_template<T>(a, i);
		t = PyObject_Repr(v);
		Py_XDECREF(v);
		PyString_ConcatAndDel(&s, t);
	}
	Py_XDECREF(comma);
	PyString_ConcatAndDel(&s, PyString_FromString("]"));
	return s;
}


extern "C" {
extern PyTypeObject Arraytype;
extern PyTypeObject LineageArraytype;
}


/// CPPONLY
template<typename T>
bool is_carrayobject_template(PyObject * /* op */)
{
	return false;
}


/// CPPONLY
template<>
bool is_carrayobject_template<GenoIterator>(PyObject * op)
{
	return op->ob_type == &Arraytype;
}


/// CPPONLY
template<>
bool is_carrayobject_template<LineageIterator>(PyObject * op)
{
	return op->ob_type == &LineageArraytype;
}


/// CPPONLY
template<typename T>
PyObject * newcarrayobject_template(T /* begin */, T /* end */)
{
	return(NULL);
}


/// CPPONLY
template<>
PyObject * newcarrayobject_template<GenoIterator>(GenoIterator begin, GenoIterator end)
{
	// create an object and copy data
	struct arrayobject_template<GenoIterator> * op;

	op = PyObject_New(struct arrayobject_template<GenoIterator>, &Arraytype);
	if (op == NULL) {
		PyObject_Del(op);
		return PyErr_NoMemory();
	}
	//
	op->ob_iter = begin;
#if PY_VERSION_HEX >= 0x030b0000
#  ifdef MUTANTALLELE
	Py_SET_SIZE(op, end.index() - begin.index());
#  else
	Py_SET_SIZE(op, end - begin);
#  endif
# else
#  ifdef MUTANTALLELE
	Py_SIZE(op) = end.index() - begin.index();
#  else
	Py_SIZE(op) = end - begin;
#  endif
#endif
	return (PyObject *)op;
}


/// CPPONLY
template<>
PyObject * newcarrayobject_template<LineageIterator>(LineageIterator begin, LineageIterator end)
{
	// create an object and copy data
	struct arrayobject_template<LineageIterator> * op;

	op = PyObject_New(struct arrayobject_template<LineageIterator>, &LineageArraytype);
	if (op == NULL) {
		PyObject_Del(op);
		return PyErr_NoMemory();
	}
	//
	op->ob_iter = begin;
#if PY_VERSION_HEX >= 0x030b0000
	Py_SET_SIZE(op, end - begin);
#else
	Py_SIZE(op) = end - begin;
#endif
	return (PyObject *)op;
}


#else  // for Python 3
/* Array object implementation */

template <typename T>
struct arrayobject_template
{
	PyObject_VAR_HEAD
	// pointer to the beginning of the genotype
	T ob_iter;
};

template <typename T>
bool is_carrayobject_template(PyObject * op);

template <typename T>
PyObject * newcarrayobject_template(T begin, T end);

template <typename T>
PyObject *
getarrayitem_template(PyObject * op, Py_ssize_t i)
{
	register struct arrayobject_template<T> * ap;

	assert(is_carrayobject_template<T>(op));
	ap = (struct arrayobject_template<T> *)op;
	assert(i >= 0 && i < Py_SIZE(ap));
	return PyInt_FromLong(*(ap->ob_iter + i) );
}


/// CPPONLY
template <typename T>
int
setarrayitem_template(struct arrayobject_template<T> * ap, ssize_t i, PyObject * v)
{
	return -1;
}


/// CPPONLY
template <>
int
setarrayitem_template<GenoIterator>(struct arrayobject_template<GenoIterator> * ap, ssize_t i, PyObject * v)
{
	// right now, the longest allele is uint16_t, but we need to be careful.
	int x;

	/* PyArg_Parse's 'b' formatter is for an unsigned char, therefore
	     must use the next size up that is signed ('h') and manually do
	     the overflow checking */
	if (!PyArg_Parse(v, "i;array item must be integer", &x))
		return -1;
	// force the value to bool to avoid a warning
#  ifdef BINARYALLELE
	*(ap->ob_iter + i) = (x != 0);
#  else
	REF_ASSIGN_ALLELE(ap->ob_iter + i, TO_ALLELE(x));
#  endif
	return 0;
}


/// CPPONLY
template <>
int setarrayitem_template<LineageIterator>(struct arrayobject_template<LineageIterator> * ap, ssize_t i, PyObject * v)
{
	// right now, the longest allele is uint16_t, but we need to be careful.
	long x;

	/* PyArg_Parse's 'b' formatter is for an unsigned char, therefore
	     must use the next size up that is signed ('h') and manually do
	     the overflow checking */
	if (!PyArg_Parse(v, "l;array item must be long", &x))
		return -1;
	*(ap->ob_iter + i) = x;
	return 0;
}


/* Methods */

template <typename T>
void
array_dealloc_template(struct arrayobject_template<T> * op)
{
	Py_TYPE(op)->tp_free((PyObject *)op);
}


template <typename T>
PyObject *
array_richcompare_template(PyObject * v, PyObject * w, int op)
{
	// will really has this case?
	if (!is_carrayobject_template<T>(v) && !is_carrayobject_template<T>(w)) {
		Py_INCREF(Py_NotImplemented);
		return Py_NotImplemented;
	}

	// both are array
	if (is_carrayobject_template<T>(v) && is_carrayobject_template<T>(w) ) {
		struct arrayobject_template<T> * va, * wa;
		PyObject * vi = NULL;
		PyObject * wi = NULL;
		int i, k;
		PyObject * res;

		va = (struct arrayobject_template<T> *)v;
		wa = (struct arrayobject_template<T> *)w;

		if (Py_SIZE(va) != Py_SIZE(wa) && (op == Py_EQ || op == Py_NE)) {
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
		for (i = 0; i < Py_SIZE(va) && i < Py_SIZE(wa); i++) {
			vi = getarrayitem_template<T>(v, i);
			wi = getarrayitem_template<T>(w, i);
			if (vi == NULL || wi == NULL) {
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

		if (k) {
			/* No more items to compare -- compare sizes */
			ssize_t vs = Py_SIZE(va);
			ssize_t ws = Py_SIZE(wa);
			int cmp;
			switch (op) {
			case Py_LT: cmp = vs < ws; break;
			case Py_LE: cmp = vs <= ws; break;
			case Py_EQ: cmp = vs == ws; break;
			case Py_NE: cmp = vs != ws; break;
			case Py_GT: cmp = vs > ws; break;
			case Py_GE: cmp = vs >= ws; break;
			default: return NULL;                                             /* cannot happen */
			}
			if (cmp)
				res = Py_True;
			else
				res = Py_False;
			Py_INCREF(res);
			return res;
		}
		/* We have an item that differs.    First, shortcuts for EQ/NE */
		if (op == Py_EQ) {
			Py_INCREF(Py_False);
			res = Py_False;
		} else if (op == Py_NE) {
			Py_INCREF(Py_True);
			res = Py_True;
		} else {
			/* Compare the final item again using the proper operator */
			res = PyObject_RichCompare(vi, wi, op);
		}
		Py_DECREF(vi);
		Py_DECREF(wi);
		return res;
	} else {
		struct arrayobject_template<T> * va;
		PyObject * wa, * res;
		bool dir;
		ssize_t vs, ws;                                                                     // direction

		// one of them is not array
		if (is_carrayobject_template<T>(v) ) {
			va = (struct arrayobject_template<T> *)v;
			wa = w;
			dir = true;
		} else {
			va = (struct arrayobject_template<T> *)w;
			wa = v;
			dir = false;
		}

		if (!PySequence_Check(wa) ) {
			// use automatic increase of size?
			PyErr_SetString(PyExc_IndexError, "only sequence can be compared");
			return NULL;
		}

		vs = Py_SIZE(va);
		ws = PySequence_Size(wa);

		if (vs != ws && (op == Py_EQ || op == Py_NE)) {
			/* Shortcut: if the lengths differ, the arrays differ */
			if (op == Py_EQ)
				res = Py_False;
			else
				res = Py_True;
			Py_INCREF(res);
			return res;
		}

		/* Search for the first index where items are different */
		PyObject * vi = NULL;
		PyObject * wi = NULL;
		int k = 1;
		for (int i = 0; i < vs && i < ws; i++) {
			vi = getarrayitem_template<T>((PyObject *)(va), i);
			wi = PySequence_GetItem(wa, i);
			if (vi == NULL || wi == NULL) {
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

		if (k) {                                                                              // if equal
			/* No more items to compare -- compare sizes */
			int cmp;
			switch (op) {
			case Py_LT: cmp = vs < ws; break;
			case Py_LE: cmp = vs <= ws; break;
			case Py_EQ: cmp = vs == ws; break;
			case Py_NE: cmp = vs != ws; break;
			case Py_GT: cmp = vs > ws; break;
			case Py_GE: cmp = vs >= ws; break;
			default: return NULL;                                             /* cannot happen */
			}
			if ((cmp && dir) || (!cmp && !dir))
				res = Py_True;
			else
				res = Py_False;
			Py_INCREF(res);
			return res;
		}

		/* We have an item that differs.    First, shortcuts for EQ/NE */
		if (op == Py_EQ) {
			Py_INCREF(Py_False);
			res = Py_False;
		} else if (op == Py_NE) {
			Py_INCREF(Py_True);
			res = Py_True;
		} else {
			/* Compare the final item again using the proper operator */
			int r = PyObject_RichCompareBool(vi, wi, op);
			if ( (r == 0 && dir) || (r != 0 && !dir) ) {           // false
				Py_INCREF(Py_False);
				res = Py_False;
			} else {
				Py_INCREF(Py_True);
				res = Py_True;
			}
		}
		Py_DECREF(vi);
		Py_DECREF(wi);
		return res;
	}
}


template <typename T>
Py_ssize_t
array_length_template(struct arrayobject_template<T> * a)
{
	return Py_SIZE(a);
}


template <typename T>
PyObject *
array_item_template(struct arrayobject_template<T> * a, Py_ssize_t i)
{
	if (i < 0 || i >= Py_SIZE(a)) {
		PyErr_SetString(PyExc_IndexError, "array index out of range");
		return NULL;
	}
	return getarrayitem_template<T>((PyObject *)a, i);
}


template <typename T>
PyObject *
array_slice_template(struct arrayobject_template<T> * a, Py_ssize_t ilow, Py_ssize_t ihigh)
{
	struct arrayobject_template<T> * np;

	if (ilow < 0)
		ilow = 0;
	else if (ilow > Py_SIZE(a))
		ilow = Py_SIZE(a);
	if (ihigh < 0)
		ihigh = 0;
	if (ihigh < ilow)
		ihigh = ilow;
	else if (ihigh > Py_SIZE(a))
		ihigh = Py_SIZE(a);
	np = (struct arrayobject_template<T> *)newcarrayobject_template<T>(a->ob_iter + ilow, a->ob_iter + ihigh);
	if (np == NULL)
		return NULL;
	return (PyObject *)np;
}


template <typename T>
int
array_ass_slice_template(struct arrayobject_template<T> * a, Py_ssize_t ilow, Py_ssize_t ihigh, PyObject * v)
{
	if (v == NULL || a == (struct arrayobject_template<T> *)v) {
		PyErr_BadArgument();
		return -1;
	}

	if (ilow < 0)
		ilow = 0;
	else if (ilow > Py_SIZE(a))
		ilow = Py_SIZE(a);
	if (ihigh < 0)
		ihigh = 0;
	if (ihigh < ilow)
		ihigh = ilow;
	else if (ihigh > Py_SIZE(a))
		ihigh = Py_SIZE(a);

	// use a single number to propagate v
	if (PyNumber_Check(v) ) {
		for (ssize_t i = ilow; i < ihigh; ++i)
			setarrayitem_template<T>(a, i, v);
		return 0;
	}
#  define b ((struct arrayobject_template<T> *)v)
	if (is_carrayobject_template<T>(v)) {                                                  /* v is of array type */
		ssize_t n = Py_SIZE(b);
		if (n != ihigh - ilow) {
			PyErr_SetString(PyExc_ValueError, "Can not extend or thrink slice");
			return -1;
		}
		for (int i = 0; i < n; ++i)
			setarrayitem_template<T>(a, i + ilow, getarrayitem_template<T>(v, i) );
		return 0;
	}
#  undef b
	/* a general sequence */
	if (PySequence_Check(v) ) {
		ssize_t n = PySequence_Size(v);
		if (n != ihigh - ilow) {
			PyErr_SetString(PyExc_ValueError, "Can not extend or thrink slice");
			return -1;
		}
		// iterator sequence
		for (int i = 0; i < n; ++i) {
			PyObject * item = PySequence_GetItem(v, i);
			setarrayitem_template<T>(a, i + ilow, item);
			Py_DECREF(item);
		}
		return 0;
	}
	PyErr_SetString(PyExc_ValueError, "Only number or list can be assigned");
	return -1;

}


template <typename T>
int
array_ass_item_template(struct arrayobject_template<T> * a, Py_ssize_t i, PyObject * v)
{
	if (i < 0 || i >= Py_SIZE(a)) {
		PyErr_SetString(PyExc_IndexError,
			"array assignment index out of range");
		return -1;
	}
	if (v == NULL)
		return array_ass_slice_template<T>(a, i, i + 1, v);
	return setarrayitem_template<T>(a, i, v);
}


template <typename T>
PyObject *
array_count_template(struct arrayobject_template<T> * self, PyObject * v)
{
	Py_ssize_t count = 0;
	Py_ssize_t i;

	for (i = 0; i < Py_SIZE(self); i++) {
		PyObject * selfi = getarrayitem_template<T>((PyObject *)self, i);
		int cmp = PyObject_RichCompareBool(selfi, v, Py_EQ);
		Py_DECREF(selfi);
		if (cmp > 0)
			count++;
		else if (cmp < 0)
			return NULL;
	}
	return PyLong_FromSsize_t(count);
}


template <typename T>
PyObject *
array_index_template(struct arrayobject_template<T> * self, PyObject * v)
{
	Py_ssize_t i;

	for (i = 0; i < Py_SIZE(self); i++) {
		PyObject * selfi = getarrayitem_template<T>((PyObject *)self, i);
		int cmp = PyObject_RichCompareBool(selfi, v, Py_EQ);
		Py_DECREF(selfi);
		if (cmp > 0) {
			return PyLong_FromLong((long)i);
		}else if (cmp < 0)
			return NULL;
	}
	PyErr_SetString(PyExc_ValueError, "array.index(x): x not in list");
	return NULL;
}


template <typename T>
PyObject *
array_tolist_template(struct arrayobject_template<T> * self, PyObject * unused)
{
	PyObject * list = PyList_New(Py_SIZE(self));
	Py_ssize_t i;

	if (list == NULL)
		return NULL;
	for (i = 0; i < Py_SIZE(self); i++) {
		PyObject * v = getarrayitem_template<T>((PyObject *)self, i);
		if (v == NULL) {
			Py_DECREF(list);
			return NULL;
		}
		PyList_SetItem(list, i, v);
	}
	return list;
}


template <typename T>
PyObject *
array_repr_template(struct arrayobject_template<T> * a)
{
	PyObject * s, * v = NULL;

	v = array_tolist_template<T>(a, NULL);
	s = PyUnicode_FromFormat("%R", v);
	Py_DECREF(v);
	return s;
}


template <typename T>
PyObject *
array_subscr_template(struct arrayobject_template<T> * self, PyObject * item)
{
	if (PyIndex_Check(item)) {
		Py_ssize_t i = PyNumber_AsSsize_t(item, PyExc_IndexError);
		if (i == -1 && PyErr_Occurred()) {
			return NULL;
		}
		if (i < 0)
			i += Py_SIZE(self);
		return array_item_template<T>(self, i);
	}else if (PySlice_Check(item)) {
		Py_ssize_t start, stop, step, slicelength;
#  if PY_VERSION_HEX >= 0x03020000
		if (PySlice_GetIndicesEx((PyObject *)item, Py_SIZE(self),
				&start, &stop, &step, &slicelength) < 0) {
			return NULL;
		}
#  else
		if (PySlice_GetIndicesEx((PySliceObject *)item, Py_SIZE(self),
				&start, &stop, &step, &slicelength) < 0) {
			return NULL;
		}
#  endif
		if (step > 1) {
			PyErr_SetString(PyExc_TypeError,
				"Slice with step > 1 is not supported for type simuPOP.array");
			return NULL;
		}

		if (slicelength <= 0)
			return newcarrayobject_template<T>(self->ob_iter, self->ob_iter);
		return newcarrayobject_template<T>(self->ob_iter + start,
		                                   self->ob_iter + stop);
	}else {
		PyErr_SetString(PyExc_TypeError,
			"array indices must be integers");
		return NULL;
	}
}


template <typename T>
int
array_ass_subscr_template(struct arrayobject_template<T> * self, PyObject * item, PyObject * value)
{
	Py_ssize_t start, stop, step, slicelength, needed;

	struct arrayobject_template<T> * other = NULL;

	if (PyIndex_Check(item)) {
		Py_ssize_t i = PyNumber_AsSsize_t(item, PyExc_IndexError);

		if (i == -1 && PyErr_Occurred())
			return -1;
		if (i < 0)
			i += Py_SIZE(self);
		if (i < 0 || i >= Py_SIZE(self)) {
			PyErr_SetString(PyExc_IndexError,
				"array assignment index out of range");
			return -1;
		}
		if (value == NULL) {
			/* Fall through to slice assignment */
			start = i;
			stop = i + 1;
			step = 1;
			slicelength = 1;
		}else
			return setarrayitem_template<T>(self, i, value);
	}else if (PySlice_Check(item)) {
#  if PY_VERSION_HEX >= 0x03020000
		if (PySlice_GetIndicesEx((PyObject *)item,
				Py_SIZE(self), &start, &stop,
				&step, &slicelength) < 0) {
			return -1;
		}
#  else
		if (PySlice_GetIndicesEx((PySliceObject *)item,
				Py_SIZE(self), &start, &stop,
				&step, &slicelength) < 0) {
			return -1;
		}
#  endif
	}else {
		PyErr_SetString(PyExc_TypeError,
			"array indices must be integer");
		return -1;
	}
	if (value == NULL) {
		other = NULL;
		needed = 0;
	}else if (is_carrayobject_template<T>(value)) {
		other = (struct arrayobject_template<T> *)value;
		needed = Py_SIZE(other);
		if (self == other) {
			/* Special case "self[i:j] = self" -- copy self first */
			int ret;
			value = array_slice_template<T>(other, 0, needed);
			if (value == NULL)
				return -1;
			ret = array_ass_subscr_template<T>(self, item, value);
			Py_DECREF(value);
			return ret;
		}
	}else if (PyLong_Check(value)) {
		for (Py_ssize_t i = 0; start + i < stop; ++i)
			setarrayitem_template(self, start + i, value);
		return 0;
	} else if (PySequence_Check(value)) {
		needed = PySequence_Size(value);
	}   else{
		PyErr_Format(PyExc_TypeError,
			"can only assign array (not \"%.200s\") to array slice",
			Py_TYPE(value)->tp_name);
		return -1;
	}
	/* for 'a[2:1] = ...', the insertion point is 'start', not 'stop' */
	if ((step > 0 && stop < start) ||
	    (step < 0 && stop > start))
		stop = start;

	if (step != 1) {
		PyErr_SetString(PyExc_BufferError,
			"Slice with step > 1 is not supported for type simuPOP.array.");
		return -1;
	}

	if (slicelength != needed) {
		PyErr_SetString(PyExc_BufferError,
			"Slice size must match.");
		return -1;
	}
	if (needed > 0) {
		// copy sequence
		if (is_carrayobject_template<T>(value))
			std::copy(other->ob_iter, other->ob_iter + stop - start, self->ob_iter + start);
		else {
			for (Py_ssize_t i = 0; start + i < stop; ++i)
				setarrayitem_template<T>(self, start + i, PySequence_GetItem(value, i));
		}
	}
	return 0;
}


template <typename T>
PyObject * array_new_template(PyTypeObject * type, PyObject * args, PyObject * kwds)
{
	return NULL;
}


extern "C" {
extern PyTypeObject Arraytype;
extern PyTypeObject LineageArraytype;
}

template<typename T>
bool is_carrayobject_template(PyObject * op)
{
	return false;
}


template<>
bool is_carrayobject_template<GenoIterator>(PyObject * op)
{
	return PyObject_TypeCheck(op, &Arraytype);
}


template<>
bool is_carrayobject_template<LineageIterator>(PyObject * op)
{
	return PyObject_TypeCheck(op, &LineageArraytype);
}


template <>
int
array_ass_subscr_template(struct arrayobject_template<GenoIterator> * self, PyObject * item, PyObject * value)
{
	Py_ssize_t start, stop, step, slicelength, needed;

	struct arrayobject_template<GenoIterator> * other = NULL;

	if (PyIndex_Check(item)) {
		Py_ssize_t i = PyNumber_AsSsize_t(item, PyExc_IndexError);

		if (i == -1 && PyErr_Occurred())
			return -1;
		if (i < 0)
			i += Py_SIZE(self);
		if (i < 0 || i >= Py_SIZE(self)) {
			PyErr_SetString(PyExc_IndexError,
				"array assignment index out of range");
			return -1;
		}
		if (value == NULL) {
			/* Fall through to slice assignment */
			start = i;
			stop = i + 1;
			step = 1;
			slicelength = 1;
		}else
			return setarrayitem_template<GenoIterator>(self, i, value);
	}else if (PySlice_Check(item)) {
#  if PY_VERSION_HEX >= 0x03020000
		if (PySlice_GetIndicesEx((PyObject *)item,
				Py_SIZE(self), &start, &stop,
				&step, &slicelength) < 0) {
			return -1;
		}
#  else
		if (PySlice_GetIndicesEx((PySliceObject *)item,
				Py_SIZE(self), &start, &stop,
				&step, &slicelength) < 0) {
			return -1;
		}
#  endif
	}else {
		PyErr_SetString(PyExc_TypeError,
			"array indices must be integer");
		return -1;
	}
	if (value == NULL) {
		other = NULL;
		needed = 0;
	}else if (is_carrayobject_template<GenoIterator>(value)) {
		other = (struct arrayobject_template<GenoIterator> *)value;
		needed = Py_SIZE(other);
		if (self == other) {
			/* Special case "self[i:j] = self" -- copy self first */
			int ret;
			value = array_slice_template<GenoIterator>(other, 0, needed);
			if (value == NULL)
				return -1;
			ret = array_ass_subscr_template<GenoIterator>(self, item, value);
			Py_DECREF(value);
			return ret;
		}
	}else if (PyLong_Check(value)) {
		for (Py_ssize_t i = 0; start + i < stop; ++i)
			setarrayitem_template(self, start + i, value);
		return 0;
	} else if (PySequence_Check(value)) {
		needed = PySequence_Size(value);
	}   else{
		PyErr_Format(PyExc_TypeError,
			"can only assign array (not \"%.200s\") to array slice",
			Py_TYPE(value)->tp_name);
		return -1;
	}
	/* for 'a[2:1] = ...', the insertion point is 'start', not 'stop' */
	if ((step > 0 && stop < start) ||
	    (step < 0 && stop > start))
		stop = start;

	if (step != 1) {
		PyErr_SetString(PyExc_BufferError,
			"Slice with step > 1 is not supported for type simuPOP.array.");
		return -1;
	}

	if (slicelength != needed) {
		PyErr_SetString(PyExc_BufferError,
			"Slice size must match.");
		return -1;
	}
	if (needed > 0) {
		// copy sequence
		if (is_carrayobject_template<GenoIterator>(value))
#  ifdef MUTANTALLELE
			simuPOP::copyGenotype(other->ob_iter, other->ob_iter + stop - start, self->ob_iter + start);
#  else
			std::copy(other->ob_iter, other->ob_iter + stop - start, self->ob_iter + start);
#  endif
		else {
			for (Py_ssize_t i = 0; start + i < stop; ++i)
				setarrayitem_template<GenoIterator>(self, start + i, PySequence_GetItem(value, i));
		}
	}
	return 0;
}


template <>
PyObject *
getarrayitem_template<GenoIterator>(PyObject * op, Py_ssize_t i)
{
	register struct arrayobject_template<GenoIterator> * ap;

	assert(is_carrayobject_template<GenoIterator>(op));
	ap = (struct arrayobject_template<GenoIterator> *)op;
	assert(i >= 0 && i < Py_SIZE(ap));
#  ifdef MUTANTALLELE
	// for read only purpose, make sure not to really insert a value using the
	// non-constant version of operator *.
	return PyInt_FromLong((ap->ob_iter + i).value());
#  else
	return PyInt_FromLong(*(ap->ob_iter + i) );
#  endif
}


/// CPPONLY
template <typename T>
PyObject * newcarrayobject_template(T begin, T end)
{
	return NULL;
}


/// CPPONLY
template <>
PyObject * newcarrayobject_template<GenoIterator>(GenoIterator begin, GenoIterator end)
{
	// create an object and copy data
	struct arrayobject_template<GenoIterator> * op;

	op = PyObject_New(struct arrayobject_template<GenoIterator>, &Arraytype);
	if (op == NULL) {
		PyObject_Del(op);
		return PyErr_NoMemory();
	}
	//
	op->ob_iter = begin;
#if PY_VERSION_HEX >= 0x030b0000
	Py_SET_SIZE(op, end - begin);
#else
	Py_SIZE(op) = end - begin;
#endif
	return (PyObject *)op;
}


/// CPPONLY
template <>
PyObject * newcarrayobject_template<LineageIterator>(LineageIterator begin, LineageIterator end)
{
	// create an object and copy data
	struct arrayobject_template<LineageIterator> * op;

	op = PyObject_New(struct arrayobject_template<LineageIterator>, &LineageArraytype);
	if (op == NULL) {
		PyObject_Del(op);
		return PyErr_NoMemory();
	}
	//
	op->ob_iter = begin;
#if PY_VERSION_HEX >= 0x030b0000
	Py_SET_SIZE(op, end - begin);
#else
	Py_SIZE(op) = end - begin;
#endif
	return (PyObject *)op;
}


#endif
