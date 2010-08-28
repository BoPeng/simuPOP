/**
 *  $File: customizedTypes.c $
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

/// CPPONLY
struct arrayobject;                                                             /* Forward */

#  if PY_VERSION_HEX < 0x02060000
#    define Py_SIZE(obj) (((PyVarObject *)(obj))->ob_size)
#  endif

/// CPPONLY
typedef struct arrayobject
{
	PyObject_VAR_HEAD
	// pointer to the beginning of the genotype
	GenoIterator ob_iter;
} arrayobject;

/// CPPONLY
bool is_carrayobject(PyObject * op);

/// CPPONLY
static PyObject *
getarrayitem(arrayobject * op, int i)
{
	register arrayobject * ap;

	assert(is_carrayobject(op));
	ap = (arrayobject *)op;
	assert(i >= 0 && i < Py_SIZE(ap));
	return PyInt_FromLong(*(ap->ob_iter + i) );
}


/// CPPONLY
static int
setarrayitem(arrayobject * ap, int i, PyObject * v)
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
	*(ap->ob_iter + i) = Allele(x);
#  endif
	return 0;
}


/// CPPONLY
static PyObject *
carray_new(PyTypeObject * type, PyObject * args, PyObject * kwds)
{
	PyErr_SetString(PyExc_TypeError,
		"Can not create carray object from python.");
	return NULL;
}


/// CPPONLY
static PyObject *
carray_init(PyTypeObject * type, PyObject * args, PyObject * kwds)
{
	PyErr_SetString(PyExc_TypeError,
		"Can not create carray object from python.");
	return NULL;
}


/// CPPONLY
PyObject * newcarrayobject(GenoIterator begin, GenoIterator end);

/// CPPONLY
static void
array_dealloc(arrayobject * op)
{
	PyObject_Del(op);
}


/// CPPONLY
static PyObject *
array_richcompare(PyObject * v, PyObject * w, int op)
{
	// will really has this case?
	if (!is_carrayobject(v) && !is_carrayobject(w)) {
		Py_INCREF(Py_NotImplemented);
		return Py_NotImplemented;
	}

	// both are array
	if (is_carrayobject(v) && is_carrayobject(w) ) {
		arrayobject * va, * wa;
		PyObject * vi = NULL;
		PyObject * wi = NULL;
		int i, k;
		PyObject * res;

		va = (arrayobject *)v;
		wa = (arrayobject *)w;

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
			vi = getarrayitem((arrayobject *)v, i);
			wi = getarrayitem((arrayobject *)w, i);
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
			int vs = Py_SIZE(va);
			int ws = Py_SIZE(wa);
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
		arrayobject * va;
		PyObject * wa, * res;
		bool dir;
		int vs, ws;                                                                     // direction

		// one of them is not array
		if (is_carrayobject(v) ) {
			va = (arrayobject *)v;
			wa = w;
			dir = true;
		} else {
			va = (arrayobject *)w;
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
			vi = getarrayitem(va, i);
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
static Py_ssize_t array_length(arrayobject * a)
{
	return Py_SIZE(a);
}


/// CPPONLY
static PyObject * array_concat(arrayobject * a, PyObject * bb)
{
	PyErr_SetString(PyExc_TypeError,
		"Can not concat carray object.");
	return NULL;
}


/// CPPONLY
static PyObject * array_repeat(arrayobject * a, Py_ssize_t n)
{
	PyErr_SetString(PyExc_TypeError,
		"Can not repeat carray object.");
	return NULL;
}


/// CPPONLY
static PyObject * array_item(arrayobject * a, Py_ssize_t i)
{
	if (i < 0 || i >= Py_SIZE(a)) {
		PyErr_SetString(PyExc_IndexError, "array index out of range");
		return NULL;
	}
	return getarrayitem(a, i);
}


/// CPPONLY
static PyObject * array_slice(arrayobject * a, Py_ssize_t ilow, Py_ssize_t ihigh)
{
	arrayobject * np;

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
	np = (arrayobject *)newcarrayobject(a->ob_iter + ilow, a->ob_iter + ihigh);
	if (np == NULL)
		return NULL;
	return (PyObject *)np;
}


/// CPPONLY
static int array_ass_slice(arrayobject * a, Py_ssize_t ilow, Py_ssize_t ihigh, PyObject * v)
{
	if (v == NULL || a == (arrayobject *)v) {
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
		for (int i = ilow; i < ihigh; ++i)
			setarrayitem(a, i, v);
		return 0;
	}
#  define b ((arrayobject *)v)
	if (is_carrayobject(v)) {                                                  /* v is of array type */
		int n = Py_SIZE(b);
		if (n != ihigh - ilow) {
			PyErr_SetString(PyExc_ValueError, "Can not extend or thrink slice");
			return -1;
		}
		for (int i = 0; i < n; ++i)
			setarrayitem(a, i + ilow, getarrayitem(b, i) );
		return 0;
	}
#  undef b
	/* a general sequence */
	if (PySequence_Check(v) ) {
		int n = PySequence_Size(v);
		if (n != ihigh - ilow) {
			PyErr_SetString(PyExc_ValueError, "Can not extend or thrink slice");
			return -1;
		}
		// iterator sequence
		for (int i = 0; i < n; ++i) {
			PyObject * item = PySequence_GetItem(v, i);
			setarrayitem(a, i + ilow, item);
			Py_DECREF(item);
		}
		return 0;
	}
	PyErr_SetString(PyExc_ValueError, "Only number or list can be assigned");
	return -1;
}


/// CPPONLY
static Py_ssize_t array_ass_item(arrayobject * a, Py_ssize_t i, PyObject * v)
{
	if (i < 0 || i >= Py_SIZE(a)) {
		PyErr_SetString(PyExc_IndexError,
			"array assignment index out of range");
		return -1;
	}
	if (v == NULL)
		return array_ass_slice(a, i, i + 1, v);
	return setarrayitem(a, i, v);
}


/// CPPONLY
static PyObject * array_count(arrayobject * self, PyObject * args)
{
	int count = 0;
	int i;
	PyObject * v;

	if (!PyArg_ParseTuple(args, "O:count", &v))
		return NULL;
	for (i = 0; i < Py_SIZE(self); i++) {
		PyObject * selfi = getarrayitem(self, i);
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
static char count_doc [] =
    "count(x)\n\
\n\
Return number of occurences of x in the array."                                          ;

/// CPPONLY
static PyObject * array_index(arrayobject * self, PyObject * args)
{
	int i;
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
		PyObject * selfi = getarrayitem(self, i);
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


static char index_doc [] =
    "index(x, [start, [stop]])\n\
\n\
Return index of first occurence of x in the array.";

/// CPPONLY
static PyObject * array_tolist(arrayobject * self, PyObject * args)
{
	PyObject * list = PyList_New(Py_SIZE(self));
	int i;

	if (!PyArg_ParseTuple(args, ":tolist"))
		return NULL;
	if (list == NULL)
		return NULL;
	for (i = 0; i < Py_SIZE(self); i++) {
		PyObject * v = getarrayitem(self, i);
		if (v == NULL) {
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
Convert array to an ordinary list with the same items."                                                          ;

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

/// CPPONLY
static PyObject * array_getattr(arrayobject * a, char * name)
{
	if (strcmp(name, "__members__") == 0) {
		PyObject * list = PyList_New(0);
		return list;
	}
	return Py_FindMethod(array_methods, (PyObject *)a, name);
}


/// CPPONLY
static int array_print(arrayobject * a, FILE * fp, int flags)
{
	int ok = 0;
	int i, len;
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
		v = getarrayitem(a, i);
		ok = PyObject_Print(v, fp, 0);
		Py_XDECREF(v);
	}
	fprintf(fp, "]");
	return ok;
}


/// CPPONLY
static PyObject *
array_repr(arrayobject * a)
{
	char buf[256];
	PyObject * s, * t, * comma, * v;
	int i, len;

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
		v = getarrayitem(a, i);
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
#  if PY_VERSION_HEX < 0x02050000
	(inquiry)array_length,                                                      /*sq_length*/
	(binaryfunc)array_concat,                                                   /*sq_concat*/
	(intargfunc)array_repeat,                                                   /*sq_repeat*/
	(intargfunc)array_item,                                                     /*sq_item*/
	(intintargfunc)array_slice,                                                 /*sq_slice*/
	(intobjargproc)array_ass_item,                                              /*sq_ass_item*/
	(intintobjargproc)array_ass_slice,                                          /*sq_ass_slice*/
#  else
	(lenfunc)array_length,                                                      /*sq_length*/
	(binaryfunc)array_concat,                                                   /*sq_concat*/
	(ssizeargfunc)array_repeat,                                                 /*sq_repeat*/
	(ssizeargfunc)array_item,                                                   /*sq_item*/
	(ssizessizeargfunc)array_slice,                                             /*sq_slice*/
	(ssizeobjargproc)array_ass_item,                                            /*sq_ass_item*/
	(ssizessizeobjargproc)array_ass_slice,                                      /*sq_ass_slice*/
#  endif
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
        "                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        ;

PyTypeObject Arraytype =
{
	PyObject_HEAD_INIT(NULL)
	0,
	"simuPOP.carray",   /* mudoule.type name */
	sizeof(arrayobject),
	0,
	(destructor)array_dealloc,  /* tp_dealloc */
	(printfunc)array_print,     /* tp_print */
	(getattrfunc)array_getattr, /* tp_getattr */
	0,                          /* tp_setattr */
	0,                          /* tp_compare */
	(reprfunc)array_repr,       /* tp_repr */
	0,                          /* tp_as _number*/
	&array_as_sequence,         /* tp_as _sequence*/
	0,                          /* tp_as _mapping*/
	0,                          /* tp_hash */
	0,                          /* tp_call */
	0,                          /* tp_str */
	0,                          /* tp_getattro */
	0,                          /* tp_setattro */
	0,                          /* tp_as_buffer*/
	Py_TPFLAGS_DEFAULT,         /* tp_flags */
	arraytype_doc,              /* tp_doc */
	0,                          /* tp_traverse */
	0,                          /* tp_clear */
	array_richcompare,          /* tp_richcompare */
	0,                          /* tp_weaklistoffset */
	0,                          /* tp_iter */
	0,                          /* tp_iternext */
	0,                          /* tp_methods */
	0,                          /* tp_members */
	0,                          /* tp_getset */
	0,                          /* tp_base */
	0,                          /* tp_dict */
	0,                          /* tp_descr_get */
	0,                          /* tp_descr_set */
	0,                          /* tp_dictoffset */
	(initproc)carray_init,      /* tp_init */
	0,                          /* tp_alloc */
	carray_new,                 /* tp_new */
};


/// CPPONLY
bool is_carrayobject(PyObject * op)
{
	return op->ob_type == &Arraytype;
}


/// CPPONLY
PyObject * newcarrayobject(GenoIterator begin, GenoIterator end)
{
	// create an object and copy data
	arrayobject * op;

	op = PyObject_New(arrayobject, &Arraytype);
	if (op == NULL) {
		PyObject_Del(op);
		return PyErr_NoMemory();
	}
	//
	op->ob_iter = begin;
	Py_SIZE(op) = end - begin;
	return (PyObject *)op;
}


/* defdict type *********************************************************/

typedef struct
{
	PyDictObject dict;
} defdictobject;

//static PyTypeObject defdict_type; /* Forward */

static PyObject *
dict_subscript(PyDictObject * mp, register PyObject * key)
{
	PyObject * v;
	long hash;
	PyDictEntry * ep;

	assert(mp->ma_table != NULL);
	if (!PyString_CheckExact(key) ||
	    (hash = ((PyStringObject *)key)->ob_shash) == -1) {
		hash = PyObject_Hash(key);
		if (hash == -1)
			return NULL;
	}
	ep = (mp->ma_lookup)(mp, key, hash);
	if (ep == NULL)
		return NULL;
	v = ep->me_value;
	if (v == NULL) {
		if (!PyDict_CheckExact(mp))
			return PyInt_FromLong(0);
		return NULL;
	} else
		Py_INCREF(v);
	return v;
}


PyDoc_STRVAR(getitem__doc__, "x.__getitem__(y) <==> x[y]");

static PyMethodDef defdict_methods[] = {
	{ "__getitem__", (PyCFunction)dict_subscript, METH_O | METH_COEXIST,
	  getitem__doc__ },
	{ NULL }
};


static PyMemberDef defdict_members[] = {
	{ NULL }
};

static int
defdict_print(defdictobject * dd, FILE * fp, int flags)
{
	int sts;

	fprintf(fp, "defdict(");
	sts = PyDict_Type.tp_print((PyObject *)dd, fp, 0);
	fprintf(fp, ")");
	return sts;
}


static int
defdict_init(PyObject * self, PyObject * args, PyObject * kwds)
{
	return PyDict_Type.tp_init(self, args, kwds);
}


PyDoc_STRVAR(defdict_doc,
	"defdict() --> dict with default value 0\n\
\n\
A defdict compares equal to a dict with the same items.\n\
");

/* See comment in xxsubtype.c */
#  define DEFERRED_ADDRESS(ADDR) 0

static void
defdict_dealloc(defdictobject * dd)
{
	PyDict_Type.tp_dealloc((PyObject *)dd);
}


static int
defdict_traverse(PyObject * self, visitproc visit, void * arg)
{
	return PyDict_Type.tp_traverse(self, visit, arg);
}


static int
defdict_tp_clear(defdictobject * dd)
{
	return PyDict_Type.tp_clear((PyObject *)dd);
}


// will change one of them to my...
static PyMappingMethods defdict_as_mapping = *PyDict_Type.tp_as_mapping;

static PyTypeObject defdict_type = {
	PyObject_HEAD_INIT(DEFERRED_ADDRESS(&PyType_Type))
	0,                              /* ob_size */
	"simupop.defdict",              /* tp_name */
	sizeof(defdictobject),          /* tp_basicsize */
	0,                              /* tp_itemsize */
	/* methods */
	(destructor)defdict_dealloc,    /* tp_dealloc */
	(printfunc)defdict_print,       /* tp_print */
	0,                              /* tp_getattr */
	0,                              /* tp_setattr */
	0,                              /* tp_compare */
	0,                              /* tp_repr */
	0,                              /* tp_as_number */
	0,                              /* tp_as_sequence */
	&defdict_as_mapping,            /* tp_as_mapping */
	0,                              /* tp_hash */
	0,                              /* tp_call */
	0,                              /* tp_str */
	PyObject_GenericGetAttr,        /* tp_getattro */
	0,                              /* tp_setattro */
	0,                              /* tp_as_buffer */
	Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC |
	Py_TPFLAGS_HAVE_WEAKREFS,       /* tp_flags */
	defdict_doc,                    /* tp_doc */
	defdict_traverse,               /* tp_traverse */
	(inquiry)defdict_tp_clear,      /* tp_clear */
	0,                              /* tp_richcompare */
	0,                              /* tp_weaklistoffset*/
	0,                              /* tp_iter */
	0,                              /* tp_iternext */
	defdict_methods,                /* tp_methods */
	defdict_members,                /* tp_members */
	0,                              /* tp_getset */
	DEFERRED_ADDRESS(&PyDict_Type), /* tp_base */
	0,                              /* tp_dict */
	0,                              /* tp_descr_get */
	0,                              /* tp_descr_set */
	0,                              /* tp_dictoffset */
	defdict_init,                   /* tp_init */
	PyType_GenericAlloc,            /* tp_alloc */
	0,                              /* tp_new */
	PyObject_GC_Del,                /* tp_free */
};


PyObject * PyDefDict_New()
{
	defdictobject * obj;

	// This should call PyDict_Type.tp_new and create an object
	obj = (defdictobject *)defdict_type.tp_new((PyTypeObject *)(&defdict_type), NULL, NULL);
	if (obj == NULL) {
		PyObject_Del(obj);
		return PyErr_NoMemory();
	}
	// initialize this object (call PyDict_Type.tp_init)
	PyObject * args = PyTuple_New(0);
	PyDict_Type.tp_init((PyObject *)obj, args, NULL);
	Py_DECREF(args);
	return (PyObject *)obj;
}


bool is_defdict(PyTypeObject * type)
{
	return type == &defdict_type;
}


// we do not import or export hings,
// carray is defined within simuPOP.
/// CPPONLY
int initCustomizedTypes(void)
{
	// this will be done in PyType_Ready() is your read this
	// from python reference manual.
	Arraytype.ob_type = &PyType_Type;
	if (PyType_Ready(&Arraytype) < 0)
		return -1;
	//
	defdict_type.ob_type = &PyType_Type;
	defdict_type.tp_base = &PyDict_Type;
	defdict_type.tp_as_mapping->mp_subscript = (binaryfunc)dict_subscript;

	if (PyType_Ready(&defdict_type) < 0)
		return -1;
	//Py_INCREF(&defdict_type);
	return 0;
}


#else  // for Python 3
/* Array object implementation */

struct arrayobject; /* Forward */

typedef struct arrayobject
{
	PyObject_VAR_HEAD
	// pointer to the beginning of the genotype
	GenoIterator ob_iter;
} arrayobject;

bool is_carrayobject(PyObject * op);

PyObject * newcarrayobject(GenoIterator begin, GenoIterator end);

static PyObject *
getarrayitem(PyObject * op, Py_ssize_t i)
{
	register arrayobject * ap;

	assert(is_carrayobject(op));
	ap = (arrayobject *)op;
	assert(i >= 0 && i < Py_SIZE(ap));
	return PyInt_FromLong(*(ap->ob_iter + i) );
}


/// CPPONLY
static int
setarrayitem(arrayobject * ap, int i, PyObject * v)
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
	*(ap->ob_iter + i) = Allele(x);
#  endif
	return 0;
}


/* Methods */

static void
array_dealloc(arrayobject * op)
{
	Py_TYPE(op)->tp_free((PyObject *)op);
}


static PyObject *
array_richcompare(PyObject * v, PyObject * w, int op)
{
	// will really has this case?
	if (!is_carrayobject(v) && !is_carrayobject(w)) {
		Py_INCREF(Py_NotImplemented);
		return Py_NotImplemented;
	}

	// both are array
	if (is_carrayobject(v) && is_carrayobject(w) ) {
		arrayobject * va, * wa;
		PyObject * vi = NULL;
		PyObject * wi = NULL;
		int i, k;
		PyObject * res;

		va = (arrayobject *)v;
		wa = (arrayobject *)w;

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
			vi = getarrayitem(v, i);
			wi = getarrayitem(w, i);
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
			int vs = Py_SIZE(va);
			int ws = Py_SIZE(wa);
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
		arrayobject * va;
		PyObject * wa, * res;
		bool dir;
		int vs, ws;                                                                     // direction

		// one of them is not array
		if (is_carrayobject(v) ) {
			va = (arrayobject *)v;
			wa = w;
			dir = true;
		} else {
			va = (arrayobject *)w;
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
			vi = PyInt_FromLong(*(va->ob_iter + i));
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


static Py_ssize_t
array_length(arrayobject * a)
{
	return Py_SIZE(a);
}


static PyObject *
array_item(arrayobject * a, Py_ssize_t i)
{
	if (i < 0 || i >= Py_SIZE(a)) {
		PyErr_SetString(PyExc_IndexError, "array index out of range");
		return NULL;
	}
	return getarrayitem((PyObject *)a, i);
}


static PyObject *
array_slice(arrayobject * a, Py_ssize_t ilow, Py_ssize_t ihigh)
{
	arrayobject * np;

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
	np = (arrayobject *)newcarrayobject(a->ob_iter + ilow, a->ob_iter + ihigh);
	if (np == NULL)
		return NULL;
	return (PyObject *)np;
}


static int
array_ass_slice(arrayobject * a, Py_ssize_t ilow, Py_ssize_t ihigh, PyObject * v)
{
	if (v == NULL || a == (arrayobject *)v) {
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
		for (int i = ilow; i < ihigh; ++i)
			setarrayitem(a, i, v);
		return 0;
	}
#  define b ((arrayobject *)v)
	if (is_carrayobject(v)) {                                                  /* v is of array type */
		int n = Py_SIZE(b);
		if (n != ihigh - ilow) {
			PyErr_SetString(PyExc_ValueError, "Can not extend or thrink slice");
			return -1;
		}
		for (int i = 0; i < n; ++i)
			setarrayitem(a, i + ilow, getarrayitem(v, i) );
		return 0;
	}
#  undef b
	/* a general sequence */
	if (PySequence_Check(v) ) {
		int n = PySequence_Size(v);
		if (n != ihigh - ilow) {
			PyErr_SetString(PyExc_ValueError, "Can not extend or thrink slice");
			return -1;
		}
		// iterator sequence
		for (int i = 0; i < n; ++i) {
			PyObject * item = PySequence_GetItem(v, i);
			setarrayitem(a, i + ilow, item);
			Py_DECREF(item);
		}
		return 0;
	}
	PyErr_SetString(PyExc_ValueError, "Only number or list can be assigned");
	return -1;

}


static int
array_ass_item(arrayobject * a, Py_ssize_t i, PyObject * v)
{
	if (i < 0 || i >= Py_SIZE(a)) {
		PyErr_SetString(PyExc_IndexError,
			"array assignment index out of range");
		return -1;
	}
	if (v == NULL)
		return array_ass_slice(a, i, i + 1, v);
	return setarrayitem(a, i, v);
}


static PyObject *
array_count(arrayobject * self, PyObject * v)
{
	Py_ssize_t count = 0;
	Py_ssize_t i;

	for (i = 0; i < Py_SIZE(self); i++) {
		PyObject * selfi = getarrayitem((PyObject *)self, i);
		int cmp = PyObject_RichCompareBool(selfi, v, Py_EQ);
		Py_DECREF(selfi);
		if (cmp > 0)
			count++;
		else if (cmp < 0)
			return NULL;
	}
	return PyLong_FromSsize_t(count);
}


PyDoc_STRVAR(count_doc,
	"count(x)\n\
\n\
Return number of occurrences of x in the array."                     );

static PyObject *
array_index(arrayobject * self, PyObject * v)
{
	Py_ssize_t i;

	for (i = 0; i < Py_SIZE(self); i++) {
		PyObject * selfi = getarrayitem((PyObject *)self, i);
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


PyDoc_STRVAR(index_doc,
	"index(x)\n\
\n\
Return index of first occurrence of x in the array."                     );

static PyObject *
array_tolist(arrayobject * self, PyObject * unused)
{
	PyObject * list = PyList_New(Py_SIZE(self));
	Py_ssize_t i;

	if (list == NULL)
		return NULL;
	for (i = 0; i < Py_SIZE(self); i++) {
		PyObject * v = getarrayitem((PyObject *)self, i);
		if (v == NULL) {
			Py_DECREF(list);
			return NULL;
		}
		PyList_SetItem(list, i, v);
	}
	return list;
}


PyDoc_STRVAR(tolist_doc,
	"tolist() -> list\n\
\n\
Convert array to an ordinary list with the same items."                             );


static PyMethodDef array_methods[] = {
	{ "count",	(PyCFunction)array_count,	METH_O,
	  count_doc },
	{ "index",	(PyCFunction)array_index,	METH_O,
	  index_doc },
	{ "tolist", (PyCFunction)array_tolist,	METH_NOARGS,
	  tolist_doc },
	{ NULL,		NULL }      /* sentinel */
};

static PyObject *
array_repr(arrayobject * a)
{
	PyObject * s, * v = NULL;
	v = array_tolist(a, NULL);
	s = PyUnicode_FromFormat("%R", v);
	Py_DECREF(v);
	return s;
}


static PySequenceMethods array_as_sequence = {
	(lenfunc)array_length,                  /*sq_length*/
	0,                                      /*sq_concat*/
	0,                                      /*sq_repeat*/
	(ssizeargfunc)array_item,               /*sq_item*/
	0,                                      /*sq_slice*/
	(ssizeobjargproc)array_ass_item,        /*sq_ass_item*/
	0,					/*sq_ass_slice*/
	0,                                      /*sq_contains*/
	0,                                      /*sq_inplace_concat*/
	0                                       /*sq_inplace_repeat*/
};

PyObject * array_new(PyTypeObject * type, PyObject * args, PyObject * kwds)
{
	return NULL;
}


PyDoc_STRVAR(arraytype_doc,
	" \n\
\n\
Methods:\n\
\n\
append() -- append a new item to the end of the array\n\
buffer_info() -- return information giving the current memory info\n\
byteswap() -- byteswap all the items of the array\n\
count() -- return number of occurrences of an object\n\
extend() -- extend array by appending multiple elements from an iterable\n\
fromfile() -- read items from a file object\n\
fromlist() -- append items from the list\n\
fromstring() -- append items from the string\n\
index() -- return index of first occurrence of an object\n\
insert() -- insert a new item into the array at a provided position\n\
pop() -- remove and return item (default last)\n\
remove() -- remove first occurrence of an object\n\
reverse() -- reverse the order of the items in the array\n\
tofile() -- write all items to a file object\n\
tolist() -- return the array converted to an ordinary list\n\
tostring() -- return the array converted to a string\n\
\n\
Attributes:\n\
\n\
itemsize -- the length in bytes of one array item\n\
"                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        );

static PyTypeObject Arraytype = {
	PyVarObject_HEAD_INIT(NULL, 0)
	"simuPOP.array",
	sizeof(arrayobject),
	0,
	(destructor)array_dealloc,                  /* tp_dealloc */
	0,                                          /* tp_print */
	0,                                          /* tp_getattr */
	0,                                          /* tp_setattr */
	0,                                          /* tp_reserved */
	(reprfunc)array_repr,                       /* tp_repr */
	0,                                          /* tp_as_number*/
	&array_as_sequence,                         /* tp_as_sequence*/
	0,                                          /* tp_as_mapping*/
	0,                                          /* tp_hash */
	0,                                          /* tp_call */
	0,                                          /* tp_str */
	PyObject_GenericGetAttr,                    /* tp_getattro */
	0,                                          /* tp_setattro */
	0,                                          /* tp_as_buffer*/
	Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,   /* tp_flags */
	arraytype_doc,                              /* tp_doc */
	0,                                          /* tp_traverse */
	0,                                          /* tp_clear */
	array_richcompare,                          /* tp_richcompare */
	0,                                          /* tp_weaklistoffset */
	0,                                          /* tp_iter */
	0,                                          /* tp_iternext */
	array_methods,                              /* tp_methods */
	0,                                          /* tp_members */
	0,                                          /* tp_getset */
	0,                                          /* tp_base */
	0,                                          /* tp_dict */
	0,                                          /* tp_descr_get */
	0,                                          /* tp_descr_set */
	0,                                          /* tp_dictoffset */
	0,                                          /* tp_init */
	PyType_GenericAlloc,                        /* tp_alloc */
	array_new,                                  /* tp_new */
	PyObject_Del,                               /* tp_free */
};


bool is_carrayobject(PyObject * op)
{
	return PyObject_TypeCheck(op, &Arraytype);
}


/// CPPONLY
PyObject * newcarrayobject(GenoIterator begin, GenoIterator end)
{
	// create an object and copy data
	arrayobject * op;

	op = PyObject_New(arrayobject, &Arraytype);
	if (op == NULL) {
		PyObject_Del(op);
		return PyErr_NoMemory();
	}
	//
	op->ob_iter = begin;
	Py_SIZE(op) = end - begin;
	return (PyObject *)op;
}


typedef struct
{
	PyDictObject dict;
} defdictobject;

//static PyTypeObject defdict_type; /* Forward */

PyDoc_STRVAR(defdict_missing_doc,
	"__missing__(key) # Called by __getitem__ for missing key; pseudo-code:\n\
  Return 0\n\
"                                                                                             );

static PyObject *
defdict_missing(defdictobject * dd, PyObject * key)
{
	return PyInt_FromLong(0);
}


PyDoc_STRVAR(defdict_copy_doc, "D.copy() -> a shallow copy of D.");

static PyObject *
defdict_copy(defdictobject * dd)
{
	/* This calls the object's class.  That only works for subclasses
	   whose class constructor has the same signature.  Subclasses that
	   define a different constructor signature must override copy().
	 */
	return PyObject_CallFunctionObjArgs((PyObject *)Py_TYPE(dd), Py_None, dd, NULL);
}


static PyObject *
defdict_reduce(defdictobject * dd)
{
	/* __reduce__ must return a 5-tuple as follows:

	   - additional state (here None)
	   - sequence iterator (here None)
	   - dictionary iterator (yielding successive (key, value) pairs

	   This API is used by pickle.py and copy.py.
	 */
	PyObject * items;
	PyObject * iter;
	PyObject * result;

	items = PyObject_CallMethod((PyObject *)dd, "items", "()");
	if (items == NULL)
		return NULL;
	iter = PyObject_GetIter(items);
	if (iter == NULL) {
		Py_DECREF(items);
		return NULL;
	}
	result = PyTuple_Pack(4, Py_TYPE(dd),
		Py_None, Py_None, iter);
	Py_DECREF(iter);
	Py_DECREF(items);
	return result;
}


static PyMethodDef defdict_methods[] = {
	{ "__missing__", (PyCFunction)defdict_missing, METH_O,
	  defdict_missing_doc },
	{ "copy",		 (PyCFunction)defdict_copy,	   METH_NOARGS,
	  defdict_copy_doc },
	{ "__copy__",	 (PyCFunction)defdict_copy,	   METH_NOARGS,
	  defdict_copy_doc },
	{ "__reduce__",	 (PyCFunction)defdict_reduce,  METH_NOARGS,
	  "" },
	{ NULL }
};

static PyMemberDef defdict_members[] = {
	{ NULL }
};

static void
defdict_dealloc(defdictobject * dd)
{
	PyDict_Type.tp_dealloc((PyObject *)dd);
}


static PyObject *
defdict_repr(defdictobject * dd)
{
	PyObject * baserepr;
	PyObject * result;

	baserepr = PyDict_Type.tp_repr((PyObject *)dd);
	if (baserepr == NULL)
		return NULL;
	result = PyUnicode_FromFormat("defdict(%U)",
		baserepr);
	Py_DECREF(baserepr);
	return result;
}


static int
defdict_traverse(PyObject * self, visitproc visit, void * arg)
{
	return PyDict_Type.tp_traverse(self, visit, arg);
}


static int
defdict_tp_clear(defdictobject * dd)
{
	return PyDict_Type.tp_clear((PyObject *)dd);
}


static int
defdict_init(PyObject * self, PyObject * args, PyObject * kwds)
{
	return PyDict_Type.tp_init(self, args, kwds);
}


PyDoc_STRVAR(defdict_doc,
	"defdict() --> dict with default value\n\
\n\
The default value is returned when an invalid key is used.\n\
"                                                                                                                );

/* See comment in xxsubtype.c */
#  define DEFERRED_ADDRESS(ADDR) 0

static PyTypeObject defdict_type = {
	PyVarObject_HEAD_INIT(DEFERRED_ADDRESS(&PyType_Type),		   0)
	"simuPOP.defaultdict",                                          /* tp_name */
	sizeof(defdictobject),                                          /* tp_basicsize */
	0,                                                              /* tp_itemsize */
	/* methods */
	(destructor)defdict_dealloc,                                    /* tp_dealloc */
	0,                                                              /* tp_print */
	0,                                                              /* tp_getattr */
	0,                                                              /* tp_setattr */
	0,                                                              /* tp_reserved */
	(reprfunc)defdict_repr,                                         /* tp_repr */
	0,                                                              /* tp_as_number */
	0,                                                              /* tp_as_sequence */
	0,                                                              /* tp_as_mapping */
	0,                                                              /* tp_hash */
	0,                                                              /* tp_call */
	0,                                                              /* tp_str */
	PyObject_GenericGetAttr,                                        /* tp_getattro */
	0,                                                              /* tp_setattro */
	0,                                                              /* tp_as_buffer */
	Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC,
	/* tp_flags */
	defdict_doc,                                                    /* tp_doc */
	defdict_traverse,                                               /* tp_traverse */
	(inquiry)defdict_tp_clear,                                      /* tp_clear */
	0,                                                              /* tp_richcompare */
	0,                                                              /* tp_weaklistoffset*/
	0,                                                              /* tp_iter */
	0,                                                              /* tp_iternext */
	defdict_methods,                                                /* tp_methods */
	defdict_members,                                                /* tp_members */
	0,                                                              /* tp_getset */
	DEFERRED_ADDRESS(&PyDict_Type),                                 /* tp_base */
	0,                                                              /* tp_dict */
	0,                                                              /* tp_descr_get */
	0,                                                              /* tp_descr_set */
	0,                                                              /* tp_dictoffset */
	defdict_init,                                                   /* tp_init */
	PyType_GenericAlloc,                                            /* tp_alloc */
	0,                                                              /* tp_new */
	PyObject_GC_Del,                                                /* tp_free */
};

bool is_defdict(PyTypeObject * type)
{
	return type == &defdict_type;
}


PyObject * PyDefDict_New()
{
	defdictobject * obj;

	// This should call PyDict_Type.tp_new and create an object
	obj = (defdictobject *)defdict_type.tp_new((PyTypeObject *)(&defdict_type), NULL, NULL);
	if (obj == NULL) {
		PyObject_Del(obj);
		return PyErr_NoMemory();
	}
	// initialize this object (call PyDict_Type.tp_init)
	PyObject * args = PyTuple_New(0);
	PyDict_Type.tp_init((PyObject *)obj, args, NULL);
	Py_DECREF(args);
	return (PyObject *)obj;
}


int initCustomizedTypes(void)
{
	Py_TYPE(&Arraytype) = &PyType_Type;
	if (PyType_Ready(&Arraytype) < 0)
		return -1;
	//
	Py_TYPE(&defdict_type) = &PyType_Type;
	defdict_type.tp_base = &PyDict_Type;

	if (PyType_Ready(&defdict_type) < 0)
		return -1;
	//Py_INCREF(&defdict_type);
	return 0;
}


#endif
