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
typedef struct arrayobject_template<GenoIterator> arrayobject;                                  /* Forward */

#  if PY_VERSION_HEX < 0x02060000
#    define Py_SIZE(obj) (((PyVarObject *)(obj))->ob_size)
#  endif


/// CPPONLY
bool is_carrayobject(PyObject * op);

/// CPPONLY
PyObject *
getarrayitem(arrayobject * op, Py_ssize_t i)
{
	return getarrayitem_template<GenoIterator>(op, i);
}


/// CPPONLY
int
setarrayitem(arrayobject * ap, Py_ssize_t i, PyObject * v)
{
	return setarrayitem_template<GenoIterator>(ap, i, v);
}


/// CPPONLY
PyObject *
carray_new(PyTypeObject * a, PyObject * b, PyObject * c)
{
	return carray_new_template<GenoIterator>(a, b, c);
}


/// CPPONLY
PyObject *
carray_init(PyTypeObject * a, PyObject * b, PyObject * c)
{
	return carray_init_template<GenoIterator>(a, b, c);
}


/// CPPONLY
PyObject * newcarrayobject(GenoIterator begin, GenoIterator end);

/// CPPONLY
void
array_dealloc(arrayobject * op)
{
	array_dealloc_template<GenoIterator>(op);
}


/// CPPONLY
PyObject *
array_richcompare(PyObject * v, PyObject * w, int op)
{
	return(array_richcompare_template<GenoIterator>(v, w, op));
}


/// CPPONLY
Py_ssize_t array_length(arrayobject * a)
{
	return(array_length_template<GenoIterator>(a));
}


/// CPPONLY
PyObject * array_concat(arrayobject * a, PyObject * o)
{
	return(array_concat_template<GenoIterator>(a, o));
}


/// CPPONLY
PyObject * array_repeat(arrayobject * a, Py_ssize_t i)
{
	return(array_repeat_template<GenoIterator>(a, i));
}


/// CPPONLY
PyObject * array_item(arrayobject * a, Py_ssize_t i)
{
	return(array_item_template<GenoIterator>(a, i));
}


/// CPPONLY
PyObject * array_slice(arrayobject * a, Py_ssize_t ilow, Py_ssize_t ihigh)
{
	return(array_slice_template<GenoIterator>(a, ilow, ihigh));
}


/// CPPONLY
int array_ass_slice(arrayobject * a, Py_ssize_t ilow, Py_ssize_t ihigh, PyObject * v)
{
	return(array_ass_slice_template<GenoIterator>(a, ilow, ihigh, v));
}


/// CPPONLY
Py_ssize_t array_ass_item(arrayobject * a, Py_ssize_t i, PyObject * v)
{
	return(array_ass_item_template<GenoIterator>(a, i, v));
}


/// CPPONLY
PyObject * array_count(arrayobject * self, PyObject * args)
{
	return(array_count_template<GenoIterator>(self, args));
}


/// CPPONLY
char count_doc [] =
    "count(x)\n\
\n\
Return number of occurences of x in the array."                                          ;

/// CPPONLY
PyObject * array_index(arrayobject * self, PyObject * args)
{
	return(array_index_template<GenoIterator>(self, args));
}


char index_doc [] =
    "index(x, [start, [stop]])\n\
\n\
Return index of first occurence of x in the array.";

/// CPPONLY
PyObject * array_tolist(arrayobject * self, PyObject * args)
{
	return(array_tolist_template<GenoIterator>(self, args));
}


char tolist_doc [] =
    "tolist() -> list\n\
\n\
Convert array to an ordinary list with the same items."                                                          ;

#pragma GCC diagnostic ignored "-Wmissing-field-initializers"
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
PyObject * array_getattr(arrayobject * a, char * name)
{
	if (strcmp(name, "__members__") == 0) {
		PyObject * list = PyList_New(0);
		return list;
	}
	return Py_FindMethod(array_methods, (PyObject *)a, name);
}


/// CPPONLY
int array_print(arrayobject * a, FILE * fp, int flags)
{
	return(array_print_template<GenoIterator>(a, fp, flags));
}


/// CPPONLY
PyObject *
array_repr(arrayobject * a)
{
	return(array_repr_template<GenoIterator>(a));
}


PySequenceMethods array_as_sequence =
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

char arraytype_doc [] =
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
	return is_carrayobject_template<GenoIterator>(op);
}


/// CPPONLY
PyObject * newcarrayobject(GenoIterator begin, GenoIterator end)
{
	return(newcarrayobject_template<GenoIterator>(begin, end));
}

/* lineage array type ***************************/

/// CPPONLY
typedef struct arrayobject_template<LineageIterator> arrayobject_lineage;                      /* Forward */


/// CPPONLY
bool is_carrayobject_lineage(PyObject * op);

/// CPPONLY
PyObject *
getarrayitem_lineage(arrayobject_lineage * op, Py_ssize_t i)
{
	return getarrayitem_template<LineageIterator>(op, i);
}


/// CPPONLY
int
setarrayitem_lineage(arrayobject_lineage * ap, Py_ssize_t i, PyObject * v)
{
	return setarrayitem_template<LineageIterator>(ap, i, v);
}


/// CPPONLY
PyObject *
carray_new_lineage(PyTypeObject * a, PyObject * b, PyObject * c)
{
	return carray_new_template<LineageIterator>(a, b, c);
}


/// CPPONLY
PyObject *
carray_init_lineage(PyTypeObject * a, PyObject * b, PyObject * c)
{
	return carray_init_template<LineageIterator>(a, b, c);
}


/// CPPONLY
PyObject * newcarrayobject_lineage(LineageIterator begin, LineageIterator end);

/// CPPONLY
void
array_dealloc_lineage(arrayobject_lineage * op)
{
	array_dealloc_template<LineageIterator>(op);
}


/// CPPONLY
PyObject *
array_richcompare_lineage(PyObject * v, PyObject * w, int op)
{
	return(array_richcompare_template<LineageIterator>(v, w, op));
}


/// CPPONLY
Py_ssize_t array_length_lineage(arrayobject_lineage * a)
{
	return(array_length_template<LineageIterator>(a));
}


/// CPPONLY
PyObject * array_concat_lineage(arrayobject_lineage * a, PyObject * o)
{
	return(array_concat_template<LineageIterator>(a, o));
}


/// CPPONLY
PyObject * array_repeat_lineage(arrayobject_lineage * a, Py_ssize_t i)
{
	return(array_repeat_template<LineageIterator>(a, i));
}


/// CPPONLY
PyObject * array_item_lineage(arrayobject_lineage * a, Py_ssize_t i)
{
	return(array_item_template<LineageIterator>(a, i));
}


/// CPPONLY
PyObject * array_slice_lineage(arrayobject_lineage * a, Py_ssize_t ilow, Py_ssize_t ihigh)
{
	return(array_slice_template<LineageIterator>(a, ilow, ihigh));
}


/// CPPONLY
int array_ass_slice_lineage(arrayobject_lineage * a, Py_ssize_t ilow, Py_ssize_t ihigh, PyObject * v)
{
	return(array_ass_slice_template<LineageIterator>(a, ilow, ihigh, v));
}


/// CPPONLY
Py_ssize_t array_ass_item_lineage(arrayobject_lineage * a, Py_ssize_t i, PyObject * v)
{
	return(array_ass_item_template<LineageIterator>(a, i, v));
}


/// CPPONLY
PyObject * array_count_lineage(arrayobject_lineage * self, PyObject * args)
{
	return(array_count_template<LineageIterator>(self, args));
}


/// CPPONLY
char count_doc_lineage [] =
    "count(x)\n\
\n\
Return number of occurences of x in the array."                                          ;

/// CPPONLY
PyObject * array_index_lineage(arrayobject_lineage * self, PyObject * args)
{
	return(array_index_template<LineageIterator>(self, args));
}


char index_doc_lineage [] =
    "index(x, [start, [stop]])\n\
\n\
Return index of first occurence of x in the array.";

/// CPPONLY
PyObject * array_tolist_lineage(arrayobject_lineage * self, PyObject * args)
{
	return(array_tolist_template<LineageIterator>(self, args));
}


char tolist_doc_lineage [] =
    "tolist() -> list\n\
\n\
Convert array to an ordinary list with the same items."                                                          ;

#pragma GCC diagnostic ignored "-Wmissing-field-initializers"
PyMethodDef array_methods_lineage[] =
{
	{
		"count", (PyCFunction)array_count_lineage, METH_VARARGS,
		count_doc_lineage
	},
	{
		"index", (PyCFunction)array_index_lineage, METH_VARARGS,
		index_doc_lineage
	},
	{
		"tolist",    (PyCFunction)array_tolist_lineage,    METH_VARARGS,
		tolist_doc_lineage
	},
	{                                                                                             /* sentinel */
		NULL,        NULL
	}
};

/// CPPONLY
PyObject * array_getattr_lineage(arrayobject_lineage * a, char * name)
{
	if (strcmp(name, "__members__") == 0) {
		PyObject * list = PyList_New(0);
		return list;
	}
	return Py_FindMethod(array_methods_lineage, (PyObject *)a, name);
}


/// CPPONLY
int array_print_lineage(arrayobject_lineage * a, FILE * fp, int flags)
{
	return(array_print_template<LineageIterator>(a, fp, flags));
}


/// CPPONLY
PyObject *
array_repr_lineage(arrayobject_lineage * a)
{
	return(array_repr_template<LineageIterator>(a));
}


PySequenceMethods array_as_sequence_lineage =
{
#  if PY_VERSION_HEX < 0x02050000
	(inquiry)array_length_lineage,                                                      /*sq_length*/
	(binaryfunc)array_concat_lineage,                                                   /*sq_concat*/
	(intargfunc)array_repeat_lineage,                                                   /*sq_repeat*/
	(intargfunc)array_item_lineage,                                                     /*sq_item*/
	(intintargfunc)array_slice_lineage,                                                 /*sq_slice*/
	(intobjargproc)array_ass_item_lineage,                                              /*sq_ass_item*/
	(intintobjargproc)array_ass_slice_lineage,                                          /*sq_ass_slice*/
#  else
	(lenfunc)array_length_lineage,                                                      /*sq_length*/
	(binaryfunc)array_concat_lineage,                                                   /*sq_concat*/
	(ssizeargfunc)array_repeat_lineage,                                                 /*sq_repeat*/
	(ssizeargfunc)array_item_lineage,                                                   /*sq_item*/
	(ssizessizeargfunc)array_slice_lineage,                                             /*sq_slice*/
	(ssizeobjargproc)array_ass_item_lineage,                                            /*sq_ass_item*/
	(ssizessizeobjargproc)array_ass_slice_lineage,                                      /*sq_ass_slice*/
#  endif
};

char arraytype_doc_lineage [] =
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

PyTypeObject LineageArraytype =
{
	PyObject_HEAD_INIT(NULL)
	0,
	"simuPOP.carray_lineage",   /* mudule.type name */
	sizeof(arrayobject_lineage),
	0,
	(destructor)array_dealloc_lineage,  /* tp_dealloc */
	(printfunc)array_print_lineage,     /* tp_print */
	(getattrfunc)array_getattr_lineage, /* tp_getattr */
	0,                          /* tp_setattr */
	0,                          /* tp_compare */
	(reprfunc)array_repr_lineage,       /* tp_repr */
	0,                          /* tp_as _number*/
	&array_as_sequence_lineage, /* tp_as _sequence*/
	0,                          /* tp_as _mapping*/
	0,                          /* tp_hash */
	0,                          /* tp_call */
	0,                          /* tp_str */
	0,                          /* tp_getattro */
	0,                          /* tp_setattro */
	0,                          /* tp_as_buffer*/
	Py_TPFLAGS_DEFAULT,         /* tp_flags */
	arraytype_doc_lineage,      /* tp_doc */
	0,                          /* tp_traverse */
	0,                          /* tp_clear */
	array_richcompare_lineage,  /* tp_richcompare */
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
	(initproc)carray_init_lineage,      /* tp_init */
	0,                          /* tp_alloc */
	carray_new_lineage,         /* tp_new */
};


/// CPPONLY
bool is_carrayobject_lineage(PyObject * op)
{
	return is_carrayobject_template<LineageIterator>(op);
}


/// CPPONLY
PyObject * newcarrayobject_lineage(LineageIterator begin, LineageIterator end)
{
	return(newcarrayobject_template<LineageIterator>(begin, end));
}


/* defdict type *********************************************************/

typedef struct
{
	PyDictObject dict;
} defdictobject;

//PyTypeObject defdict_type; /* Forward */

PyObject *
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

PyMethodDef defdict_methods[] = {
	{ "__getitem__", (PyCFunction)dict_subscript, METH_O | METH_COEXIST,
	  getitem__doc__ },
	{ NULL }
};


PyMemberDef defdict_members[] = {
	{ NULL }
};

int
defdict_print(defdictobject * dd, FILE * fp, int /* flags */)
{
	int sts;

	fprintf(fp, "defdict(");
	sts = PyDict_Type.tp_print((PyObject *)dd, fp, 0);
	fprintf(fp, ")");
	return sts;
}


int
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

void
defdict_dealloc(defdictobject * dd)
{
	PyDict_Type.tp_dealloc((PyObject *)dd);
}


int
defdict_traverse(PyObject * self, visitproc visit, void * arg)
{
	return PyDict_Type.tp_traverse(self, visit, arg);
}


int
defdict_tp_clear(defdictobject * dd)
{
	return PyDict_Type.tp_clear((PyObject *)dd);
}


// will change one of them to my...
PyMappingMethods defdict_as_mapping = *PyDict_Type.tp_as_mapping;

PyTypeObject defdict_type = {
	PyObject_HEAD_INIT(DEFERRED_ADDRESS(&PyType_Type))
	0,                              /* ob_size */
	"simuPOP.defdict",              /* tp_name */
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
int initCustomizedTypes(PyObject * m)
{
	// this will be done in PyType_Ready() is your read this
	// from python reference manual.
	Arraytype.ob_type = &PyType_Type;
	LineageArraytype.ob_type = &PyType_Type;
	if (PyType_Ready(&Arraytype) < 0 || PyType_Ready(&LineageArraytype) < 0)
		return -1;
	//
	defdict_type.ob_type = &PyType_Type;
	defdict_type.tp_base = &PyDict_Type;
	defdict_type.tp_as_mapping->mp_subscript = (binaryfunc)dict_subscript;

	if (PyType_Ready(&defdict_type) < 0)
		return -1;

	Py_INCREF(&defdict_type);
	if (PyModule_AddObject(m, "defdict", (PyObject *)&defdict_type) < 0)
		return -1;
	return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////

#else  // for Python 3

/* Array object implementation */

typedef struct arrayobject_template<GenoIterator> arrayobject;

bool is_carrayobject(PyObject * op);

PyObject * newcarrayobject(GenoIterator begin, GenoIterator end);

PyObject *
getarrayitem(PyObject * op, Py_ssize_t i)
{
	return getarrayitem_template<GenoIterator>(op, i);
}


/// CPPONLY
int
setarrayitem(arrayobject * ap, int i, PyObject * v)
{
	return setarrayitem_template<GenoIterator>(ap, i, v);
}

/* Methods */

void
array_dealloc(arrayobject * op)
{
	array_dealloc_template<GenoIterator>(op);
}


PyObject *
array_richcompare(PyObject * v, PyObject * w, int op)
{
	return array_richcompare_template<GenoIterator>(v, w, op);

}


Py_ssize_t
array_length(arrayobject * a)
{
	return array_length_template<GenoIterator>(a);
}


PyObject *
array_item(arrayobject * a, Py_ssize_t i)
{
	return array_item_template<GenoIterator>(a, i);
}


PyObject *
array_slice(arrayobject * a, Py_ssize_t ilow, Py_ssize_t ihigh)
{
	return array_slice_template<GenoIterator>(a, ilow, ihigh);
}


int
array_ass_slice(arrayobject * a, Py_ssize_t ilow, Py_ssize_t ihigh, PyObject * v)
{
	return array_ass_slice_template<GenoIterator>(a, ilow, ihigh, v);
}


int
array_ass_item(arrayobject * a, Py_ssize_t i, PyObject * v)
{
	return array_ass_item_template<GenoIterator>(a, i, v);
}


PyObject *
array_count(arrayobject * self, PyObject * v)
{
	return array_count_template<GenoIterator>(self, v);
}


PyDoc_STRVAR(count_doc,
	"count(x)\n\
\n\
Return number of occurrences of x in the array."                     );

PyObject *
array_index(arrayobject * self, PyObject * v)
{
	return array_index_template<GenoIterator>(self, v);
}


PyDoc_STRVAR(index_doc,
	"index(x)\n\
\n\
Return index of first occurrence of x in the array."                     );

PyObject *
array_tolist(arrayobject * self, PyObject * unused)
{
	return array_tolist_template<GenoIterator>(self, unused);
}


PyDoc_STRVAR(tolist_doc,
	"tolist() -> list\n\
\n\
Convert array to an ordinary list with the same items."                             );


PyMethodDef array_methods[] = {
	{ "count",	(PyCFunction)array_count,	METH_O,
	  count_doc },
	{ "index",	(PyCFunction)array_index,	METH_O,
	  index_doc },
	{ "tolist", (PyCFunction)array_tolist,	METH_NOARGS,
	  tolist_doc },
	{ NULL,		NULL }      /* sentinel */
};

PyObject *
array_repr(arrayobject * a)
{
	return array_repr_template<GenoIterator>(a);
}


PySequenceMethods array_as_sequence = {
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


PyObject*
array_subscr(arrayobject* self, PyObject* item)
{
	return array_subscr_template<GenoIterator>(self, item);
}


int
array_ass_subscr(arrayobject* self, PyObject* item, PyObject* value)
{
	return array_ass_subscr_template<GenoIterator>(self, item, value);
}

PyMappingMethods array_as_mapping = {
	(lenfunc)array_length,
	(binaryfunc)array_subscr,
	(objobjargproc)array_ass_subscr
};

PyObject * array_new(PyTypeObject * type, PyObject * args, PyObject * kwds)
{
	return array_new_template<GenoIterator>(type, args, kwds);
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

PyTypeObject Arraytype = {
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
	&array_as_mapping,                          /* tp_as_mapping*/
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
	return is_carrayobject_template<GenoIterator>(op);
}


/// CPPONLY
PyObject * newcarrayobject(GenoIterator begin, GenoIterator end)
{
	return newcarrayobject_template<GenoIterator>(begin, end);
}


/** lineage array type ******************* */

typedef struct arrayobject_template<LineageIterator> arrayobject_lineage;

bool is_carrayobject_lineage(PyObject * op);

PyObject * newcarrayobject_lineage(LineageIterator begin, LineageIterator end);

PyObject *
getarrayitem_lineage(PyObject * op, Py_ssize_t i)
{
	return getarrayitem_template<LineageIterator>(op, i);
}


/// CPPONLY
int
setarrayitem_lineage(arrayobject_lineage * ap, int i, PyObject * v)
{
	return setarrayitem_template<LineageIterator>(ap, i, v);
}

/* Methods */

void
array_dealloc_lineage(arrayobject_lineage * op)
{
	array_dealloc_template<LineageIterator>(op);
}


PyObject *
array_richcompare_lineage(PyObject * v, PyObject * w, int op)
{
	return array_richcompare_template<LineageIterator>(v, w, op);

}


Py_ssize_t
array_length_lineage(arrayobject_lineage * a)
{
	return array_length_template<LineageIterator>(a);
}


PyObject *
array_item_lineage(arrayobject_lineage * a, Py_ssize_t i)
{
	return array_item_template<LineageIterator>(a, i);
}


PyObject *
array_slice_lineage(arrayobject_lineage * a, Py_ssize_t ilow, Py_ssize_t ihigh)
{
	return array_slice_template<LineageIterator>(a, ilow, ihigh);
}


int
array_ass_slice_lineage(arrayobject_lineage * a, Py_ssize_t ilow, Py_ssize_t ihigh, PyObject * v)
{
	return array_ass_slice_template<LineageIterator>(a, ilow, ihigh, v);
}


int
array_ass_item_lineage(arrayobject_lineage * a, Py_ssize_t i, PyObject * v)
{
	return array_ass_item_template<LineageIterator>(a, i, v);
}


PyObject *
array_count_lineage(arrayobject_lineage * self, PyObject * v)
{
	return array_count_template<LineageIterator>(self, v);
}


PyDoc_STRVAR(count_doc_lineage,
	"count(x)\n\
\n\
Return number of occurrences of x in the array."                     );

PyObject *
array_index_lineage(arrayobject_lineage * self, PyObject * v)
{
	return array_index_template<LineageIterator>(self, v);
}


PyDoc_STRVAR(index_doc_lineage,
	"index(x)\n\
\n\
Return index of first occurrence of x in the array."                     );

PyObject *
array_tolist_lineage(arrayobject_lineage * self, PyObject * unused)
{
	return array_tolist_template<LineageIterator>(self, unused);
}


PyDoc_STRVAR(tolist_doc_lineage,
	"tolist() -> list\n\
\n\
Convert array to an ordinary list with the same items."                             );


PyMethodDef array_methods_lineage[] = {
	{ "count",	(PyCFunction)array_count_lineage,	METH_O,
	  count_doc_lineage },
	{ "index",	(PyCFunction)array_index_lineage,	METH_O,
	  index_doc_lineage },
	{ "tolist", (PyCFunction)array_tolist_lineage,	METH_NOARGS,
	  tolist_doc_lineage },
	{ NULL,		NULL }      /* sentinel */
};

PyObject *
array_repr_lineage(arrayobject_lineage * a)
{
	return array_repr_template<LineageIterator>(a);
}


PySequenceMethods array_as_sequence_lineage = {
	(lenfunc)array_length_lineage,                  /*sq_length*/
	0,                                      /*sq_concat*/
	0,                                      /*sq_repeat*/
	(ssizeargfunc)array_item_lineage,               /*sq_item*/
	0,                                      /*sq_slice*/
	(ssizeobjargproc)array_ass_item_lineage,        /*sq_ass_item*/
	0,					/*sq_ass_slice*/
	0,                                      /*sq_contains*/
	0,                                      /*sq_inplace_concat*/
	0                                       /*sq_inplace_repeat*/
};


PyObject*
array_subscr_lineage(arrayobject_lineage* self, PyObject* item)
{
	return array_subscr_template<LineageIterator>(self, item);
}


int
array_ass_subscr_lineage(arrayobject_lineage* self, PyObject* item, PyObject* value)
{
	return array_ass_subscr_template<LineageIterator>(self, item, value);
}

PyMappingMethods array_as_mapping_lineage = {
	(lenfunc)array_length_lineage,
	(binaryfunc)array_subscr_lineage,
	(objobjargproc)array_ass_subscr_lineage
};

PyObject * array_new_lineage(PyTypeObject * type, PyObject * args, PyObject * kwds)
{
	return array_new_template<LineageIterator>(type, args, kwds);
}


PyDoc_STRVAR(arraytype_doc_lineage,
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

PyTypeObject LineageArraytype = {
	PyVarObject_HEAD_INIT(NULL, 0)
	"simuPOP.array_lineage",
	sizeof(arrayobject_lineage),
	0,
	(destructor)array_dealloc_lineage,                  /* tp_dealloc */
	0,                                          /* tp_print */
	0,                                          /* tp_getattr */
	0,                                          /* tp_setattr */
	0,                                          /* tp_reserved */
	(reprfunc)array_repr_lineage,                       /* tp_repr */
	0,                                          /* tp_as_number*/
	&array_as_sequence_lineage,                         /* tp_as_sequence*/
	&array_as_mapping_lineage,                          /* tp_as_mapping*/
	0,                                          /* tp_hash */
	0,                                          /* tp_call */
	0,                                          /* tp_str */
	PyObject_GenericGetAttr,                    /* tp_getattro */
	0,                                          /* tp_setattro */
	0,                                          /* tp_as_buffer*/
	Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,   /* tp_flags */
	arraytype_doc_lineage,                              /* tp_doc */
	0,                                          /* tp_traverse */
	0,                                          /* tp_clear */
	array_richcompare_lineage,                          /* tp_richcompare */
	0,                                          /* tp_weaklistoffset */
	0,                                          /* tp_iter */
	0,                                          /* tp_iternext */
	array_methods_lineage,                              /* tp_methods */
	0,                                          /* tp_members */
	0,                                          /* tp_getset */
	0,                                          /* tp_base */
	0,                                          /* tp_dict */
	0,                                          /* tp_descr_get */
	0,                                          /* tp_descr_set */
	0,                                          /* tp_dictoffset */
	0,                                          /* tp_init */
	PyType_GenericAlloc,                        /* tp_alloc */
	array_new_lineage,                                  /* tp_new */
	PyObject_Del,                               /* tp_free */
};


bool is_carrayobject_lineage(PyObject * op)
{
	return is_carrayobject_template<LineageIterator>(op);
}


/// CPPONLY
PyObject * newcarrayobject_lineage(LineageIterator begin, LineageIterator end)
{
	return newcarrayobject_template<LineageIterator>(begin, end);
}


/**  defdict type ******************************/
typedef struct
{
	PyDictObject dict;
} defdictobject;

//PyTypeObject defdict_type; /* Forward */

PyDoc_STRVAR(defdict_missing_doc,
	"__missing__(key) # Called by __getitem__ for missing key; pseudo-code:\n\
  Return 0\n\
"                                                                                             );

PyObject *
defdict_missing(defdictobject * dd, PyObject * key)
{
	return PyInt_FromLong(0);
}


PyDoc_STRVAR(defdict_copy_doc, "D.copy() -> a shallow copy of D.");

PyObject *
defdict_copy(defdictobject * dd)
{
	/* This calls the object's class.  That only works for subclasses
	   whose class constructor has the same signature.  Subclasses that
	   define a different constructor signature must override copy().
	 */
	return PyObject_CallFunctionObjArgs((PyObject *)Py_TYPE(dd), Py_None, dd, NULL);
}


PyObject *
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
	result = PyTuple_Pack(5, Py_TYPE(dd),
		PyTuple_New(0), Py_None, Py_None, iter);
	Py_DECREF(iter);
	Py_DECREF(items);
	return result;
}


PyMethodDef defdict_methods[] = {
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

PyMemberDef defdict_members[] = {
	{ NULL }
};

void
defdict_dealloc(defdictobject * dd)
{
	PyDict_Type.tp_dealloc((PyObject *)dd);
}


PyObject *
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


int
defdict_traverse(PyObject * self, visitproc visit, void * arg)
{
	return PyDict_Type.tp_traverse(self, visit, arg);
}


int
defdict_tp_clear(defdictobject * dd)
{
	return PyDict_Type.tp_clear((PyObject *)dd);
}


int
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

PyTypeObject defdict_type = {
	PyVarObject_HEAD_INIT(DEFERRED_ADDRESS(&PyType_Type),		   0)
	"simuPOP.defdict",                                          /* tp_name */
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


int initCustomizedTypes(PyObject * m)
{
#if PY_VERSION_HEX >= 0x030b0000
	Py_SET_TYPE(&Arraytype, &PyType_Type);
	if (PyType_Ready(&Arraytype) < 0)
		return -1;
	//
	Py_SET_TYPE(&defdict_type, &PyType_Type);
	defdict_type.tp_base = &PyDict_Type;
#else
	Py_TYPE(&Arraytype) = &PyType_Type;
	if (PyType_Ready(&Arraytype) < 0)
		return -1;
	//
	Py_TYPE(&defdict_type) = &PyType_Type;
	defdict_type.tp_base = &PyDict_Type;
#endif
	if (PyType_Ready(&defdict_type) < 0)
		return -1;

	Py_INCREF(&defdict_type);
	PyModule_AddObject(m, "defdict", (PyObject *)&defdict_type);
	return 0;
}


#endif
