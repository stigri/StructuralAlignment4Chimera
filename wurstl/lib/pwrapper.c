/*
 * pwrapper.c
 *
 *  Created on: Dec 10, 2014
 *      Author: stine
 */
#include <Python.h>

static PyObject* structural_alignment(PyObject* self, PyObject* args)
{
	char* s;

	if (!PyArg_ParseTuple(args, "s", &s)) {
		return NULL;
	}

    return Py_BuildValue("s", s);
}

static char structural_alignment_docs[] =
    "structural_alignment( ): Any message you want to put here!!\n";

static PyMethodDef wurstl_funcs[] = {
    {"structural_alignment", (PyCFunction)structural_alignment,
     METH_VARARGS, structural_alignment_docs},
    {NULL}
};

void initwurstl(void)
{
    Py_InitModule3("wurstl", wurstl_funcs,
                   "WurstL extension module");
}

