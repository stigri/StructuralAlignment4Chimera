/*
 * pwrapper.c
 *
 *  Created on: Dec 10, 2014
 *      Author: stine
 */
#include <Python.h>

static PyObject* structualignment(PyObject* self, PyObject* args)
{
	char* s;

	if (!PyArg_ParseTuple(args, "s", &s)) {
		return NULL;
	}

    return Py_BuildValue("s", s);
}

static char structualignment_docs[] =
    "structualignment( ): Any message you want to put here!!\n";

static PyMethodDef libwurstl_funcs[] = {
    {"structualignment", (PyCFunction)structualignment,
     METH_VARARGS, structualignment_docs},
    {NULL}
};

void initlibwurstl(void)
{
    Py_InitModule3("libwurstl", libwurstl_funcs,
                   "Extension module example!");
}

