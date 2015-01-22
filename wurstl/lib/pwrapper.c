/*
 * pwrapper.c
 *
 *  Created on: Dec 10, 2014
 *      Author: stine
 */
#include <Python.h>
#include <string.h>

#include "coord.h"
#include "pair_set.h"
#include "yesno.h"
#include "pdbin.h"
#include "pw_alignment.h"

static struct allatoms get_all_atoms(PyObject* obj)
{
	/* Vars for all atoms */
	struct allatoms allatoms;
	struct atompdb atom = nothing_atom();
	struct atompdb *atoms;

	/* Vars needed for atomS list conversion */
    PyObject* seq;
    int i, len;

    /* Vars needed for individual atom conversion */
    PyObject* item;

    seq = PySequence_Fast(obj, "expected a sequence");
    len = PySequence_Size(obj);

	atoms = malloc(len * sizeof(nothing_atom()));
	allatoms.atoms = atoms;
	allatoms.n = len;

	/* http://www.wwpdb.org/documentation/format33/sect9.html#ATOM */

	for (i = 0; i < len; i++) {
    	item = (PyObject*) PyList_GET_ITEM(seq, i);
    	if ((long)PyList_Size(item) == 9) {

    		atom = nothing_atom();

    		atom.at_num = i; // just map atom location to the position in the list
    		strncpy(atom.at_name, PyString_AsString(PyList_GET_ITEM(item,0)), 5);
    		//ID Atom Type not used right now!!!
    		//idam_type = PyString_AsString(PyList_GET_ITEM(item,1));
    		strncpy(atom.res_name, PyString_AsString(PyList_GET_ITEM(item,2)), 4);
    		atom.icode = PyString_AsString(PyList_GET_ITEM(item,3))[0];
    		atom.chain_id = PyString_AsString(PyList_GET_ITEM(item,4))[0];
    		atom.res_num = (int) PyInt_AsLong(PyList_GET_ITEM(item,5)); // Correct mapping?
    		atom.x = PyFloat_AsDouble(PyList_GET_ITEM(item,6));
    		atom.y = PyFloat_AsDouble(PyList_GET_ITEM(item,7));
    		atom.z = PyFloat_AsDouble(PyList_GET_ITEM(item,8));

    		atoms[i] = atom;
    	}
    }
    Py_DECREF(seq);
    return allatoms;
}

static PyObject* structural_alignment(PyObject* self, PyObject* args)
{
	PyObject* obj;

	/* params for pw_align */
	struct pair_set* set_alg;
    float rmsd = 0;
	struct algnm_param * params;
    size_t n = 0;
    char *pbs;

    /* set_alg to PyObject conversion */
	PyObject* result = NULL;
	size_t idx;
	int **p;

	/* model A */
	struct allatoms allatomsA;
	struct coord *cA = NULL;
	char* compndA;
    PyObject* atomsA;

	/* model B */
	struct allatoms allatomsB;
	struct coord *cB = NULL;
    char* compndB;
    PyObject* atomsB;

	/* Immediately flush stdout (helpful during debugging) */
    //setbuf(stdout, NULL);

    /* convert python parameter types into C counterparts */
	if (!PyArg_ParseTuple(args, "O", &obj)) {
		return NULL;
	}

	/*
	 * WARNING: This code only uses two models stored in the dynamic structure obj,
	 * even though the python layer could pass many more models
	 */

	/* for first model */
	compndA = PyString_AsString(PyList_GetItem(obj, 0));
    atomsA = PyList_GET_ITEM(obj, 1);

    allatomsA = get_all_atoms(atomsA);
    cA = atoms2mdl(allatomsA);
    cA->compnd = compndA;

    /* for second model */
	compndB = PyString_AsString(PyList_GetItem(obj, 2));
    atomsB = PyList_GET_ITEM(obj, 3);

    allatomsB = get_all_atoms(atomsB);
    cB = atoms2mdl(allatomsB);
    cB->compnd = compndB;

    /*
     * At this point all python objects have been copied to local
     * variables. Thus, we can release Python's global interpreter lock (GIL)
     * to allow the calling Chimera extension to do its stuff (e.g. painting the UI).
     *
     * Also note the Py_END_ALLOW_THREADS macro below. It means we will re-acquire
     * GIL.
     *
     * see https://docs.python.org/2/c-api/init.html#releasing-the-gil-from-extension-code
     */
    Py_BEGIN_ALLOW_THREADS;

    /*get and set parameters*/
    params = init_algnm_param();
    params->alg_type = 1; /* default smith waterman */
    params->rmsd_thresh = 3; /* default 3 Angstrom */

    /*
     * Here's where the beef is.
     * (we ignore rmsd, n and pbs)
     */
    set_alg = pw_algnt(cA, cB, params, &rmsd, &n, &pbs);

    /*
     * Free stuff that has been allocated by *us* earlier on.
     */
    free(allatomsA.atoms);
    free(allatomsB.atoms);

    Py_DECREF(atomsA);
    Py_DECREF(atomsB);

    /*
     * Convert wurstl structs to python-consumable data.
     * WARNING: This code (incorrectly) assumes that the
     * m-dimension of pair_set is fixed to 2, meaning there
     * are only to sequences. This assumption is borrowed
     * from other code passages that indicate that wurstl
     * in general cannot handle more than two sequences yet.
     */
    p = set_alg->indices;
    result = PyList_New(set_alg->n);
   	for (idx = 0; idx < set_alg->n; idx++) {
   		PyList_SetItem(result, idx, Py_BuildValue("(ii)", p[idx][0], p[idx][1]));
   	}

   	/* Destroy the pair_set after it has been converted to the result */
   	pair_set_destroy(set_alg);

   	Py_END_ALLOW_THREADS;

   	return result;
}

/*
 * see https://docs.python.org/2/extending/extending.html#writing-extensions-in-c
 * see https://docs.python.org/2/c-api/arg.html
 * see https://docs.python.org/2/c-api/concrete.html
 * see http://www.tutorialspoint.com/python/python_further_extensions.htm
 */
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

