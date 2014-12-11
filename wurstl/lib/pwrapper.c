/*
 * pwrapper.c
 *
 *  Created on: Dec 10, 2014
 *      Author: stine
 */
#include <Python.h>
#include <string.h>

#include "align_i.h"
#include "e_malloc.h"
#include "coord.h"
#include "coord_i.h"
#include "read_ac_strct_i.h"
#include "prob_vec_i.h"
#include "score_mat_i.h"
#include "read_seq_i.h"
#include "score_probvec.h"
#include "pair_set_i.h"
#include "pair_set_p_i.h"
#include "pair_set.h"
#include "yesno.h"
#include "pdbin.h"
#include "pdbin_i.h"
#include "pdbout_i.h"
#include "lsqf.h"
#include "pair_set_sel.h"
#include "score_alnm.h"

const float sw1_pgap_open = 3.25;
const float sw1_pgap_widen = 0.8942;
const float mat_shift = -0.1;


/* This function is an almost exact copy of tools/main.c (except for exception handling) */
static int structural_alignment_internal(struct coord *coord1, struct coord *coord2) {

    char *filename = "/home/stine/src/StructuralAlignment4Chimera/wurstl/examples/classfile";
    int alg_type;
    int res;
    float rmsd_thresh;

    FILE *classfile;
    struct aa_strct_clssfcn *class;
    struct prob_vec *pvec1, *pvec2;
    struct seq *seq1, *seq2;
    struct score_mat  *crap = NULL;
    struct pair_set *set_alg, *set_alg_tmp;
    struct score_struct *scores;
    struct score_mat  *matrix;
    float rmsd;
    char *pbs;

    alg_type = 1;                   /* default smith waterman */
    rmsd_thresh = 3;                /* default 3 Angstrom     */

    class = aa_strct_clssfcn_read (filename, 0.4);
    if(class == NULL) {
    	PyErr_SetString(PyExc_IOError,"Opening classifier file failed");
    	return EXIT_FAILURE;
    }

    pvec1 = strct_2_prob_vec(coord1, class, YES );
    pvec2 = strct_2_prob_vec(coord2, class, YES );
    seq1 = coord_get_seq(coord1);
    seq2 = coord_get_seq(coord2);

    /* DP Matrix computation  */
    matrix = score_mat_new(seq_size(seq1), seq_size(seq2));
    matrix = score_mat_shift(matrix, mat_shift);
    res = score_pvec(matrix, pvec1, pvec2);

    /*  Alignment computation */
    set_alg = score_mat_sum_full(&crap, matrix,
                             sw1_pgap_open, sw1_pgap_widen,
                             sw1_pgap_open, sw1_pgap_widen,
                             NULL, NULL, alg_type, NULL);

    /* Alignment evaluation  */
    scores = get_scores(set_alg, coord1, coord2, sw1_pgap_open, sw1_pgap_widen);

    /*superimposed residues selection*/
    pbs = pair_set_sel_geti(set_alg);
    set_alg_tmp = selected_pair_set_get(set_alg, pbs);
    if (set_alg_tmp == NULL) {
    	//TODO error/exception handling
    	//err_printf(argv[0], "Can not get any selected pairs from alignment.\n");
    	return EXIT_FAILURE;
    }
    do {
        if (coord_rmsd(set_alg_tmp, coord1, coord2, 0, &rmsd, &coord1, &coord2) == EXIT_SUCCESS) {
            /*remove max rmsd paar. */
            pbs = pair_set_sel_delmaxdistance(coord1, coord2, set_alg, pbs);
            pair_set_destroy(set_alg_tmp);
            set_alg_tmp = selected_pair_set_get(set_alg, pbs);
        } else {
        	//TODO error/exception handling
            //err_printf(argv[0], "coord_rmsd was not successful. \n");
            break;
        }
    } while ((rmsd > rmsd_thresh) && (set_alg_tmp->n > 50));

    /* Just for debug/verification purposes
    mprintf ("score_pvec = %d scores evaluation = %f RMSD = %f\n",
             res, scores->scr_tot, rmsd);
    mprintf ("%ld pairs took part at the superposition. (uppercase) \n",
             set_alg_tmp->n);
    pair_set_print_prepare(seq1, seq2, set_alg, pbs);
    mprintf ("%s \n", pair_set_pretty_string(set_alg, seq1, seq2, NULL, NULL));
	*/

    score_mat_destroy(crap);
    pair_set_destroy(set_alg);
    coord_destroy(coord1);
    coord_destroy(coord2);
    /*pair_set_destroy(set_alg_tmp);*/
    free(scores);
    aa_strct_clssfcn_destroy(class);
    score_mat_destroy(matrix);
    prob_vec_destroy(pvec1);
    prob_vec_destroy(pvec2);
}

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

	for (i = 0; i < len; i++) {
    	item = (PyObject*) PyList_GET_ITEM(seq, i);
    	if ((long)PyList_Size(item) == 9) {

    		atom = nothing_atom();

    		atom.at_num = i; // just map atom location to the position in the list
    		strncpy(atom.at_name, PyString_AsString(PyList_GET_ITEM(item,0)), 5);
    		//ID Atom Type not used right now!!!
    		//idam_type = PyString_AsString(PyList_GET_ITEM(item,1));
    		strncpy(atom.res_name, PyString_AsString(PyList_GET_ITEM(item,2)), 4);
    		atom.icode = (char) PyString_AsString(PyList_GET_ITEM(item,3));
    		atom.chain_id = PyString_AsString(PyList_GET_ITEM(item,4));
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

    structural_alignment_internal(cA, cB	);

    return Py_BuildValue("s", compndA);
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

