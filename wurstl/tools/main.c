/*
 * main.c
 *
 *  Created on: Dec 8, 2010
 *      Author: ibondarenko
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <strings.h>
#include <limits.h>
#include <unistd.h>
#include <getopt.h>

#include "align_i.h"
#include "e_malloc.h"
#include "mprintf.h"
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
#include "pdbin_i.h"
#include "pdbout_i.h"
#include "lsqf.h"
#include "pair_set_sel.h"
#include "score_alnm.h"

const float sw1_pgap_open = 3.25;
const float sw1_pgap_widen = 0.8942;
const float mat_shift = -0.1;
const unsigned char mag_num_asrest = 50;
const int N_AND_W = 0;


/* ---------------- usage  ------------------------------------
 */
static void
usage(char *progname)
{
    static const char *errmsg =
        "pdb1 pdb2 classfile pdbout [-a alg_type] [-r rmsd_threshold]\n\
    pdb1, pdb2 file names of coordinate files\n\
    classfile - classification file\n\
    alg_type - 0 means smith and waterman, 1 means needleman & wunsch (default 0)\n\
    rmsd_threshold in Angstrom (default 3.0)\n";
    err_printf (progname, errmsg);
}

/* ---------------- pdb_read_with_err -------------------------
 * Wrapper around pdb_read, to print out a decent error message
 * At the moment, we do not allow one to pick which chain is read
 * from the file. This is fundamentally wrong and should be fixed.
 */
static struct coord *
pdb_read_with_err (const char *fname)
{
    static const char *acq_c = "";
    const char chain = '_';
    static const char *this_sub = "pdb_read_with_err";
    struct coord *r = pdb_read (fname, acq_c, chain);
    if (r == NULL)
        err_printf (this_sub, "Fail coordinate file %s\n", fname);
    return r;
}

/* ---------------- main   ------------------------------------
 */
int
main (int argc, char *argv[])
{
    char *path1;
    char *path2;
    char *filename;
    char *pdbout;
    int alg_type;
    int res;
    float rmsd_thresh;

    FILE *classfile;
    struct coord *coord1, *coord2;
    struct aa_strct_clssfcn *class;
    struct prob_vec *pvec1, *pvec2;
    struct seq *seq1, *seq2;
    struct score_mat  *crap = NULL;
    struct pair_set *set_alg, *set_alg_tmp;
    struct score_struct *scores;
    struct score_mat  *matrix;
    float rmsd;
    char *pbs;
    int optn;
    
    if (argc < 5) {
        usage(argv[0]);
        return EXIT_FAILURE;
    }

    path1= argv[1];
    path2 = argv[2];
    filename = argv[3];
    pdbout = argv[4];
    alg_type = 1;                   /* default smith waterman */
    rmsd_thresh = 3;                /* default 3 Angstrom     */
    
    
    while ((optn = getopt(argc, argv, "a:r:")) != -1) {
        switch (optn) {
        case 'a':
            alg_type = atoi(optarg); 
            break;
        case 'r': 
            rmsd_thresh = atof(optarg);
            break;
        }
    }
    
    mprintf("alg_type is %d, rmsd_threshold is %f\n", alg_type, rmsd_thresh);
    /* read coordinates from pdb-file.*/
    if ((coord1 = pdb_read_with_err(path1)) == NULL)
        return EXIT_FAILURE;
    if ((coord2 = pdb_read_with_err(path2)) == NULL)
        return EXIT_FAILURE;

    /* read the classification structure */
    classfile = fopen(filename, "r");
    if (classfile == NULL) {
        err_printf(argv[0], "Open fail on classfile %s\n", filename);
        return EXIT_FAILURE;
    }
    fclose(classfile);
    
    class = aa_strct_clssfcn_read (filename, 0.4);

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
    	err_printf(argv[0], "Can not get any selected pairs from alignment.\n");
    	return EXIT_FAILURE;
    }
    do {
        if (coord_rmsd(set_alg_tmp, coord1, coord2, 0, &rmsd, &coord1, &coord2) == EXIT_SUCCESS) {
            /*remove max rmsd paar. */
            pbs = pair_set_sel_delmaxdistance(coord1, coord2, set_alg, pbs);
            pair_set_destroy(set_alg_tmp);
            set_alg_tmp = selected_pair_set_get(set_alg, pbs);
        } else {
            err_printf(argv[0], "coord_rmsd was not successful. \n");
            break;
        }
    } while ((rmsd > rmsd_thresh) && (set_alg_tmp->n > 50));

    if (coord_2_pdb(pdbout, coord1, seq1) == EXIT_FAILURE) {
        err_printf (argv[0], "Fail writing coordinates to %s\n", pdbout);
        return EXIT_FAILURE;
    }
    mprintf ("score_pvec = %d scores evaluation = %f RMSD = %f\n",
             res, scores->scr_tot, rmsd);
    mprintf ("%ld pairs took part at the superposition. (uppercase) \n",
             set_alg_tmp->n);
//     mprintf ("S1: %s \n", seq_print(seq1));
    /*converts characters up-/lowcase (superposited/not) */
    pair_set_print_prepare(seq1, seq2, set_alg, pbs);
    if (alg_type == N_AND_W) {
        pair_set_extend(set_alg, seq_size(seq1), seq_size(seq2), EXT_LONG);
    }

    mprintf ("%s \n", pair_set_pretty_string(set_alg, seq1, seq2, NULL, NULL));
//     mprintf ("S2: %s \n", seq_print(seq2));

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

    return EXIT_SUCCESS;
}
