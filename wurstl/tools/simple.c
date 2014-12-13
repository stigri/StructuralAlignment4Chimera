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
#include "pair_set_i.h"
#include "pair_set_p_i.h"
#include "pair_set.h"
#include "pdbin_i.h"
#include "pdbout_i.h"
#include "lsqf.h"
#include "read_seq_i.h"
#include "score_alnm.h"
#include "pw_alignment.h"
#include "pair_set_sel.h"

const int N_AND_W = 0;


/* ---------------- usage  ------------------------------------
 */
static void
usage(char *progname)
{
    static const char *errmsg =
        "pdb1 pdb2 pdbout [-a alg_type] [-r rmsd_threshold]\n\
    pdb1, pdb2 file names of coordinate files\n\
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
    char *pdbout;
    int alg_type;
    float rmsd_thresh;

    struct coord *coord1, *coord2;
    struct seq *seq1, *seq2;
    struct pair_set *set_alg;
    struct score_struct *scores;
    float rmsd = 0;
    struct algnm_param * params;
    int optn;
    size_t n = 0;
    char *pbs;
    
    if (argc < 8) {
        usage(argv[0]);
        return EXIT_FAILURE;
    }

    path1= argv[1];
    path2 = argv[2];
    pdbout = argv[3];
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

    /*get and set parameters*/
    params = init_algnm_param();
    params->alg_type = alg_type;
    params->rmsd_thresh = rmsd_thresh;
    /*  Alignment computation */
    set_alg = pw_algnt(coord1, coord2, params, &rmsd, &n, &pbs);
    mprintf("simple>> pbs is %s\n", pbs);

    /* Alignment evaluation  */
    scores = get_scores(set_alg, coord1, coord2, params->sw1_pgap_open, params->sw1_pgap_widen);

    seq1 = coord_get_seq(coord1);
    if (coord_2_pdb(pdbout, coord1, seq1) == EXIT_FAILURE) {
        err_printf (argv[0], "Fail writing coordinates to %s\n", pdbout);
        return EXIT_FAILURE;
    }
    mprintf ("scores evaluation = %f RMSD = %f\n", scores->scr_tot, rmsd);
    mprintf ("%ld pairs took part at the superposition. (uppercase) \n",
             n);
//     mprintf ("S1: %s \n", seq_print(seq1));
    /*converts characters up-/lowcase (superposited/not) */
    seq2 = coord_get_seq(coord2);
    pair_set_print_prepare(seq1, seq2, set_alg, pbs);
    if (alg_type == N_AND_W) {
        pair_set_extend(set_alg, seq_size(seq1), seq_size(seq2), EXT_LONG);
    }

    mprintf ("%s \n", pair_set_pretty_string(set_alg, seq1, seq2, NULL, NULL));
//     mprintf ("S2: %s \n", seq_print(seq2));

    pair_set_destroy(set_alg);
    coord_destroy(coord1);
    coord_destroy(coord2);
    free(scores);
    bit_set_free(pbs);
    
    return EXIT_SUCCESS;
}
