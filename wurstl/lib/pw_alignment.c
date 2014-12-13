/*
 * main.c
 *
 *  Created on: Dec 11, 2014
 *      Author: ibondarenko
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "align_i.h"
#include "e_malloc.h"
#include "mprintf.h"
#include "coord.h"
#include "coord_i.h"
#include "pair_set_i.h"
#include "pair_set_p_i.h"
#include "pair_set.h"
#include "score_probvec.h"
#include "score_mat_i.h"
#include "pair_set_sel.h"
#include "yesno.h"
#include "read_ac_strct_i.h"
#include "read_seq_i.h"
#include "prob_vec_i.h"
#include "lsqf.h"
#include "const_clss.h"
#include "pw_alignment.h"
#include "bit_set.h"

struct algnm_param def_params = {
    3.25,      /*sw1_pgap_open*/
    0.8942,    /*sw1_pgap_widen*/
    -0.1,      /*mat_shift*/
    50,        /*mag_num_asrest*/
    3,         /*rmsd_thresh*/
    0          /*N_AND_W*/
};


struct algnm_param *
init_algnm_param (void) {
    return &def_params;
}


struct pair_set *
pw_algnt(struct coord *coord1, struct coord *coord2, struct algnm_param* param, float *rmsd_out, size_t *num, char **pbs) {
    
    struct aa_strct_clssfcn *class;
    struct prob_vec *pvec1, *pvec2;
    struct seq *seq1, *seq2;
    struct score_mat  *crap = NULL;
    struct pair_set *set_alg, *set_alg_tmp;
    struct score_mat  *matrix;
    
    class = paa_strct_clss6();

    pvec1 = strct_2_prob_vec(coord1, class, YES );
    pvec2 = strct_2_prob_vec(coord2, class, YES );
    seq1 = coord_get_seq(coord1);
    seq2 = coord_get_seq(coord2);

    /* DP Matrix computation  */
    matrix = score_mat_new(seq_size(seq1), seq_size(seq2));
    matrix = score_mat_shift(matrix, param->mat_shift);
    if (score_pvec(matrix, pvec1, pvec2) == EXIT_FAILURE) {
        err_printf("pw_algnt", "Fail to calculate a score matrix.\n");
        return 0;
    }

    /*  Alignment computation */
    set_alg = score_mat_sum_full(&crap, matrix,
                             param->sw1_pgap_open, param->sw1_pgap_widen,
                             param->sw1_pgap_open, param->sw1_pgap_widen,
                             NULL, NULL, param->alg_type, NULL);

    /*superimposed residues selection*/
    *pbs = pair_set_sel_geti(set_alg);
    set_alg_tmp = selected_pair_set_get(set_alg, *pbs);
    if (set_alg_tmp == NULL) {
        err_printf("pw_algnt", "Can not get any selected pairs from alignment.\n");
        return 0;
    }
    do {
        if (coord_rmsd(set_alg_tmp, coord1, coord2, 0, rmsd_out, &coord1, &coord2) == EXIT_SUCCESS) {
            /*remove max rmsd paar. */
            *pbs = pair_set_sel_delmaxdistance(coord1, coord2, set_alg, *pbs);
            pair_set_destroy(set_alg_tmp);
            set_alg_tmp = selected_pair_set_get(set_alg, *pbs);
        } else {
            err_printf("pw_algnt", "coord_rmsd was not successful. \n");
            break;
        }
    } while ((*rmsd_out > param->rmsd_thresh) && (set_alg_tmp->n > param->mag_num_asrest));
    
    *num = set_alg_tmp->n;
    score_mat_destroy(crap);
    aa_strct_clssfcn_destroy(class);
    score_mat_destroy(matrix);
    prob_vec_destroy(pvec1);
    prob_vec_destroy(pvec2);
    pair_set_destroy(set_alg_tmp);
    mprintf("about leave pw_align pbs is %s\n", *pbs);
    
    return set_alg;
}