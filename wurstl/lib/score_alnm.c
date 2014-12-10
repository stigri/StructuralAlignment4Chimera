/*
 * 11 Juni 2014
 * This defines the structures for annotating an scores of alignment, especially 
 * for pairwise alignments like in Yana project, wurstl tool or in ost project
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "coord.h"
#include "coord_i.h"
#include "align_i.h"
#include "prob_vec_i.h"
#include "e_malloc.h"
#include "pair_set.h"
#include "score_alnm.h"

struct score_struct*
get_scores (struct pair_set *set, struct coord *a, struct coord *b,
            const float sw1_pgap_open, const float sw1_pgap_widen)
{
    int cover=0,i ;
    float cover_float;
    char *pcover1;
    char *pcover2;
    struct seq* seq_a;
    struct score_struct *scores;
    float open_cost = 0, widen_cost = 0;

    scores = E_MALLOC(sizeof (*scores));

    pair_set_coverage(set, coord_size(a), coord_size(b),
                            &pcover1, &pcover2);
    for (i=0; i < coord_size(a); i++) {
        if(pcover1[i] == '1') cover++;
    }

    free(pcover1);
    free(pcover2);

    seq_a = coord_get_seq(a);
    cover_float = (float)cover/(float)seq_size(seq_a);

    if (cover_float >= 0.05) {
        pair_set_gap(set, &open_cost, &widen_cost, 1, 1);
    }

    scores->scr_tot = set->smpl_score + sw1_pgap_open*open_cost + sw1_pgap_widen*widen_cost;

    seq_destroy(seq_a);
    scores->cvr = cover_float;

    return scores;
}