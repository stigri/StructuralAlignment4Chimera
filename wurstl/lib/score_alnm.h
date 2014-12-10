/*
 * 11 Juni 2014
 * This defines the structures for annotating an scores of alignment, especially 
 * for pairwise alignments like in Yana project, wurstl tool or in ost project
 */
#ifndef SCORE_ALNM_H
#define SCORE_ALNM_H

struct score_struct{
    float scr_tot;
    float cvr;
};

struct score_struct*
get_scores (struct pair_set *set, struct coord *a, struct coord *b, const float sw1_pgap_open,
             const float sw1_pgap_widen);

#endif