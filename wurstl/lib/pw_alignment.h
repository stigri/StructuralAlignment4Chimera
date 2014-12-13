/*
 * This is an addition file with the function to build a pair-wise alignment fast and simple
 */

#ifndef PW_ALIGNMENT_H
#define PW_ALIGNMENT_H

struct algnm_param {
    float sw1_pgap_open;           /* open gup penalty */
    float sw1_pgap_widen;           /* number of attributes per class */
    float mat_shift;               /* log of attribute probability */
    unsigned char mag_num_asrest;  /* Normalised weight of each class */
    float rmsd_thresh;
    int alg_type;
};


/* initialize the algnm_param to default parameters,
   you can change them afterwards if you need    */
struct algnm_param *
init_algnm_param (void);

/* gets alignment back and transformed coord of the coord1
   and return 0 if fails oer a pointer to an alignment,
   real rmsd in rmsd_out,
   num is the number of atoms which were superimposed
   pbs is a pointer to the bitset of the index of atoms which were superimposed*/
struct pair_set *
pw_algnt(struct coord *coord1, struct coord *coord2, struct algnm_param* param, float *rmsd_out, size_t *num, char **pbs);


#endif