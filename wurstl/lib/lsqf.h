/* "$Id: lsqf.h,v 1.4 2011/09/07 14:36:03 margraf Exp $" */
#ifndef LSQF_H 
#define LSQF_H

struct pair_set;

void
coord_trans(const struct RPoint *trans, struct coord *structure,
            const size_t size);

struct RPoint
calc_CM(const int nr_atoms, const struct RPoint *r1);

int
lsq_fit(const int nr_atoms, const struct RPoint *r1, const struct RPoint *r2,
        float R[3][3]);

void
copy_coord_elem(struct coord *dst, const struct coord *src, const int dst_ndx,
                const int src_ndx);

int get_rmsd(struct pair_set *pairset, struct coord *r1,
             struct coord *r2, float *rmsd_ptr, int *count);

int coord_rmsd (struct pair_set *pairset, struct coord *coord1,
                struct coord *coord2, int sub_flag, float *rmsd,
                struct coord **c1_new, struct coord **c2_new);

double
tm_score (struct pair_set *pairset, struct coord *coord1, struct coord *coord2);

double
tm_score_s (struct pair_set *pairset, struct coord *coord1, struct coord *coord2);

void coord_centre (struct coord *c);

#endif
