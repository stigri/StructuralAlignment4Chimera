/*
 * 10 Jan 2002
 * Just calculation of a dihedral angle.
 *
 * $Id: dihedral.c,v 1.2 2008/06/27 08:58:07 schenk Exp $
 */

#include <math.h>
#include <stdlib.h>

#include "coord.h"
#include "dihedral.h"
#include "vec_i.h"

/* ---------------- dihedral ----------------------------------
 * A helper function for the calculation of the dihedral angles
 * (phi and psi) in a residue. This returns the included torsion
 * angle with the correct sign.
 * The functions are taken from van Gunsteren et al, "Biomolecular
 * Simulation: The GROMOS96 Manual and User Guide", page II-26.
 */
float
dihedral(const struct RPoint ii, const struct RPoint jj,
         const struct RPoint kk, const struct RPoint ll)
{
    const struct RPoint *i = &ii,
                        *j = &jj,
                        *k = &kk,
                        *l = &ll;
    float tmp, theta;
    struct RPoint r_im, r_ln, tmp_vec;
    struct RPoint r_ij, r_kj, r_kl;

    vector_difference(&r_ij, i, j);
    vector_difference(&r_kj, k, j);
    vector_difference(&r_kl, k, l);

    tmp = scalar_product(&r_ij, &r_kj);
    tmp = tmp / vector_sqr_length(&r_kj);
    tmp_vec.x = tmp * r_kj.x;
    tmp_vec.y = tmp * r_kj.y;
    tmp_vec.z = tmp * r_kj.z;
    vector_difference(&r_im, &r_ij, &tmp_vec);

    tmp = scalar_product(&r_kl, &r_kj);
    tmp = tmp / vector_sqr_length(&r_kj);
    tmp_vec.x = tmp * r_kj.x;
    tmp_vec.y = tmp * r_kj.y;
    tmp_vec.z = tmp * r_kj.z;
    vector_difference(&r_ln, &tmp_vec, &r_kl);

    tmp = scalar_product(&r_im, &r_ln);
    theta = tmp / (vector_length(&r_im) * vector_length(&r_ln));
    if (theta > 1)
        theta = 1.0;  /* Numerical errors can catch us here */
    else if (theta < -1)
        theta = -1.0;

    theta = acos(theta);

    vector_product(&tmp_vec, &r_kj, &r_kl);
    tmp = scalar_product(&r_ij, &tmp_vec);
    if (tmp >= 0)
        return (theta);
    else
        return (- theta);
}

/* ---------------- bondangle ----------------------------------
 * A helper function for the calculation of the bond angles
 * in a residue. This returns the included 
 * angle with the correct sign.
 */
float
bondangle(const struct RPoint ii, const struct RPoint jj,
          const struct RPoint kk)
{
    const struct RPoint *i = &ii,
                        *j = &jj,
                        *k = &kk;
    float tmp, theta;
    struct RPoint tmp_vec;
    struct RPoint r_ji, r_jk;

    vector_difference(&r_ji, j, i);
    vector_difference(&r_jk, j, k);

    tmp = scalar_product(&r_ji, &r_jk);

    tmp /= vector_length(&r_ji) * vector_length(&r_jk);
    
    if (tmp > 1)
        tmp = 1.0;  /* Numerical errors can catch us here */
    else if (tmp < -1)
        tmp = -1.0;

    theta = acos(tmp);

    vector_product(&tmp_vec, &r_ji, &r_jk);
    tmp = vector_length(&tmp_vec) / (vector_length(&r_ji) 
                                    * vector_length(&r_jk));

    if (tmp > 1)
        tmp = 1.0;  /* Numerical errors can catch us here */
    else if (tmp < -1)
        tmp = -1.0;

    tmp = asin(tmp);

    if (tmp >= 0)
        return (theta);
    else
        return (- theta);
}

