.POSIX:

EXAMPL = examples

TEMP    = temp
TEST  = ./tools

LIBS  = src/wurstsrc
LIBD    = lib

INCS  = src/wurstsrc
INCSS = gsldir
INCD  = inc

GSLS  = src/wurstsrc
GSLSS = gsldir
GSLD  = gsl_lib

GSLS1 = $(GSLSS)/matrix.c

GSLS2 = blas.c eval.h gsl_message.h gsl_precision.h init_source.c matrix_init_source.c \
  matrix_view_source.c source_gemv_r.h templates_on.h block_init_source.c \
  gsl_complex.h gsl_nan.h gsl_sf_exp.h lu.c matrix_rowcol_source.c source_trmv_r.h \
  view.h erfc.c gsl_complex_math.h gsl_permute_double.h gsl_sys.h \
  matrix_source.c permutation_init.c source_trsv_r.h view_source.c cheb_eval.c error.c \
  gsl_machine.h gsl_permute_vector_double.h gsl_vector_complex.h matrix_swap_source.c stream.c \
  chebyshev.h error.h gsl_math.h gsl_pow_int.h matrix_view.h \
  source_gemm_r.h templates_off.h

GSLS3 = e_malloc.c mprintf.c
GSLS_ = $(GSLS1) $(GSLS2) $(GSLS3)

INC1 = gsl_blas.h gsl_blas_types.h gsl_block_double.h gsl_cblas.h gsl_check_range.h \
       gsl_errno.h gsl_linalg.h gsl_matrix_double.h gsl_mode.h gsl_permutation.h gsl_sf_erf.h \
       gsl_sf_result.h gsl_types.h gsl_vector_double.h
INC2 = e_malloc.h mprintf.h
INC_ = $(INC1) $(INC2)

LIB1 = align.c align_i.h amino_a.c amino_a.h bad_angle.h binver.c binver.h bit_set.h bit_set.c\
       class_model.c class_model.h classifyStructure.c classifyStructure.h common.c \
       coord.c coord.h coord_i.h dihedral.c dihedral.h fio.c fio.h lsqf.c lsqf.h \
       matrix.c matrix.h mgc_num.c mgc_num.h pair_set.c pair_set.h pair_set_i.h pair_set_p.c \
       pair_set_sel.c pair_set_sel.h\
       pair_set_p_i.h pdbin.c pdbin.h pdbin_i.h pdbout.c pdbout_i.h prob_vec.c  prob_vec.h prob_vec_i.h \
       read_ac.c read_ac.h read_ac_i.h read_ac_strct.c read_ac_strct.h read_ac_strct_i.h \
       read_seq.c read_seq_i.h scor_set.h score_mat.c score_mat.h score_mat_i.h score_probvec.c \
       score_probvec.h scratch.c scratch.h sec_s.c sec_s.h sec_s_i.h seq.h seq_i.h seqprof.h \
       str.c str.h vec.c vec_i.h yesno.h
LIB_ = $(LIB1)