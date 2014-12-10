.POSIX:

# 22 March 2011
# To build a library for the main executable
# The dependencies on the include files are missing. They should be put
# in.

OUTLIB=libwurstl.a
OUTLIBSTAT=libwurstl.a
OUTLIBDYN=libwurstl.so

INC = ../inc
GSL = /local/bondarenko/prg/wurstl/gsl_lib
GSLLIB = gslw
MODE=static

CFILES = align.c\
         amino_a.c\
         bit_set.c\
         binver.c\
         class_model.c\
         classifyStructure.c\
         common.c\
         coord.c\
         dihedral.c\
         fio.c\
         lsqf.c\
         matrix.c\
         mgc_num.c\
         pair_set.c\
         pair_set_p.c\
         pair_set_sel.c\
         pdbin.c\
         pdbout.c\
         prob_vec.c\
         read_ac.c\
         read_ac_strct.c\
         read_seq.c\
         score_mat.c\
         score_alnm.c\
         score_probvec.c\
         scratch.c\
         sec_s.c\
         str.c\
         vec.c\
         pwrapper.c
         
OFILES = $(CFILES:.c=.o)      

#$(OUTLIB) :   
all: $(MODE)
	@echo "wurstl has been built in $(MODE) mode."

dynam: CFLAGS  = $(CFLAGSD)
dynam: LDFLAGS = $(LDFLAGSD)
dynam: $(OUTLIBDYN)

static: CFLAGS  = $(CFLAGSS)
static: LDFLAGS = $(LDFLAGSS)
static: $(OUTLIBSTAT)
   
#dynamic lib   
$(OUTLIBDYN):$(OFILES)
	$(CC) $(CFLAGS) -Wl,-rpath,$(GSLLIB) -o $@ $(OFILES) $(LDFLAGS)
	ln -sf $(OUTLIBDYN) $(OUTLIBDYN).1
	
#static lib   
$(OUTLIBSTAT):$(OFILES) 
	ar csvq $(OUTLIBSTAT) $(OFILES)
	@echo OK
	
#CFLAGSS  = -I $(INC) -Dinline= -DHIDE_INLINE_STATIC -fPIC -g
CFLAGSS  = -I$(INC) -I/usr/include/python2.7 -Dinline=  -fPIC -g
LDFLAGSS =  -l$(GSLLIB) -L$(GSL) -l/usr/include/python2.7 -lm
CFLAGSD  = -I$(INC) -I/usr/include/python2.7 -Dinline=  -fPIC -g
LDFLAGSD =  -shared -l$(GSLLIB) -L$(GSL) -L/usr/include/python2.7 -lpython2.7 -lm


clean:
	-rm -f *.o  *~
	-rm -f *.a *.so*
