.POSIX:

# 22 March 2011
# To build a library for the main executable
# The dependencies on the include files are missing. They should be put
# in.

OUTLIB=libgslw.a
OUTLIBSTAT=libgslw.a
OUTLIBDYN=libgslw.so

INC = ../inc
MODE=static

CFILES = blas.c\
         block_init_source.c\
         cheb_eval.c\
         e_malloc.c\
         erfc.c\
         error.c\
         init_source.c\
         matrix_init_source.c\
         matrix_rowcol_source.c\
         matrix_swap_source.c\
         matrix_view_source.c\
         matrix.c\
         mprintf.c\
         permutation_init.c\
         stream.c\
         view_source.c\
         lu.c

OFILES = $(CFILES:.c=.o) 

all: msgs $(MODE) 
	@echo "gsl_lib buld success!"
   
msgs: 
	@echo "Start to build in $(MODE) mode..."
   

dynam: CFLAGS  = $(CFLAGSD)   
dynam: LDFLAGS = $(LDFLAGSD)
dynam: $(OUTLIBDYN)_dyn

static: CFLAGS  = $(CFLAGSS) 
static: LDFLAGS = $(LDFLAGSS)
static: $(OUTLIBSTAT)_stat
   
#dynamic lib   
$(OUTLIBDYN)_dyn:$(OFILES)
	$(CC) $(CFLAGS) -Wl,-soname,$(OUTLIBDYN).1 -o $(OUTLIBDYN).1  $(CFILES:.c=.o) $(LDFLAGS) 
# 	ln -sf $(OUTLIBDYN) $(OUTLIBDYN).0
	ln -sf $(OUTLIBDYN).1 $(OUTLIBDYN)
	
#static lib   
$(OUTLIBSTAT)_stat:$(CFILES:.c=.o) 
	ar cvsq $(OUTLIBSTAT) $(CFILES:.c=.o)
	@echo OK
	
CFLAGSS  = -I $(INC) -Dinline=  -fPIC -g
LDFLAGSS = -lm 
#CFLAGSD  = -I $(INC) -Dinline= -DHIDE_INLINE_STATIC -fPIC -g
CFLAGSD  = -I$(INC) -Dinline=  -fPIC -g
LDFLAGSD = -lm -shared

clean:
	-rm -f *.o *.a *.so* *~
	#rm -f *.a
