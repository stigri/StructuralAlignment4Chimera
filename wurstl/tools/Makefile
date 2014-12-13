.POSIX:
# This is for the main executable which will be called wurstl, or
# whatever is given by EXE_NAME

EXE_NAME1  = struct_pwa
EXE_NAME2  = struct_pwa_simple
LIB	  = ../lib
INC   = ../inc
GSL   = ../gsl_lib
IRGX = ../mingw-libgnurx-2.5.1-dev/include
LRGX = ./
MODE=static

CFLAGS  = $(COPT) -I$(LIB) -I$(GSL) -I$(INC) -I$(IRGX)

LDFLAGS = -L$(GSL) -L$(LIB) -lwurstl -lgslw -lm


all: $(MODE)

libgsl.a:
	cd $(GSL) && $(MAKE) MODE='$(MODE)' COPT='$(COPT)' CC=$(CC) all

libwurstl.a:
	cd $(LIB) && $(MAKE) MODE='$(MODE)' COPT='$(COPT)' CC=$(CC) all

static:$(EXE_NAME1) $(EXE_NAME2)

dynam: $(EXE_NAME1)_dyn $(EXE_NAME2)_dyn

$(EXE_NAME1)_dyn: libgsl.a libwurstl.a main.o
	$(CC) $(CFLAGS) -Wl,-rpath,$(GSL) -Wl,-rpath,$(LIB) -o $(EXE_NAME1) main.o $(LDFLAGS)
	
$(EXE_NAME1): libgsl.a libwurstl.a main.o
	$(CC) $(CFLAGS) -o $@  main.o $(LDFLAGS) 

$(EXE_NAME2)_dyn: libgsl.a libwurstl.a simple.o
	$(CC) $(CFLAGS) -Wl,-rpath,$(GSL) -Wl,-rpath,$(LIB) -o $(EXE_NAME2) simple.o $(LDFLAGS)
	
$(EXE_NAME2): libgsl.a libwurstl.a simple.o
	$(CC) $(CFLAGS) -o $@  simple.o $(LDFLAGS) 
	
clean: 
	 @echo Cleaning...
	rm -f *.o *~
	rm -f $(EXE_NAME1)*
	cd $(LIB) && make clean
	cd $(GSL) && make clean
