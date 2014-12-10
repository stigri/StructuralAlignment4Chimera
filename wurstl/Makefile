.POSIX:
# For building wurstl which simply does pairwise alignments of proteins.


# These flags will be passed down to the subdirectories
CDEBUG = -O2
# gcc
#CWARN  = -ansi -pedantic -Wall
CWARN  =  -pedantic -Wall
# suns cc
#CWARN  = -Xc -xc99=none
CMISC  = -D_XOPEN_SOURCE=600
CMISC  = -D_POSIX_SOURCE -D_XOPEN_SOURCE

COPT = $(CDEBUG) $(CWARN) $(CMISC)

REP = /smallfiles/bm/cvs
#BRANCH = wurstl_branch
#BRANCH = release120425
#BRANCH = fixing_wurstl
#BRANCH = rev_20131220
BRANCH = rev_20140617
#BRANCH_LIB = wurstl_lib_branch1
#BRANCH_LIB = fixing_wurstl
#BRANCH_LIB = rev_20131220
PRJROOT = wurst
DIRS = gsl_lib, lib, tools

include files.mk
# LIB_ = $(LIB1)

MODE=static


all: build

help:
	@echo use: 'make checkout or co' or\
	'make build' or\
	'make cob' for checkout und build

cob: checkout build

checkout: gsl $(GSLS_) inc $(INC_) lib $(LIB_)
co : checkout

#checkin, ci: 
#	echo it's not done yet.
build:
	@echo "I'm about to start wurstl in $(MODE) mode....."
	cd $(TEST) && $(MAKE) MODE='$(MODE)' COPT='$(COPT)' all

gsl: 
	mkdir -p $(GSLD)

$(GSLS1):
	cvs -d $(REP) co -d $(TEMP)/$(GSLD)/1 $(PRJROOT)/$(GSLS)/$@
	cvs -d $(REP) update -r $(BRANCH) $(TEMP)/$(GSLD)/1/$(@:$(GSLSS)/%=%)
	ln -s -f ../$(TEMP)/$(GSLD)/1/$(@:$(GSLSS)/%=%) ./$(GSLD)/$(@:$(GSLSS)/%=%)
# 	ln -s -f ../$(TEMP)/$(GSLD)/1/$(@F) ./$(GSLD)/$(@F)

$(GSLS2):
	cvs -d $(REP) co -d $(TEMP)/$(GSLD)/2 $(PRJROOT)/$(GSLS)/$(GSLSS)/$@
	cvs -d $(REP) update -r $(BRANCH) $(TEMP)/$(GSLD)/2/$@
	ln -s -f ../$(TEMP)/$(GSLD)/2/$@ ./$(GSLD)/$@

$(GSLS3): 
	cvs -d $(REP) co -d $(TEMP)/$(GSLD)/3 $(PRJROOT)/$(GSLS)/$@
	cvs -d $(REP) update -r $(BRANCH) $(TEMP)/$(GSLD)/3/$@
	ln -s -f ../$(TEMP)/$(GSLD)/3/$@ ./$(GSLD)/$@

inc:
	mkdir -p $(INCD)

$(INC1):
	cvs -d $(REP) co -d $(TEMP)/$(INCD)/1 $(PRJROOT)/$(INCS)/$(INCSS)/$@
	cvs -d $(REP) update -r $(BRANCH) $(TEMP)/$(INCD)/1/$@
	ln -s -f ../$(TEMP)/$(INCD)/1/$@ ./$(INCD)/$@

$(INC2):
	cvs -d $(REP) co -d $(TEMP)/$(INCD)/2 $(PRJROOT)/$(INCS)/$@
	cvs -d $(REP) update -r $(BRANCH) $(TEMP)/$(INCD)/2/$@
	ln -s -f ../$(TEMP)/$(INCD)/2/$@ ./$(INCD)/$@

lib:
	mkdir -p $(LIBD)

$(LIB1):
	cvs -d $(REP) co -d $(TEMP)/$(LIBD)/1 $(PRJROOT)/$(LIBS)/$@
	cvs -d $(REP) update -r $(BRANCH) $(TEMP)/$(LIBD)/1/$@  #upadtes back into the branch
	ln -s -f ../$(TEMP)/$(LIBD)/1/$@ ./$(LIBD)/$@

clean: 
	 @echo Cleaning...
	-rm *.o *.a *~
	for d in $(DIRS) ; do  \
             cd $$d && $(MAKE) clean ; \
        done
