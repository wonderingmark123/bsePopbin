CMPLR = gfortran
# FFLAGS = -O3 -g -fbacktrace -fcheck=all -ffpe-trap=zero,overflow -finit-real=nan

FFLAGS = -O3 -g -fbacktrace -fcheck=all -finit-real=nan 
# haotian  
#  the setting above is explained bt following website. 
# https://gcc.gnu.org/onlinedocs/gfortran/Debugging-Options.html
# remove the -ffpe-trap=zero,overflow  to avoid stopping progress.
LFLAGS = const_bse.h zdata.h

all:	sse bse popbin betest

.f.o:
	$(CMPLR) -c $(FFLAGS) $<

be.o:	be.f be.h $(shell find bedata -type f)

SRCE1 = \
sse.f deltat.f evolv1.f hrdiag.f kick.f mlwind.f mrenv.f \
ran3.f star.f zcnsts.f zfuncs.f

OBJT1 = $(SRCE1:.f=.o) be.o

sse:    $(OBJT1) $(LFLAGS)
	$(CMPLR) $(FFLAGS) $(OBJT1) -o sse

SRCE2 = \
bse.f comenv.f corerd.f deltat.f dgcore.f evolv2.f gntage.f \
hrdiag.f instar.f kick.f mix.f mlwind.f mrenv.f ran3.f rl.f \
star.f zcnsts.f zfuncs.f myfun.f

OBJT2 = $(SRCE2:.f=.o) be.o

bse:    $(OBJT2) $(LFLAGS)
	$(CMPLR) $(FFLAGS) $(OBJT2) -o bse

SRCE3 = \
popbin_mine.f comenv.f corerd.f deltat.f dgcore.f evolv2.f gntage.f \
hrdiag.f instar.f kick.f mix.f mlwind.f mrenv.f ran3.f rl.f \
star.f zcnsts.f zfuncs.f myfun.f

OBJT3 = $(SRCE3:.f=.o) be.o

popbin: $(OBJT3) $(LFLAGS)
	$(CMPLR) $(FFLAGS) $(OBJT3) -o popbin

betest.o: betest.f90
	$(CMPLR) -c $(FFLAGS) betest.f90

betest: betest.o be.o
	$(CMPLR) $(FFLAGS) betest.o be.o -o betest

.PHONY: clean
clean:
	-rm -f sse bse popbin betest *.o
