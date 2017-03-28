include ../Make.inc

SRCS =  spline1.f90  

OBJS = $(SRCS:%.c=%.o)
DOBJS = $(SRCS:%.c=%.do)

ifneq ($(MAKECMDGOALS),debug)
XOBJS = $(OBJS)
else
XOBJS = $(DOBJS)
endif

CFLAGS-add += -DMATHLIB_STANDALONE

default: release

%.o: %.f90
	$(FC) $(FFLAGS-add) $(FFLAGS) $(fPIC) -std=gnu99  -DNDEBUG -c $< -o $@
%.do: %.f90
	$(FC) $(FFLAGS-add) $(FFLAGS) $(fPIC) -std=gnu99  -g  -c $< -o $@

release debug: fortran-code.$(SHLIB_EXT)

fortran-code.$(SHLIB_EXT): $(XOBJS)
	rm -rf $@
ifeq ($(USE_DSFMT),1)
	$(FC) $(fPIC) $(LDFLAGS-add) $(LDFLAGS) -shared -o $@ $^ -L$(DSFMT_libdir) -ldSFMT
else
	$(FC) $(fPIC) $(LDFLAGS-add) $(LDFLAGS) -shared -o $@ $^
endif

clean:
	rm -f *.o *.do *.a *.$(SHLIB_EXT) core* *~ *#
