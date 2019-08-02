SHELL := bash

REGIONAL_GRID=regional_grid

.PHONY: all
all : $(REGIONAL_GRID)

$(REGIONAL_GRID): pkind.o pietc.o pmat.o pmat4.o pmat5.o psym2.o gen_schmidt.o hgrid_ak.o regional_grid.o
	$(FCMP) $(FFLAGS) ${LIBS} -o $@ $^

.SUFFIXES:
.SUFFIXES: .f90 .o

.f90.o:
	$(FCMP) $(FFLAGS) -c $<

.PHONY: clean
clean:
	rm -f *.o *.mod $(REGIONAL_GRID)

install:
	cp $(REGIONAL_GRID) ../../../../exec
