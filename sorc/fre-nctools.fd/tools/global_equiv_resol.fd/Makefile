SHELL := bash

MAKEFLAGS += --warn-undefined-variables

EXEC = global_equiv_resol

.PHONY: all
all : $(EXEC)

$(EXEC): global_equiv_resol.o
	$(FC) $(FFLAGS) $(LIBS) -o $@ $^

.SUFFIXES:
.SUFFIXES: .f90 .o

.f90.o:
	$(FC) $(FFLAGS) -c $<

.PHONY: clean
clean:
	rm -f *.o *.mod $(EXEC)

install:
	cp $(EXEC) ../../../../exec
