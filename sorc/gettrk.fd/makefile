SHELL=  /bin/sh
LDFLAGS= 

gettrk:      gettrk_main.f gettrk_modules.o module_waitfor.o cwaitfor.o
	@echo " "
	@echo "  Compiling the main tracking program and subroutines....."
	$(FCOMP) $(FFLAGS) $(LDFLAGS) gettrk_modules.o module_waitfor.o cwaitfor.o gettrk_main.f $(LIBS) -o gettrk
	@echo " "

cwaitfor.o: cwaitfor.c
	@echo " "
	@echo "  Compiling the waitfor C routine...."
	$(CCOMP) $(CFLAGS) -c cwaitfor.c -o cwaitfor.o

module_waitfor.o: module_waitfor.f
	@echo " "
	@echo "  Compiling the waitfor fortran module...."
	$(FCOMP) $(FFLAGS) -c module_waitfor.f -o module_waitfor.o

gettrk_modules.o:    gettrk_modules.f
	@echo " "
	@echo "  Compiling the regular tracker fortran modules....."
	$(FCOMP) $(FFLAGS) -c gettrk_modules.f -o gettrk_modules.o
	@echo " "

CMD =   gettrk

clean:
	-rm -f  *.o  *.mod $(CMD)

install:
	mv $(CMD) ../../exec/$(CMD)
