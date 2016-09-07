.SILENT:
SHELL = /bin/sh
NAME=./bin/negf-pt


DFT=$(NAME)

F77     = gfortran

#FFLAGS  = -check bounds -debug all -module objmod -traceback -ftrapuv -O0 -fpe0
#FFLAGS  = -O0 -fno-backtrace -fbounds-check  -J./objmod 
#O0 is the optimization level 0(without optimization)     for debug

FFLAGS  = -O2 -fintrinsic-modules-path objmod -J./objmod

#FOPTS   = -static

# for huanghe
# #LIBDIR  =  /share/apps/intel/Compiler/11.1/073/mkl/lib/em64t
# # for alps
# #LIBDIR  = /share/apps/intel/composerxe-2011.2.137/mkl/lib/intel64
# #LIBS    =  -Wl,--start-group $(LIBDIR)/libmkl_intel_lp64.a $(LIBDIR)/libmkl_intel_thread.a $(LIBDIR)/libmkl_core.a -Wl,--end-group -openmp -lpthread
#
LIBDIR  = /usr/lib64
LIBS    =  $(LIBDIR)/liblapack.so.3 $(LIBDIR)/libblas.so.3




all:;
	make lodestar

lodestar:; 
	make DIR=module        obj    SHELL=/bin/sh  FF=$(F77)
	make DIR=main          obj    SHELL=/bin/sh  FF=$(F77)
	make DIR=rmatrix       obj    SHELL=/bin/sh  FF=$(F77)
	make $(DFT)


$(DFT): objects/*.o 
	rm -f $(DFT)
	$(F77) $(FOPTS)  objects/*.o $(LIBS) -o $(DFT)
	echo Compiling --$(DFT)-- done !

obj:;
	ext=".f90";\
	for d in  $$DIR/*.f90; \
	do \
	  filename=`basename $$d $$ext`; \
	  make FILE=$$filename DIR=$$DIR FF=$$FF objects/$$filename.o; \
	done 

objects/$(FILE).o: $(DIR)/$(FILE).f90 ;
	$(FF) $(FFLAGS) -c $(DIR)/$(FILE).f90;\
	mv $(FILE).o objects;\
	echo $(FF) $(FFLAGS) -c $(DIR)/$(FILE).f90
	echo update $(DIR)/$(FILE)

clean:
	rm -f $(DFT) objects/*.o objmod/*.mod
	echo delete objects/*.o objmod/*.mod $(DFT)

tags: */*.f90 
	@ctags  */*.f90 
	@echo tags updated!
