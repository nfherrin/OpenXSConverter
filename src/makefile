#-----------------------------------------------------------------------------------
#Compile Options
#-----------------------------------------------------------------------------------
FC := gfortran
EXE := OpenXSConverter
EXT := .exe

all:			TYPE :=
MPparallel:    TYPE := _MP
debug:    TYPE := _Debug
catch:    TYPE := _Catch
coverage: TYPE := _Coverage

all:      OPT := -O3
MPparallel: OPT := -O3 -fopenmp
debug:    OPT := -O0 -g -fbacktrace -ffpe-trap=zero,overflow,underflow -fcheck=all
catch:    OPT := -O0 -g -Wall -W -Wsurprising -Werror
coverage: OPT := -O0 -g --coverage

#-----------------------------------------------------------------------------------
#Directory Paths
#-----------------------------------------------------------------------------------
CONTRIB := ../../../contrib

#-----------------------------------------------------------------------------------
#Add source files as necessary
#If the files are not compiled using the generic rules, add commands under their
#rule declaration. Add these items to FILTER
#-----------------------------------------------------------------------------------
SRC := \
				main.f90 \
				globals.f90 \
        precisions.f90 \
				outfuncs.f90 \
				string_module.f90 \
				infuncs.f90

OBJ := $(SRC:.f90=.o)
MOD := $(OBJ:.o=.mod)

FILTER := \

OBJ_FILTER := $(FILTER:.f90=.o)
MOD_FILTER := $(FILTER:.f90=.mod)

#h5 datafiles inclusion
LIB=usr/lib/x86_64-linux-gnu/hdf5/serial
INC=usr/include/hdf5/serial

#-----------------------------------------------------------------------------------
#Phony targets for cleaning and building
#-----------------------------------------------------------------------------------
.PHONY: all clean reset

all: $(EXE)

MPparallel: $(EXE)

debug: $(EXE)

catch: $(EXE)

coverage: $(EXE)

#Intended to clean up compilation artifacts but leave executable & coverage
clean:
	rm -f $(OBJ) $(CONT_OBJ)
	rm -f $(MOD) $(CONT_MOD)
	rm -f $(COMP_DEP)

#Intended to reset directory to fresh state with no exe or artifacts
reset: clean
	rm -f *.gcno *.gcda
	rm -f -r $(EXE)*.dSYM
	rm -f ../$(EXE)*

#-----------------------------------------------------------------------------------
#Generics for source files
#-----------------------------------------------------------------------------------
$(filter-out $(OBJ_FILTER), $(OBJ)): %.o:	%.f90
	$(FC) -c $(OPT) $< -I/$(INC) -L/$(LIB) -lhdf5_fortran -lhdf5

$(filter-out $(MOD_FILTER), $(MOD)):	%.mod:	%.f90
	$(FC) -c $(OPT) $< -I/$(INC) -L/$(LIB) -lhdf5_fortran -lhdf5

$(EXE): $(OBJ)
	$(FC) -o $@$(TYPE)$(EXT) $(OPT) $(OBJ) -I/$(INC) -L/$(LIB) -lhdf5_fortran -lhdf5
	mv ./$(EXE)$(TYPE)$(EXT) ../

#-----------------------------------------------------------------------------------
#Rules for entries in COMP_DEP. Don't forget to add them to make clean / reset
#-----------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------
#Dependency List
#Use [gfortran -M -cpp *.f90] repeatedly until clean compile to update rules below
#-----------------------------------------------------------------------------------
globals.o globals.mod: globals.f90
infuncs.o infuncs.mod: infuncs.f90 globals.mod string_module.mod
main.o: main.f90 globals.mod infuncs.mod outfuncs.mod
outfuncs.o outfuncs.mod: outfuncs.f90 globals.mod
precisions.o precisions.mod: precisions.f90
string_module.o string_module.mod: string_module.f90 precisions.mod