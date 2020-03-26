############################################################
#                     MAKEFILE FOR F95                     #
############################################################

# Zero axiom: the file name of a makefile must be mandatorily
#
#              makefile    OR    Makefile

# To start the compilation, type in terminal: make

# The names for aliases are completely free!

# FUNDAMENTAL WARNINGS
#
# The actions $(CMP) and $(LNK) must be NECESSARILY
# preceded by a DAMNED! TAB character
#
# The DAMNED! TAB character is needed also in continuation
# line of actions

############################################################
#
# Definition of VARIABLES (usually in CAPITALS)
#
# Compiler and linker
#
# The following definitions are for Fortran 95 only, but they
# are completely customizable if another compiler is used.

# LNK alias to be used as $(LNK) to compile and link main files .f90
# CMP alias to be used as $(CMP) to compile the source files .f90
# OPT alias to be used as $(OPT) to select the compilation options

# NAG Compiler
# ============
#
# CMP = f95 -c
# LNK = f95
# OPT = -kind=byte ### -kind=byte
# OPT = -gline ### -kind=byte  # Line number is indicated for errors
                               # occurring in runtime (for debugging)
# OPT = -O3    ### -kind=byte  # Execution optimization


# Intel Compiler
# ==============
#
# CMP = ifort  -c  -heap-arrays 100
# LNK = ifort
# OPT =
# OPT = -g  -DD  -debug -traceback  -fpe0  -check all
# OPT = -fast
# OPT = -O3


# GNU fortran Compiler
# ===================
#
CMP = gfortran  -c
LNK = gfortran
#OPT = -g -fcheck=all -Wall -fno-automatic  -fbounds-check  -mcmodel=medium -ggdb3
#OPT = -g -fbacktrace
#OPT = -ggdb -g -O0 -pedantic -Wall -Wconversion -fcheck=all -std=f2008ts -fall-intrinsics -finit-real=snan -ffpe-trap=invalid -finit-integer=-9999 -Wintrinsics-std -fbacktrace
OPT = -O3 -fopenmp
#OPT =

 #OPT = -fno-automatic  -O3

############################################################
#
# Objects: list of all objects *.o

OBJS = mod_initialize.o\
			 mod_types.o\
			 mod_further_beyond.o\
			 mod_geometry.o\
			 mod_aero_coeff.o\
			 mod_velocity.o\
			 mod_check.o\
			 mod_advectVortices.o\
          mod_flapping.o\
			 mod_output.o



# OBJS alias to be used as $(OBJS)


############################################################
# Executable generation by the linker / the effect of command
#                                     / make starts from here

vortex_particle.exe:  vortex_particle.o  $(OBJS)
	$(LNK) $(OPT)   vortex_particle.o  $(OBJS)  \
               -o vortex_particle.exe  # Name chosen for program


	export OMP_NUM_THREADS=4

################# PER PARALLELIZZAZIONEEEEEEEEEEEEE!!!!!!!!!!!!!!!!!!!!! #############################
### Sostituire OMP_NUM_THREADS=8 con il numero di threads del pc. Se non parallelizza allora cancella
### quella riga dal makefile ed eseguila nella bash di ubuntu prima di compilare nuovamente. Poi compila

### NOMORE #  ^
### NOMORE #  |
### NOMORE # A DAMNED! TAB character is leading the previous line
### NOMORE # (see WARNINGS)

############################################################
# Objects generation

vortex_particle.o:       vortex_particle.f90  $(OBJS)
	    $(CMP) $(OPT)  vortex_particle.f90

mod_types.o:	mod_types.f90 mod_initialize.f90
		 $(CMP) $(OPT) mod_types.f90

mod_check.o:	mod_check.f90 mod_geometry.f90
		 $(CMP) $(OPT) mod_check.f90

mod_initialize.o:   mod_initialize.f90
		$(CMP) $(OPT)  mod_initialize.f90

mod_further_beyond.o:   mod_further_beyond.f90
		$(CMP) $(OPT)  mod_further_beyond.f90

mod_geometry.o: mod_geometry.f90 mod_initialize.f90
		$(CMP) $(OPT) mod_geometry.f90

mod_aero_coeff.o: mod_aero_coeff.f90
		$(CMP) $(OPT) mod_aero_coeff.f90

mod_velocity.o: mod_velocity.f90 mod_types.f90 mod_geometry.f90 mod_aero_coeff.f90  mod_types.f90
		$(CMP) $(OPT) mod_velocity.f90

mod_advectVortices.o: mod_advectVortices.f90 mod_types.f90
		$(CMP) $(OPT) mod_advectVortices.f90

mod_output.o: mod_output.f90 mod_initialize.f90 mod_geometry.f90
		$(CMP) $(OPT) mod_output.f90

mod_flapping.o: mod_flapping.f90 mod_initialize.f90
		$(CMP) $(OPT) mod_flapping.f90





############################################################
# Cleaning command to remove all objects *.o, files *.mod and
# the executables *.exe

clean:
	@echo Cleaning objects, modules and executables
	rm   *.o *.mod *.exe

cleandata:
	@echo cleaning video files
	rm  solution/blade_coord/*.dat
	rm  solution/wake_panel_coord/*.dat
	rm  solution/particle_coord/*.dat
	rm  solution/CT.dat

# The command "make clean" deletes all the indicated files

############################################################

cleanmod:
	@echo Cleaning objects, modules
	rm build/*.o build/*.mod
