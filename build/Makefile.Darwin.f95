#^CFG COPYRIGHT UM

SHELL=/bin/sh

# Fortran language related part of Makefile.conf: Makefile.Darwin.f95
#
#	Space Weather Modeling Framework (SWMF) 
#	NAG f95 Fortran 90/95 Compiler
#

COMPILE.f77     = ${CUSTOMPATH_F}f95
COMPILE.f90     = ${CUSTOMPATH_F}f95
LINK.f90	= ${CUSTOMPATH_MPI}mpif90
AR = ar -rs

SINGLEPREC =
DOUBLEPREC = -r8
PRECISION  = ${DOUBLEPREC}

MPILIB = 
#MPILIB = -L${LIBDIR} -lNOMPI

# This is the search path for used modules
# SEARCH_EXTRA should be set in the individual Makefiles

SEARCH = -I${SHAREDIR} ${SEARCH_EXTRA}

DEBUGFLAG = -C -gline -nan
DEBUG     = 

OPT0 = -O0
OPT1 = -O1
OPT2 = -O2
OPT3 = -O3
OPT4 = -O4

CFLAG = ${SEARCH} -c -w ${DEBUG}

Cflag0  = ${CFLAG} ${PRECISION} ${OPT0}
Cflag1  = ${CFLAG} ${PRECISION} ${OPT1}
Cflag2  = ${CFLAG} ${PRECISION} ${OPT2}
Cflag3  = ${CFLAG} ${PRECISION} ${OPT3}
Cflag4  = ${CFLAG} ${PRECISION} ${OPT4}

# RCM compilation flags
# To allow RCM to compile as double precision, add PRECISION flag
CFLAGS = ${CFLAG} -save

Lflag1  = ${PRECISION} ${MPILIB} ${DEBUG}
Lflag2  = ${PRECISION} ${DEBUG}

# BLAS and LAPACK libraries
LBLAS =
BLAS  = lapack.o blas.o


#
#       General rules
#

.SUFFIXES:
.SUFFIXES: .f90 .F90 .f .for .ftn .o

.f90.o:
	${COMPILE.f90} ${Cflag4} $<

.F90.o:
	cpp -P -DsysDarwin -DcompNAGF95 $*.F90 | grep -v '^#' > $*.f95
	${COMPILE.f90} -o $*.o ${Cflag4} $*.f95
	rm -f $*.f95

.f.o:
	${COMPILE.f77} ${Cflag2} -132 $<

.for.o:
	${COMPILE.f77} ${Cflag2} -132 $<

.ftn.o:
	${COMPILE.f77} ${Cflag2} -132 $<

clean:	
	rm -f *~ core *.o *.mod fort.* a.out *.exe *.a *.so *.protex


# keep this line
