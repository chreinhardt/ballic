#########################
#                       #
#   Ballic Makefile.    #
#                       #
#########################

CC              = gcc
FC              = gfortran

OBJECTS         = ballic.o ballic.single.o ballic.multi.o modelsolve.o
TIPSY_OBJECTS   = tipsy.o 
TILL_OBJECTS    = tillotson/tillotson.o tillotson/tillinitlookup.o tillotson/tillsplint.o tillotson/interpol/brent.o
FORTRAN_OBJECTS = icosahedron.o
INC_TIPSYDEFS   = -I./tipsydefs
TIPSYSRC        = ./tipsydefs/tipsy.c
BALLIC          = ballic

EXE             = ballic ballic.single ballic.multi modelsolve

# GNU Science library (uncomment if not needed)
GSL_LIB          = -lgsl -lgslcblas
# tirpc library (needed if glibc >= 2.32)
RPC_LIB          = -ltirpc
LIBS            ?= -lm $(GSL_LIB) $(RPC_LIB)

CFLAGS          ?= -O3 -march=native
FFLAGS          ?= $(CFLAGS)

# NOTE:
# Call --> make ballic.single CFLAGS+=-DDOUBLE_TIPSY <-- from command line for double precision tipsy files
# Before switching precision you have to invoke a make cleanall to clean the tipsy.o file

default:
	@echo "Please specify, which tool to make (e.g., ballic, single, multi)."

all:
	default

# Standard version of ballic that reads an equilibrium model from a file and generates an particle distribution,
# that matches the desired density profile. Currently not reliably working for multi component models.

ballic: ballic.o $(TIPSY_OBJECTS) $(FORTRAN_OBJECTS)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

# Single component version of ballic. It first calculates an equilibrium model for a given material and then
# generates a particle representation.

single: ballic.single.o $(TIPSY_OBJECTS) $(TILL_OBJECTS) $(FORTRAN_OBJECTS)
	$(CC) $(CFLAGS) $^ -o $(BALLIC).$@ $(LIBS)

# Two component version of ballic. It first calculates an equilibrium model for a given material and then
# generates a particle representation.

multi: ballic.multi.o $(TIPSY_OBJECTS) $(TILL_OBJECTS) $(FORTRAN_OBJECTS)
	$(CC) $(CFLAGS) $^ -o $(BALLIC).$@ $(LIBS)

# A two component model, that has a low density atmosphere of the mantle material

multi.atm: ballic.multi.atm.o $(TIPSY_OBJECTS) $(TILL_OBJECTS) $(FORTRAN_OBJECTS)
	$(CC) $(CFLAGS) $^ -o $(BALLIC).$@ $(LIBS)

# Calculates equilibrium models for a given material and different initial densities and internal energies.

modelsolve: modelsolve.o $(TILL_OBJECTS)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

tipsy.o:
	$(CC) $(CFLAGS) $(INC_TIPSYDEFS) $^ -o $@ -c $(TIPSYSRC)


clean:
	rm -f $(EXE) $(OBJECTS) 

cleanall:
	rm -f $(EXE) $(OBJECTS) $(TIPSY_OBJECTS) $(TILL_OBJECTS) $(FORTRAN_OBJECTS)
