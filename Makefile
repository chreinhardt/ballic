
objects = ballic.o ballic.single.o ballic.multi.o
tipsy_objects = tipsy.o 
till_objects = tillotson/tillotson.o tillotson/tillinitlookup.o tillotson/tillsplint.o tillotson/interpol/coeff.o tillotson/interpol/interpol.o tillotson/interpol/brent.o tillotson/nr/nrcubicspline.o tillotson/nr/nrutil.o
fortran_objects = icosahedron.o
FC := gfortran

exe = ballic ballic.single ballic.multi

CFLAGS ?= -O3 -march=native
FFLAGS ?= $(CFLAGS)

default:
	@echo "Please specify, which tool to make (e.g., ballic, single, multi)."

all:
	default

# Standard version of ballic that reads an equilibrium model from a file and generates an particle distribution,
# that matches the desired density profile. Currently not reliably working for multi component models.
ballic: ballic.o $(tipsy_objects) $(fortran_objects)
	cc -o ballic ballic.o $(tipsy_objects) $(fortran_objects) -lm

# Single component version of ballic. It first calculates an equilibrium model for a given material and then
# generates a particle representation.
single: ballic.single.o $(tipsy_objects) $(till_objects) $(fortran_objects)
	cc -o ballic.single ballic.single.o $(tipsy_objects) $(till_objects) $(fortran_objects) -lm

# Two component version of ballic. It first calculates an equilibrium model for a given material and then
# generates a particle representation.
multi: ballic.multi.o $(tipsy_objects) $(till_objects) $(fortran_objects)
	cc -o ballic.multi ballic.multi.o $(tipsy_objects) $(till_objects) $(fortran_objects) -lm

# A two component model, that has a low density atmosphere of the mantle material
multi.atm: ballic.multi.atm.o $(tipsy_objects) $(till_objects) $(fortran_objects)
	cc -o ballic.multi.atm ballic.multi.atm.o $(tipsy_objects) $(till_objects) $(fortran_objects) -lm

clean:
	rm -f $(exe) $(objects) 

cleanall:
	rm -f $(exe) $(objects) $(tipsy_objects) $(till_objects) $(fortran_objects)
