
objects = ballic.o ballic.single.o ballic.multi.o
tipsy_objects = tipsy.o 
till_objects = tillotson/tillotson.o tillotson/tillinitlookup.o tillotson/tillsplint.o tillotson/interpol/coeff.o tillotson/interpol/interpol.o tillotson/interpol/brent.o tillotson/nr/nrcubicspline.o tillotson/nr/nrutil.o
fortran_objects = icosahedron.o

exe = ballic ballic.single ballic.multi

CFLAGS ?= -O3

default:
	@echo "Please specify, which tool to make (e.g., ballic, single, multi)."

all:
	default

ballic: ballic.o $(tipsy_objects)
	cc -o ballic ballic.o $(tipsy_objects) $(fortran_objects) -lm

single: ballic.single.o $(tipsy_objects) $(till_objects)
	cc -o ballic.single ballic.single.o $(tipsy_objects) $(till_objects) $(fortran_objects) -lm

clean:
	rm $(exe) $(objects)

cleanall:
	rm  $(exe) $(objects) $(tipsy_objects) $(till_objects)
