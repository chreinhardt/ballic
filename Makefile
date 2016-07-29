
objects = ballic.o tipsy.o tillotson/tillotson.o tillotson/tillinitlookup.o tillotson/tillsplint.o tillotson/interpol/coeff.o tillotson/interpol/interpol.o tillotson/interpol/brent.o tillotson/nr/nrcubicspline.o tillotson/nr/nrutil.o
fortran_objects = icosahedron.o

CFLAGS ?= -O3

ballic: $(objects)
	cc -o ballic $(objects) $(fortran_objects) -lm

all: ballic

clean:
	rm $(objects)

