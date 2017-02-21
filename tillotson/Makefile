
#objects = tillotson.o interpol/coeff.o interpol/interpol.o interpol/brent.o nr/tridag.o nr/spline.o nr/splint.o nr/nrutil.o
#objects = tillotson.o tillwoolfson.o tillinitlookup.o tillsplint.o interpol/coeff.o interpol/interpol.o interpol/brent.o nr/nrcubicspline.o nr/nrutil.o
objects = tillotson.o tillinitlookup.o tillsplint.o interpol/brent.o nr/nrcubicspline.o nr/nrutil.o

exe = table pressureoldnew lookup lookup_cold testu1 testspline testsplint testnewsplint testsplint2 testsplinerho testsplintrho testsplinev testsplintv testcubicintrho testlookupucold testudrho testudv testgrid testpolyv printderiv printpress pressneg testisintable testisbelowcoldcurve testrhomin testoutofbounds testsolvebc calcisentrope testrhoptemp calcpressure

defs = -DTILL_PRESS_NP

CFLAGS ?= -O3 $(defs)

default:
	@echo "Please specify which tool you want to make."

all:
	default

# Debug the function tillLookupU() by comparing the results of the lookup to a
# direct calculation using tillCalcU() which simply solves the ODE using a RK4
# integrator.
table: table.o $(objects)
	cc -o table table.o $(objects) -lm

# Compare the results of tillPressureSoundold() that contains the old code we used in
# Gasoline and tillPressureSound()
pressureoldnew: pressureoldnew.o $(objects)
	cc -o pressureoldnew pressureoldnew.o $(objects) -lm

# Generate a lookup table for a given material
lookup: lookup.o $(objects)
	cc -o lookup lookup.o $(objects) -lm

# Make a lookup table for the cold curve
lookup_cold: lookup_cold.o $(objects)
	cc -o lookup_cold lookup_cold.o $(objects) -lm

testu1: testu1.o $(objects)
	cc -o testu1 testu1.o $(objects) -lm

testspline: testspline.o $(objects)
	cc -o testspline testspline.o $(objects) -lm

# Code for debugging the interpolation function tillCubicInt()
testsplint: testsplint.o $(objects)
	cc -o testsplint testsplint.o $(objects) -lm

# Pretty much the same but it compares the old with the new interpolator 
testnewsplint: testnewsplint.o $(objects)
	cc -o testnewsplint testnewsplint.o $(objects) -lm

testsplint2: testsplint2.o $(objects)
	cc -o testsplint2 testsplint2.o $(objects) -lm

testsplinerho: testsplinerho.o $(objects)
	cc -o testsplinerho testsplinerho.o $(objects) -lm

testsplintrho: testsplintrho.o $(objects)
	cc -o testsplintrho testsplintrho.o $(objects) -lm

testsplinev: testsplinev.o $(objects)
	cc -o testsplinev testsplinev.o $(objects) -lm

testsplintv: testsplintv.o $(objects)
	cc -o testsplintv testsplintv.o $(objects) -lm

testcubicintrho: testcubicintrho.o $(objects)
	cc -o testcubicintrho testcubicintrho.o $(objects) -lm

testlookupucold: testlookupucold.o $(objects)
	cc -o testlookupucold testlookupucold.o $(objects) -lm

testudrho: testudrho.o $(objects)
	cc -o testudrho testudrho.o $(objects) -lm

testudv: testudv.o $(objects)
	cc -o testudv testudv.o $(objects) -lm

testgrid: testgrid.o $(objects)
	cc -o testgrid testgrid.o $(objects) -lm

testpolyv: testpolyv.o $(objects)
	cc -o testpolyv testpolyv.o $(objects) -lm

printderiv: printderiv.o $(objects)
	cc -o printderiv printderiv.o $(objects) -lm

printpress: printpress.o $(objects)
	cc -o printpress printpress.o $(objects) -lm

pressneg: pressneg.o $(objects)
	cc -o pressneg pressneg.o $(objects) -lm

testisintable: testisintable.o $(objects)
	cc -o testisintable testisintable.o $(objects) -lm

testisbelowcoldcurve: testisbelowcoldcurve.o $(objects)
	cc -o testisbelowcoldcurve testisbelowcoldcurve.o $(objects) -lm

testrhomin: testrhomin.o $(objects)
	cc -o testrhomin testrhomin.o $(objects) -lm

testoutofbounds: testoutofbounds.o $(objects)
	cc -o testoutofbounds testoutofbounds.o $(objects) -lm

testsolvebc: testsolvebc.o $(objects)
	cc -o testsolvebc testsolvebc.o $(objects) -lm

calcisentrope: calcisentrope.o $(objects)
	cc -o calcisentrope calcisentrope.o $(objects) -lm

testrhoptemp: testrhoptemp.o $(objects)
	cc -o testrhoptemp testrhoptemp.o $(objects) -lm

calcpressure: calcpressure.o $(objects)
	cc -o calcpressure calcpressure.o $(objects) -lm

clean:
	rm $(objects)

cleanall:
	rm $(exe) $(objects) *.o
