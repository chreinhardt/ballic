
#objects = tillotson.o interpol/coeff.o interpol/interpol.o interpol/brent.o nr/tridag.o nr/spline.o nr/splint.o nr/nrutil.o

objects = tillotson.o tillinitlookup.o tillsplint.o interpol/coeff.o interpol/interpol.o interpol/brent.o nr/nrcubicspline.o nr/nrutil.o

exe = table cold coldlookup pressureoldnew lookup lookup_cold testu1 testspline testsplint testsplint2 testsplinerho testsplintrho testsplinev testsplintv testcubicintrho testlookupucold testudrho testudv testgrid testpolyv printderiv printpress pressneg testisintable testrhomin testoutofbounds testsolvebc calcisentrope

CFLAGS ?= -O3

default:
	@echo "Please specify which tool you want to make."

all:
	default

table: table.o $(objects)
	cc -o table table.o $(objects) -lm

cold: cold.o $(objects)
	cc -o cold cold.o $(objects) -lm

coldlookup: coldlookup.o $(objects)
	cc -o coldlookup coldlookup.o $(objects) -lm

pressureoldnew: pressureoldnew.o $(objects)
	cc -o pressureoldnew pressureoldnew.o $(objects) -lm

lookup: lookup.o $(objects)
	cc -o lookup lookup.o $(objects) -lm

lookup_cold: lookup_cold.o $(objects)
	cc -o lookup_cold lookup_cold.o $(objects) -lm

testu1: testu1.o $(objects)
	cc -o testu1 testu1.o $(objects) -lm

testspline: testspline.o $(objects)
	cc -o testspline testspline.o $(objects) -lm

testsplint: testsplint.o $(objects)
	cc -o testsplint testsplint.o $(objects) -lm

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

testrhomin: testrhomin.o $(objects)
	cc -o testrhomin testrhomin.o $(objects) -lm

testoutofbounds: testoutofbounds.o $(objects)
	cc -o testoutofbounds testoutofbounds.o $(objects) -lm

testsolvebc: testsolvebc.o $(objects)
	cc -o testsolvebc testsolvebc.o $(objects) -lm

calcisentrope: calcisentrope.o $(objects)
	cc -o calcisentrope calcisentrope.o $(objects) -lm

clean:
	rm $(objects)

cleanall:
	rm $(exe) $(objects) *.o
