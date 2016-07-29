"""
This script produces a plot of the cold curve for the Tillotson EOS.
"""
from matplotlib import *
rcParams['font.family'] = 'serif'
from matplotlib.patches import *
from numpy import *
from matplotlib.pyplot import *

nTable = 1000
vmax = 1200.0
n = 1000

#i = range(0,nTable,1)
#v = range(0,nTable,1)
#v = zeros(nTable)
v_old = 0.0

print "# i          v         dv"
for i in range(0,nTable,1):
	v = vmax/((nTable-1)**n)*i**n
	dv = abs(v-v_old)
	print str(i)+"   "+str(v)+"   "+str(dv)
	v_old = v


