"""
This script produces a plot of the cold curve for the Tillotson EOS.
"""
from matplotlib import *
rcParams['font.family'] = 'serif'
from matplotlib.patches import *
from numpy import *
from matplotlib.pyplot import *

data = loadtxt('testsplintrho.txt')
# Load the whole lookup table
data2 = loadtxt('lookup.txt')

rhocold = data2[:,0]
ucold = data2[:,1]

rho = data[:,0]
u = data[:,1]
u_cubicint = data[:,2]

# Values for granite
us = 3.5
us2 = 18.0
rho0 = 7.33

xmax = 25
ymax = 25
xlim(0,xmax)
ylim(0,ymax)

xlabel('Density')
ylabel('Internal energy')

# Plot the cold curve
plot(rhocold,ucold,'-',color='red',linewidth=2,label='Cold curve (T=0)')

# Plot different isentropes (Lookup(i, j) = Lookup(rho,v))
for i in arange(2,1000,1):
		plot(data2[:,0],data2[:,i],'-',color='green',label='Table')

plot(rho,u,'.',color='blue',markersize=4,label='Lookup')
plot(rho,u_cubicint,'.',color='red',markersize=2,label='Lookup')
show()

#savefig('testsplintrho.pdf')
savefig('testsplintrho.png')

