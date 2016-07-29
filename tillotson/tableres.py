"""
This script produces a plot of the cold curve for the Tillotson EOS.
"""
#import matplotlib as mpl; mpl.rcParams['font.family'] = 'serif'
from matplotlib import *
rcParams['font.family'] = 'serif'
from matplotlib.patches import *
from numpy import *
from matplotlib.pyplot import *

#isentrope = loadtxt('/home/ics/creinh/code/ic/adiabat.txt')
coldcurve = loadtxt('/home/ics/creinh/code/condensed/coldcurve/coldcurve.txt')
data1 = loadtxt('table30.txt')
data2 = loadtxt('table100.txt')
data3 = loadtxt('table1000.txt')

#rhoi = isentrope[:,0]
#ui = isentrope[:,1]
rhocold = coldcurve[:,0]
ucold = coldcurve[:,1]

rho1 = data1[:,0]
u1 = data1[:,1]

rho2 = data2[:,0]
u2 = data2[:,1]


# Values for granite
us = 3.5
us2 = 18.0
rho0 = 7.33

xmax = 25
ymax = 25
xlim(0,xmax)
ylim(0,ymax)

#title(r'The internal energy profile of the target')
xlabel('Density')
ylabel('Internal energy')

#plot(rhoi,ui,'.',color='green',markersize=0.1)
plot(rhocold,ucold,'-',color='red',linewidth=2,label='Cold curve (T=0)')
#scatter(rho,u,s=5,color='blue',label='Tillotson.c')
plot(rho,u,'-',color='blue',linewidth=1,label='Tillotson.c')

plot([0,rho0],[us,us],'r--')
plot([0,rho0],[us2,us2],'r--')
plot([rho0,rho0],[0,25],'r--')

legend()

#savefig('target100kcondensedc2.old.ic.pdf')
#savefig('target100kcondensedc2.old.ic.png')
savefig('table.png')

