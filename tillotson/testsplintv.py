"""
This script produces a plot of the cold curve for the Tillotson EOS.
"""
from matplotlib import *
rcParams['font.family'] = 'serif'
from matplotlib.patches import *
from numpy import *
from matplotlib.pyplot import *

# This might be useful later to debug our look up table.
#isentrope = loadtxt('/home/ics/creinh/code/ic/adiabat.txt')
coldcurve = loadtxt('/home/ics/creinh/code/condensed/coldcurve/coldcurve.txt')

data = loadtxt('testsplintv.txt')
# Load the whole lookup table
data2 = loadtxt('lookup.txt')

rhocold = coldcurve[:,0]
ucold = coldcurve[:,1]

rho = data[:,0]
v = data[:,1]
u = data[:,2]

# Values for granite
us = 3.5
us2 = 18.0
rho0 = 7.33

xmax = 25
ymax = 25
#xlim(0,xmax)
#ylim(0,ymax)

#title('The Tillotson equation of state')
xlabel('Density')
ylabel('Internal energy')

plot(rhocold,ucold,'-',color='red',linewidth=2,label='Cold curve (T=0)')

# Plot different isentropes (Lookup(i, j) = Lookup(rho,v))
#for i in arange(1,1000,100):
#		plot(data2[:,0],data2[:,i],'.',color='green',markersize=1,label='Table')

for i in arange(1,size(data2[0,:]),1):
		plot(data2[:,0],data2[:,i],'.',color='green',markersize=1,label='Table')

#rho = zeros(size(v))
#for i in arange(0,size(v),1):
#		rho[i]=rho0

plot(rho,u,'.',color='blue',markersize=1,label='Lookup')

#plot(rho,u,'-',color='blue',linewidth=1,label='Lookup table')

#fill_between(rhocold,ucold,color='orange')

plot([0,rho0],[us,us],'r--')
plot([0,rho0],[us2,us2],'r--')
plot([rho0,rho0],[0,25],'r--')

#xticks([0,rho0,2*rho0,3*rho0],[r'$0$',r'$\rho_0$',r'$2\rho_0$',r'$3\rho_0$'])
#xticks([0,rho0],[r'$0$',r'$\rho_0$'],size='large')
#yticks([0,us,us2],[r'$0$',r'$u_{IV}$',r'$u_{CV}$'],size='large')

#savefig('lookup.pdf')
savefig('testsplintv.png')

