"""
This script produces a plot of the cold curve for the Tillotson EOS.
"""
from matplotlib import *
rcParams['font.family'] = 'serif'
from matplotlib.patches import *
from numpy import *
from matplotlib.pyplot import *

# Load the whole lookup table
data = loadtxt('lookup.txt')
# Load first derivative
data2 = loadtxt('printpress.txt')

#v = data[:,0]
#u = data[:,1]

# Plot the lookup table
#for i in range(1,size(data[:,0]),1):
#		plot(data[:,0],data[:,i],'-',color='red',markersize=1,label='Table')

for i in range(1,size(data2[0,:]),1):
		plot(data2[:,0],data2[:,i],'-',color='green',markersize=1,label='P(rho,u)')

#plot(data2[:,0],data2[:,99],'-',color='blue',linewidth=1,label='u1(v=vmax)')
title('The Tillotson equation of state')
xlabel('Density')
ylabel('Pressure')
#legend(loc='best')

#annotate(r'dv=0.252525, drho=0.0250084', xy = (1,1), xytext = (100, 25*0.1))

#xmax = 25
#ymax = 25
#xlim(0,xmax)
#ylim(0,ymax)

#plot(rho,u,'.',color='blue',markersize=1,label='Lookup')



#fill_between(rhocold,ucold,color='orange')

#plot([0,rho0],[us,us],'r--')
#plot([0,rho0],[us2,us2],'r--')
#plot([rho0,rho0],[0,25],'r--')

#xticks([0,rho0,2*rho0,3*rho0],[r'$0$',r'$\rho_0$',r'$2\rho_0$',r'$3\rho_0$'])
#xticks([0,rho0],[r'$0$',r'$\rho_0$'],size='large')
#yticks([0,us,us2],[r'$0$',r'$u_{IV}$',r'$u_{CV}$'],size='large')

#savefig('lookup.pdf')
savefig('testsplint.png')

