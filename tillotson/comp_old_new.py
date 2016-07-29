"""
This script produces a plot of the cold curve for the Tillotson EOS.
"""
from matplotlib import *
rcParams['font.family'] = 'serif'
from matplotlib.patches import *
from numpy import *
from matplotlib.pyplot import *

# Load the whole lookup table
data1 = loadtxt('lookup.txt')
data3 = loadtxt('../tillotson.prasenjit.backup/lookup.txt')
# Load interpolated values
data2 = loadtxt('testsplint.new.txt')
data4 = loadtxt('../tillotson.prasenjit.backup/testsplint.old.txt')

#v = data[:,0]
#u = data[:,1]

subplot(1,2,1)

# Plot the lookup table
for i in range(1,size(data1[:,0]),1):
		plot(data1[:,0],data1[:,i],'-',color='red',markersize=1,label='Table (new)')

for i in range(1,size(data2[0,:]),1):
		plot(data2[:,0],data2[:,i],'-',color='green',markersize=1,label='Lookup (new)')

subplot(1,2,2)

# Plot the lookup table
for i in range(1,size(data3[:,0]),1):
		plot(data3[:,0],data3[:,i],'-',color='cyan',markersize=1,label='Table (old)')

for i in range(1,size(data4[0,:]),1):
		plot(data4[:,0],data4[:,i],'-',color='blue',markersize=1,label='Lookup (new)')

annotate(r'dv=0.252525, drho=0.0250084', xy = (1,1), xytext = (100, 25*0.1))

title('The Tillotson equation of state')
xlabel('Density')
ylabel('Internal energy')

#xmax = 25
#ymax = 25
#xlim(0,xmax)
#ylim(0,ymax)

#plot(rho,u,'.',color='blue',markersize=1,label='Lookup')

#plot(rho,u,'-',color='blue',linewidth=1,label='Lookup table')

#fill_between(rhocold,ucold,color='orange')

#plot([0,rho0],[us,us],'r--')
#plot([0,rho0],[us2,us2],'r--')
#plot([rho0,rho0],[0,25],'r--')

#xticks([0,rho0,2*rho0,3*rho0],[r'$0$',r'$\rho_0$',r'$2\rho_0$',r'$3\rho_0$'])
#xticks([0,rho0],[r'$0$',r'$\rho_0$'],size='large')
#yticks([0,us,us2],[r'$0$',r'$u_{IV}$',r'$u_{CV}$'],size='large')

#savefig('lookup.pdf')
savefig('testsplint.compare.png')

