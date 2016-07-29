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
# Load interpolated values
data2 = loadtxt('testsplint2.txt')

rho1_start = data2[:,0]
u1_start = data2[:,1]
rho2_start = data2[:,2]
u2_start = data2[:,3]

rho1_end = data2[:,4]
u1_end = data2[:,5]
rho2_end = data2[:,6]
u2_end = data2[:,7]

# Plot the lookup table
for i in arange(1,size(data[:,0])+1,1):
#		plot(data[:,0],data[:,i],'-',color='blue',markersize=1,label='Table')
	plot(data[:,0],data[:,i],'-',color='cyan',alpha=0.2,markersize=0.1)

for i in arange(1,size(data[:,0])+1,10):
#		plot(data[:,0],data[:,i],'-',color='blue',markersize=1,label='Table')
		plot(data[:,0],data[:,i],'-',color='magenta',alpha=0.5,markersize=1)

scatter(rho1_start,u1_start,s=16,color='green',alpha=1,label='(rho0,v)')
scatter(rho2_start,u2_start,s=16,color='red',alpha=1,label='(rho2,u2)')
scatter(rho2_end,u2_end,s=9,color='blue',alpha=0.5,label='(rho1,u1)')

xmax = 25
ymax = 25
xlim(0,xmax)
ylim(0,ymax)

#title('The Tillotson equation of state')
xlabel('Density')
ylabel('Internal energy')

#plot(rho,u,'-',color='blue',linewidth=1,label='Lookup table')
#fill_between(rhocold,ucold,color='orange')

#plot([0,rho0],[us,us],'r--')
#plot([0,rho0],[us2,us2],'r--')
#plot([rho0,rho0],[0,25],'r--')

#xticks([0,rho0,2*rho0,3*rho0],[r'$0$',r'$\rho_0$',r'$2\rho_0$',r'$3\rho_0$'])
#xticks([0,rho0],[r'$0$',r'$\rho_0$'],size='large')
#yticks([0,us,us2],[r'$0$',r'$u_{IV}$',r'$u_{CV}$'],size='large')

#savefig('lookup.pdf')
savefig('testsplint2.png')

