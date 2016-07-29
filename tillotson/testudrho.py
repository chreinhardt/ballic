"""
This script produces a plot of the cold curve for the Tillotson EOS.
"""
from matplotlib import *
rcParams['font.family'] = 'serif'
from matplotlib.patches import *
from numpy import *
from matplotlib.pyplot import *

data = loadtxt('testudrho.txt')

# Load data
i = data[:,0]
j = data[:,1]
rho = data[:,2]
u = data[:,3]
dudrho = data[:,4]
dudrho_a = data[:,5]

delta = (abs(dudrho - dudrho_a))/max(abs(dudrho - dudrho_a))
#delta = (abs(dudrho - dudrho_a))

print "delta max:"+str(max(delta))+" at "+str(where(delta == max(delta)))

print where(delta>0.01)

print "Plotting"

scatter(rho,u,s=1,c=delta,cmap='rainbow',linewidth=0)
#scatter(rho[where(delta>0.01)],u[where(delta>0.01)],s=1,c=delta[where(delta>0.01)],cmap='rainbow',linewidth=0)
colorbar()

#savefig('lookup.pdf')
savefig('testudrho.png')

