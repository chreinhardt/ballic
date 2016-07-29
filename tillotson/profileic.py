"""
This script produces a plot of the sound speed for the Tillotson EOS.
"""
#import matplotlib as mpl; mpl.rcParams['font.family'] = 'serif'
from matplotlib import *
rcParams['font.family'] = 'serif'
from matplotlib.patches import *
from numpy import *
from matplotlib.pyplot import *

cs = loadtxt('cshot.txt')

rho = cs[:,0]
u = cs[:,1]
csold = cs[:,2]
cslinold = cs[:,3]
cslinnew = cs[:,4]

#xmax = 25
#ymax = 25
#xlim(0,xmax)
#ylim(0,ymax)

title(r'The Sound Speed for the Tillotson EOS (u=25)')
xlabel('Density')
ylabel(r'$c_s^2$')

#plot(rho,u,'.',color='green',markersize=0.1)

plot(rho,csold,'-',color='red',linewidth=1,label='Gasoline')
plot(rho,cslinold,'-',color='green',linewidth=1,label='Lin. interpolation')
plot(rho,cslinnew,'-',color='blue',linewidth=1,label='tillotson.c')

legend(loc='lower right')

#savefig('target100kcondensedc2.old.ic.pdf')
#savefig('target100kcondensedc2.old.ic.png')
savefig('target100kcondensedc2.linc2.ic.png')

