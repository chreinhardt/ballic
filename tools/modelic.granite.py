"""
This script produces a u(rho) plot for a planetary model.
"""
from matplotlib import *
rcParams['font.family'] = 'serif'
from numpy import *
from matplotlib.pyplot import *
from sys import *

if len(argv) != 2:
		print "Usage: modelic.py model"
		exit(1)

modelfile = argv[1]
model = loadtxt(modelfile)
i
"""
Probably you will have to change this path!
"""
# The cold curve for granite
coldcurve = loadtxt('~/code/ballic/tools/cold_granite.txt')

# Load cold curve
rhocold = coldcurve[:,0]
ucold = coldcurve[:,1]

# Load model
R = model[:,0]
rho = model[:,1]
u = model[:,3]


# Values for granite
us = 3.5
us2 = 18.0
rho0 = 7.33

# First plot u(rho)
xmax = max(rho)*1.1
ymax = max(u)*1.1

xlim(0,xmax)
ylim(0,ymax)

xlabel('Density')
ylabel('Internal energy')

plot(rhocold,ucold,'-',color='red',linewidth=2,label='Cold curve (T=0)')
plot(rho,u,'-',color='blue',linewidth=2,label='Model')

plot([0,rho0],[us,us],'r--')
plot([0,rho0],[us2,us2],'r--')
plot([rho0,rho0],[0,ymax],'r--')

legend(loc='best')


savefig(modelfile+'.ic.png')
close()

print "Done."
