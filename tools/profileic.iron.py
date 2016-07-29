"""
This script produces a u(rho) plot for a tipsy binary file.
"""
from matplotlib import *
rcParams['font.family'] = 'serif'
from numpy import *
from matplotlib.pyplot import *
from sys import *

if len(argv) != 3:
		print "Usage: profileden.iron.py tipsyfile model"
		exit(1)

tipsyfile = argv[1]
modelfile = argv[2]

# The cold curve for granite
coldcurve = loadtxt('~/code/ballic/tools/cold_iron.txt')
model = loadtxt(modelfile)

# Load cold curve
rhocold = coldcurve[:,0]
ucold = coldcurve[:,1]

# Load model
rhomodel = model[:,1]
umodel = model[:,3]

# Values for iron (in code units!)
rho0=21.331
us=1.42525
us2=8.4515

# Load profile
profile = loadtxt(tipsyfile)

rho = profile[:,1]
u = profile[:,2]

# First plot u(rho)
xmax = max(rho)*1.1
ymax = max(u)*1.1

xlim(0,xmax)
ylim(0,ymax)

# title(r'Target (ISPH, modified density')
xlabel('Density')
ylabel('Internal energy')

plot(rhocold,ucold,'-',color='red',linewidth=2,label='Cold curve (T=0)')
plot(rhomodel,umodel,'-',color='blue',linewidth=2,label='Model')
scatter(rho,u,s=5,color='magenta',label='Target (SPH)')

plot([0,rho0],[us,us],'r--')
plot([0,rho0],[us2,us2],'r--')
plot([rho0,rho0],[0,ymax],'r--')

legend(loc='best')

savefig(tipsyfile+'.png')
close()

print "Done."
