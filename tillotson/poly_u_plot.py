"""
This script produces a plot of the cold curve for the Tillotson EOS.
"""
from matplotlib import *
rcParams['font.family'] = 'serif'
from matplotlib.patches import *
from numpy import *
from matplotlib.pyplot import *

# Load the data
data = loadtxt('poly_u.txt')

i = data[:,0]
v = data[:,1]
dv = data[:,2]

plot(i, dv,'-',color='red')

#title('The Tillotson equation of state')
#xlabel('Density')
#ylabel('Internal energy')

savefig('poly_u_plot.png')

