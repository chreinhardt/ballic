"""
This script produces an omega(rho) plot to check its behavior in the limit where rho goes to zero.
"""
from matplotlib import *
rcParams['font.family'] = 'serif'
from matplotlib.patches import *
from numpy import *
from matplotlib.pyplot import *

# Make an array
rho_start = 5.0
rho_end = 0.0
drho = 0.1

rho = np.array(arange(rho_end,rho_start,drho),dtype='double')

# Internal energy assumed constant (but of course u = u(rho))
u = 8.6

# Some material constant for the Tillotson EOS
rho0 = 7.33
u0 = 16

eta = rho/rho0

omega0 = u/(u0*eta*eta)+1.0

plot(rho, omega0,'-',color='blue')

xmax = rho_start
xlim(0,xmax)

#title('The Tillotson equation of state')
xlabel('Density')
ylabel(r'$\omega_0$')

savefig('lim_omega.png')

