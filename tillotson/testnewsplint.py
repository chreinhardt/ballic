"""
This script produces a plot of the cold curve for the Tillotson EOS.
"""
from matplotlib import *
rcParams['font.family'] = 'serif'
from numpy import *
from matplotlib.pyplot import *

# Load the whole lookup table
data1 = loadtxt('lookup.old.txt')
# Load interpolated values
data2 = loadtxt('testnewsplint.old.txt')

# Load the whole lookup table
data3 = loadtxt('lookup.txt')
# Load interpolated values
data4 = loadtxt('testnewsplint.txt')
#v = data[:,0]
#u = data[:,1]

xmax = 25
ymax = 25
xlim(0,xmax)
ylim(0,ymax)

#title('The Tillotson equation of state')
xlabel('Density')
ylabel('Internal energy')


# Plot the lookup table for the old interpolator
for i in range(1,size(data1[:,0]),1):
		plot(data1[:,0],data1[:,i],'-',color='red',markersize=1,label='Table')

for i in range(1,size(data2[0,:]),1):
		plot(data2[:,0],data2[:,i],'-',color='green',markersize=1,label='Lookup (old)')

# Plot the lookup table for the new interpolator
for i in range(1,size(data3[:,0]),1):
		plot(data3[:,0],data3[:,i],'-',color='black',markersize=1,label='Table')

for i in range(1,size(data4[0,:]),1):
		plot(data4[:,0],data4[:,i],'--',color='orange',markersize=1,label='Lookup (new)')

show()

savefig('testnewsplint.png')

