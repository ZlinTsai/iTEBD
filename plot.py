import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.mathtext
from   matplotlib.ticker import MultipleLocator

# set picture parameter
fig, ax = plt.subplots()
plt.tick_params(which='major', direction='in', length=7)
plt.tick_params(top='on',bottom='on',left='on',right='on')
#mlx = MultipleLocator(0.5)
#mly = MultipleLocator(0.05)
#ax.xaxis.set_minor_locator(mlx)
#ax.yaxis.set_minor_locator(mly)
plt.tick_params(which='minor', direction='in', length=4)

# set dataframe
df = pd.read_csv('database/bondDim_30_J1.csv')
df2 = pd.read_csv('database/bondDim_30_J-1.csv')

# plot
#plt.plot(df['h'], df['En'], 'r-o', label=r'Bond dimension = 30')
#plt.plot(df['h'], 0.5*abs(df['Sz_0'] + df['Sz_1']), 'b-o', label=r'J = 1, bond dimension = 30')
#plt.plot(df2['h'], 0.5*abs(df2['Sz_0'] + df2['Sz_1']), 'g-o', label=r'J = -1, bond dimension = 30')
plt.plot(df['h'], 0.5*abs(df['Sx_0'] + df['Sx_1']), 'r-o', label=r'Bond dimension = 30')

#
plt.xlabel(r'Transverse magnetic field (h)', fontsize=16) #$\Delta$
#plt.ylabel(r'Ground state energy per site', fontsize=18)
#plt.ylabel(r'$|Mz| = 0.5*|<\sigma^z_i> + <\sigma^z_{i+1}>|$', fontsize=14)
plt.ylabel(r'$|Mx| = 0.5*|<\sigma^x_i> + <\sigma^x_{i+1}>|$', fontsize=14)
#plt.title(r'Energy v.s. Transverse magnetic field', fontsize=20)
#plt.title(r'Ferromagnetic', fontsize=20)
plt.title(r'paramagnetic', fontsize=20)
plt.legend() #show legend
plt.grid(linestyle='--') #grid line
#plt.show()

fig.savefig('fig_6.PNG', format='png', dpi=2000)