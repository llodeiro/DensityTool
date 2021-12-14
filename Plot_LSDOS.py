import numpy as np
import matplotlib.pyplot as plt
from pylab import *

### VARIABLES NEEDED FROM PREVIOUS VASP CALCULATION
Ni          = 			# Number of points along r_i
Lbox        = 		 	# Box length along r_i
Eminplotted = 			# Minimum energy range of computed LDOS (EMIN)
Emaxplotted =			# Maximum energy range of computed LDOS (EMAX)
Eshift      = 			# Energy shift, added to energies (typically the Fermi energy)

### FIGURE VARIABLES
cm = "bwr"				# Coloring type
fig, ax = plt.subplots(figsize = (12,6))# Figure size plotted
xlabel('z (''$\AA$'')', fontsize=21)	# X axis name
ylabel(r'$E$ (eV)', fontsize=21)	# Y axis name
xticks(size=21)
yticks(size=21)
vmin = 				# Minimum of color scale
vmax = 				# Maximum of color scale
ymin = 				# Y axis minimum
ymax = 				# Y axis maximum
xmin = 				# X axis minimum
xmax = 				# X axis maximum

z = np.loadtxt('LSDOS.R3.all.dat') 
#########################################
for i in range(len(z[:,2])):
    z[i,2] = max(z[i,2],0.000000001)
to_plot = np.log(z[:,2])
large = z.shape[0]			
nE    = large/Ni			
RVminEplotted = Eminplotted+Eshift
RVmaxEplotted = Emaxplotted+Eshift	
to_plot.shape=(nE,Ni)
x=np.linspace(0, Lbox,Ni)
y=np.linspace(RVminEplotted, RVmaxEplotted, nE)
bc=ax.pcolor(x, y, to_plot,cmap=cm,vmin=vmin, vmax=vmax)
cbar=plt.colorbar(bc,ticks=[vmin,vmax], shrink=0.7, pad = 0.1)
cbar.ax.set_yticklabels([vmin,vmax],size=21)
cbar.ax.set_ylabel(r"$\log(LSDOS)$", rotation=270,size=21)
ax.set_ylim(ymin,ymax)
ax.set_xlim(xmin,xmax)
plt.show()

