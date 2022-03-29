import numpy as np
import matplotlib.pyplot as plt

ne = 240
nz = 40
ef = 5.962
d = 3
data = np.zeros((3,ne,nz))
for ie in range(ne):
    fname = "LDOSMAG.R"+str(d)+"."+str(ie+1).zfill(4)+".dat"
    with open(fname, 'r') as f:
        text = f.readlines()
    # number of lines
    nlin=len(text)
    for il in range(nlin):
        text[il] = text[il].split()
        data[0,ie,il] = float(text[il][0])
        data[1,ie,il] = float(text[il][1])-ef
        data[2,ie,il] = max(float(text[il][2]),0.0000000001)

# Plot
cm = "viridis"
fig, ax1 = plt.subplots(figsize = (6,6))

vmin = 0.0
vmax = 0.05
bc1=ax1.pcolormesh(data[0],data[1],data[2],cmap=cm,vmin=vmin, vmax=vmax)
cbar1=fig.colorbar(bc1,ticks=[vmin,vmax], shrink=0.7, orientation = 'horizontal',ax=ax1)
cbar1.ax.set_xticklabels([vmin,vmax],size=20)
cbar1.ax.set_xlabel(r"$\mathrm{LDOS}$",size=20)

ax1.set_xticks([])
ax1.set_ylim([-1.15,1.24])
ticks = [-1.0,0.0,1.0]
ax1.set_yticks(ticks)
ax1.set_yticklabels(ticks,size=20)
ax1.set_ylabel(r'$E-E_{\mathrm{F}}$ (eV)',size=20)

# export
fig.tight_layout()
#plt.show()
fig.savefig("ldos.png")
