import numpy as np
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

dn = "./data/"
nl = np.load(dn+"vortex_nl.npz",allow_pickle=True)

M,N = 24,24

fig,ax = plt.subplots(figsize=(6,3))

im = ax.imshow(nl['Vxy'][:,:,-1],interpolation="bicubic",cmap="seismic",origin="lower")
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='2%', pad=0.2)
fig.colorbar(im, cax=cax, orientation='vertical')

ax.set_ylabel(r'$\zeta(y,t)$',fontsize=12)
ax.set_xlabel(r'$t$',fontsize=12)

ax.set_xticks([0,M-1,2*M-2])
ax.set_yticks([0,N-1,2*N-2])
ax.set_xticklabels([r'$0$',r'$2\pi$',r'$4\pi$'],fontsize=10)
ax.set_yticklabels([r'$0$',r'$\pi$',r'$2\pi$'],fontsize=10)

plt.savefig(dn+'vortex_nl.png',bbox_inches='tight',dpi=512)
