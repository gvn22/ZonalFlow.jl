import numpy as np
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from mpl_toolkits.axes_grid1 import make_axes_locatable

dn = "./data/"
nl = np.load(dn+"jets_nl.npz",allow_pickle=True)

M,N = 16,16

# Plot spectra
fig,ax = plt.subplots(1,3,figsize=(15,6))

ax[0].set_title(r'$\sqrt{\varepsilon}\widehat{Q}(k_x,k_y)$',fontsize=12)
im = ax[0].imshow(nl['F'],interpolation="none",cmap="Greys",origin="lower",norm=LogNorm(vmin=1e-6,vmax=1e-1))
divider = make_axes_locatable(ax[0])
cax = divider.append_axes('right', size='5%', pad=0.2)
fig.colorbar(im, cax=cax, orientation='vertical')

ax[1].set_title(r'$\widehat{E}(k_x,k_y,t = t_\infty)$',fontsize=12)
im = ax[1].imshow((nl['Emn'][:,:,-1]),interpolation="none",cmap="Greys",origin="lower",norm=LogNorm(vmin=1e-6,vmax=1e-1))
divider = make_axes_locatable(ax[1])
cax = divider.append_axes('right', size='5%', pad=0.2)
fig.colorbar(im, cax=cax, orientation='vertical')

ax[2].set_title(r'$\zeta(x,y,t = t_\infty)$',fontsize=12)
im = ax[2].imshow((nl['Vxy'][:,:,-1]),interpolation="bicubic",cmap="seismic",origin="lower",vmin=-2,vmax=2)
divider = make_axes_locatable(ax[2])
cax = divider.append_axes('right', size='5%', pad=0.2)
fig.colorbar(im, cax=cax, orientation='vertical')

for a in ax[0:2]:
    a.set_xticks([0,M-1,2*M-2])
    a.set_yticks([0,N-1,2*N-2])
    a.set_xticklabels([r'$-N_x$',r'$0$',r'$N_x$'],fontsize=10)
    a.set_yticklabels([r'$-N_y$',r'$0$',r'$N_y$'],fontsize=10)

ax[2].set_xticks([0,M-1,2*M-2])
ax[2].set_yticks([0,N-1,2*N-2])
ax[2].set_xticklabels([r'$-\pi$',r'$0$',r'$\pi$'],fontsize=10)
ax[2].set_yticklabels([r'$-\pi$',r'$0$',r'$\pi$'],fontsize=10)

plt.savefig(dn+'jets_nl_spectra.png',bbox_inches='tight',dpi=512)

fig,ax = plt.subplots(figsize=(6,3))

im = ax.imshow(nl['Vyt'],interpolation="bicubic",cmap="seismic",origin="lower",vmin=-2,vmax=2)
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='2%', pad=0.2)
fig.colorbar(im, cax=cax, orientation='vertical')

ax.set_ylabel(r'$\zeta(y,t)$',fontsize=12)
ax.set_xlabel(r'$t$',fontsize=12)

ax.set_xticks([0,25,50,75,100])
ax.set_yticks([0,N-1,2*N-2])
ax.set_xticklabels([r'$0$',r'$250$',r'$500$',r'$750$',r'$1000$'],fontsize=10)
ax.set_yticklabels([r'$-\pi$',r'$0$',r'$\pi$'],fontsize=10)

plt.savefig(dn+'jets_nl_hoevmoeller.png',bbox_inches='tight',dpi=512)
