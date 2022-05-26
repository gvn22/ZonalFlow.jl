import numpy as np
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

dn = "./data/"
M,N = 12,12

nl = np.load(dn+"nl_pj.npz",allow_pickle=True)
ql = np.load(dn+"ql_pj.npz",allow_pickle=True)
ce2 = np.load(dn+"ce2_pj.npz",allow_pickle=True)
gql = np.load(dn+"gql_pj.npz",allow_pickle=True)
gce2 = np.load(dn+"gce2_pj.npz",allow_pickle=True)

fig,ax = plt.subplots(1,5,sharey='row',figsize=(15,3))

ax[0].set_title(r'$NL$',fontsize=10)
ax[0].plot(nl['t'],nl['Emt'][0],'k',label='$0$')
for i,x in enumerate(nl['Emt'][1:]):
    ax[0].plot(nl['t'],x,label=r'${i+1}$')

ax[1].set_title(r'$QL$',fontsize=10)
ax[1].plot(ql['t'],ql['Emt'][0],'k',label='$0$')
for i,x in enumerate(ql['Emt'][1:]):
    ax[1].plot(ql['t'],x,label=r'${i+1}$')

ax[2].set_title(r'$CE2$',fontsize=10)
ax[2].plot(ce2['t'],ce2['Emt'][0],'k',label='$0$')
for i,x in enumerate(ce2['Emt'][1:]):
    ax[2].plot(ce2['t'],x,label=r'${i+1}$')

ax[3].set_title(r'$GQL(1)$',fontsize=10)
ax[3].plot(gql['t'],gql['Emt'][0],'k',label='$0$')
for i,x in enumerate(gql['Emt'][1:]):
    ax[3].plot(gql['t'],x,label=r'${i+1}$')

ax[4].set_title(r'$GCE2(1)$',fontsize=10)
ax[4].plot(gce2['t'],gce2['Emt'][0],'k',label='$0$')
for i,x in enumerate(gce2['Emt'][1:]):
    ax[4].plot(gce2['t'],x,label=r'${i+1}$')

for a in ax:
    a.set_xlabel(r'$t$',fontsize=10)
    a.set_yscale('log')
    a.set_ylim(1e-6,2e2)

ax[0].set_ylabel(r'$E_m$',fontsize=10)

plt.subplots_adjust(wspace=0.1)

plt.savefig(dn+'pointjet_zonal_energy.png',bbox_inches='tight',dpi=512)
