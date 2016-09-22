import glob
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import cPickle as pickle

matplotlib.rcParams['font.family'] = 'Times New Roman'
#matplotlib.rcParams['text.usetex'] = True

plot_vs_total_mass = True  ## True to plot vs m_tot, False to plot vs m_min


allowed_regions = pickle.load(open("allowed_regions.pkl","rb"))

if( plot_vs_total_mass):
    xlim = [1e-2,1]
    idx_to_use = [2,3]
    xticklist = [1e-2, 1e-1]
    suffix = '_mtot'
    xlab = "$m_{\mathrm{tot}}$ [eV]"
else:
    xlim = [1e-4,1]
    idx_to_use = [0,1]
    xticklist = [1e-4, 1e-3, 1e-2, 1e-1]
    suffix = '_mmin'
    xlab = "$m_{\mathrm{min}}$ [eV]"
    
fig,(ax1,ax2) = plt.subplots(1,2,sharey=True)
for f in allowed_regions[idx_to_use[0]]:
    ax1.loglog( f[0], f[1], 'r'+f[2], linewidth=1.5 )

ax1.set_xlim(xlim)
ax1.set_ylim([1e-3,1])
ax1.set_xticks(xticklist)

for f in allowed_regions[idx_to_use[1]]:
    ax2.loglog( f[0], f[1], 'b'+f[2], linewidth=1.5 )
    
fig.subplots_adjust(wspace=0,bottom=0.12,top=0.95,right=0.95)

ax2.set_xlim(xlim)
ax2.set_ylim([1e-3,1])
ax2.set_xticks(xticklist + [1,])

ax1.set_xlabel( xlab ) 
ax2.set_xlabel( xlab ) 
ax1.set_ylabel( r"$<m_{\beta\beta}>$ [eV]" ) 

fig.set_size_inches(6,4.5)
plt.savefig("sens_allowed_only"+suffix+".pdf")

if( plot_vs_total_mass ):
    ax1.fill_between([0.23, 1],[1e-3,1e-3],[1,1], color='k', alpha= 0.15, edgecolor='none')
    ax2.fill_between([0.23, 1],[1e-3,1e-3],[1,1], color='k', alpha= 0.15, edgecolor='none')
    ax2.text(0.4,1.2e-3, "Planck,\narXiv:1502.01589",rotation="vertical",ha='left',va='bottom')
    
exo200_curr = 1.1e25 ## half-life limit from nature paper (yr)
exo200_proj = 5.7e25 ## half-life sensitivity from Phase1+2 projection (yr)
nexo_notag = 6.2e27 ## half-life limit in 5 yrs evaluated for McKeown/2015 (yr)
nexo_batag = 28.3e27 ## half-life limit in 5 yrs evaluated with Ba-tag for McKeown/2015 (yr)

mmin = 1.55 ## minimum matrix element, PRC 87, 064302 (2013)
mmax = 4.20 ## maximum matrix element, PRL 105, 252503 (2010)
Gxe = 1./6.88e24

lmin_curr, lmax_curr = np.sqrt(1./(mmin**2 * Gxe * exo200_curr)), np.sqrt(1./(mmax**2 * Gxe * exo200_curr))
lmin_proj, lmax_proj = np.sqrt(1./(mmin**2 * Gxe * exo200_proj)), np.sqrt(1./(mmax**2 * Gxe * exo200_proj))
lmin_notag, lmax_notag = np.sqrt(1./(mmin**2 * Gxe * nexo_notag)), np.sqrt(1./(mmax**2 * Gxe * nexo_notag))
lmin_batag, lmax_batag = np.sqrt(1./(mmin**2 * Gxe * nexo_batag)), np.sqrt(1./(mmax**2 * Gxe * nexo_batag))

print lmin_curr, lmin_proj, lmin_notag, lmin_batag

## exo200 current
col=np.array([102.,178.,255.])/256.
ax1.fill_between(xlim,[lmax_curr, lmax_curr],[lmin_curr, lmin_curr],edgecolor=col,facecolor=col,alpha=0.5)
ax2.fill_between(xlim,[lmax_curr, lmax_curr],[lmin_curr, lmin_curr],edgecolor=col,facecolor=col,alpha=0.5)
col=np.array([255.,153.,204.])/256.
ax1.fill_between(xlim,[lmax_proj, lmax_proj],[lmin_proj, lmin_proj],edgecolor=col,facecolor=col,alpha=0.5)
ax2.fill_between(xlim,[lmax_proj, lmax_proj],[lmin_proj, lmin_proj],edgecolor=col,facecolor=col,alpha=0.5)

## nexo
col=np.array([0.,150.,150.])/256.
ax1.fill_between(xlim,[lmax_notag, lmax_notag],[lmin_notag, lmin_notag],edgecolor=col,facecolor=col,alpha=0.5)
ax2.fill_between(xlim,[lmax_notag, lmax_notag],[lmin_notag, lmin_notag],edgecolor=col,facecolor=col,alpha=0.5)
col=np.array([180.,50.,150.])/256.
ax1.fill_between(xlim,[lmax_batag, lmax_batag],[lmin_batag, lmin_batag],edgecolor=col,facecolor=col,alpha=0.5)
ax2.fill_between(xlim,[lmax_batag, lmax_batag],[lmin_batag, lmin_batag],edgecolor=col,facecolor=col,alpha=0.5)

#ax1.text( xlim[0]*1.15, 0.85*np.mean([lmin_proj, lmin_curr]), r"EXO-200 current"+"\n"+r"($T_{1/2}=1.1\times10^{25}$ yr)", ha="left", va="center", fontsize=12)
#ax1.text( xlim[0]*1.15, 0.85*np.mean([lmax_proj, lmax_curr]), r"EXO-200 projected"+"\n"+r"($T_{1/2}=5.7\times10^{25}$ yr)", ha="left", va="center", fontsize=12)
ax1.text( xlim[0]*1.15, 0.85*np.mean([lmin_curr, lmax_curr]), r"EXO-200", ha="left", va="center", fontsize=12)
ax1.text( (xlim[0]*23.0,xlim[0]*5.20)[plot_vs_total_mass], 0.83*np.mean([lmin_curr, lmax_curr]), r"Nature 510, 229 (2014)", ha="left", va="center", fontsize=9, fontstyle='italic')
ax1.text( xlim[0]*1.15, 0.85*np.mean([lmin_proj, lmax_proj]), r"EXO-200 Phase-II", ha="left", va="center", fontsize=12)
ax1.text( xlim[0]*1.15, 0.95*np.mean([lmin_notag, lmax_notag]), r"nEXO 5 Years", ha="left", va="center", fontsize=12)
ax1.text( xlim[0]*1.15, 0.95*np.mean([lmax_notag, lmax_batag]), r"nEXO 5 Years w/ Ba-tagging", ha="left", va="center", fontsize=12)

plt.savefig("sens_nexo"+suffix+"_batag_v0.pdf")


plt.show()
