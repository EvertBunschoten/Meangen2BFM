import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import os

fsize = 15
ticksize = 13
P_t_in = 105325
T_t_in = 288.7
fracfactor = 1.3
tickformat='%.2f'
figsize=(6, 5)

save_direct = os.getcwd()+'/Images/'

z_blade, p_blade, T_blade, pt_blade, Tt_blade, ds_blade, alpha_blade, Mach_blade, flux_blade = np.loadtxt(save_direct+"../axialData_blade.txt", unpack=True, delimiter='\t', skiprows=1)
z_bfm, p_bfm, T_bfm, pt_bfm, Tt_bfm, ds_bfm, alpha_bfm, Mach_bfm, flux_bfm = np.loadtxt(save_direct+"../axialData_BFM.txt", unpack=True, delimiter='\t', skiprows=1)


fig = plt.figure(0, figsize=figsize)
ax = fig.add_subplot(111)
plt.title("Mass flow averaged pressure plot", fontsize=fsize)
ax.plot(z_blade, p_blade/P_t_in, 'r', label="Blade computation")
ax.plot(z_bfm, p_bfm/P_t_in, 'b', label="BFM")
ax.plot(z_blade, 1.01*p_blade/P_t_in, 'k--', label="+- 1%")
ax.plot(z_blade, 0.99*p_blade/P_t_in, 'k--')
plt.xlabel("x[m]", fontsize=fsize)
plt.ylabel(r"$\frac{p}{p_{t, in}}$", fontsize=fracfactor*fsize)
ax.legend(fontsize=fsize, bbox_to_anchor=(1.04, 1), loc="upper left")
ax.grid()
ax.tick_params(axis='both', which='both', labelsize=ticksize)
ax.xaxis.set_major_formatter(FormatStrFormatter(tickformat))
ax.yaxis.set_major_formatter(FormatStrFormatter(tickformat))
fig.savefig(save_direct + 'P_stat.eps', bbox_inches='tight', format='eps')

fig = plt.figure(1, figsize=figsize)
ax = fig.add_subplot(111)
plt.title("Mass flow averaged Mach plot", fontsize=fsize)
ax.plot(z_blade, Mach_blade, 'r', label="Blade computation")
ax.plot(z_bfm, Mach_bfm, 'b', label="BFM")
ax.plot(z_blade, 1.01*Mach_blade, 'k--', label="+- 1%")
ax.plot(z_blade, 0.99*Mach_blade, 'k--')
plt.xlabel("x[m]", fontsize=fsize)
plt.ylabel(r"$M$", fontsize=fsize)
ax.legend(fontsize=fsize, bbox_to_anchor=(1.04, 1), loc="upper left")
ax.grid()
ax.tick_params(axis='both', which='both', labelsize=ticksize)
ax.xaxis.set_major_formatter(FormatStrFormatter(tickformat))
ax.yaxis.set_major_formatter(FormatStrFormatter(tickformat))
fig.savefig(save_direct + 'Mach.eps', bbox_inches='tight', format='eps')


fig = plt.figure(2, figsize=figsize)
ax = fig.add_subplot(111)
plt.title("Mass flux plot", fontsize=fsize)
ax.plot(z_blade, flux_blade, 'r', label="Blade computation")
ax.plot(z_bfm, flux_bfm, 'b', label="BFM")
ax.plot(z_blade, 1.01*flux_blade, 'k--', label="+- 1%")
ax.plot(z_blade, 0.99*flux_blade, 'k--')
plt.xlabel("x[m]", fontsize=fsize)
plt.ylabel(r"$\rho u_x[kg m^{-2}s^{-1}]$", fontsize=fsize)
ax.legend(fontsize=fsize, bbox_to_anchor=(1.04, 1), loc="upper left")
ax.grid()
ax.tick_params(axis='both', which='both', labelsize=ticksize)
ax.xaxis.set_major_formatter(FormatStrFormatter(tickformat))
ax.yaxis.set_major_formatter(FormatStrFormatter(tickformat))
fig.savefig(save_direct + 'flux.eps', bbox_inches='tight', format='eps')

plt.show()
