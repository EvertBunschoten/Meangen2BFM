import numpy as np
import matplotlib.pyplot as plt

fsize = 16
ticksize = 13

z, rho, f_p, tds_dz, dp_dz, dmom_dz = np.loadtxt("outputData.txt", unpack=True, delimiter='\t', skiprows=1)
figsize = (6, 5)

fig = plt.figure(1, figsize=figsize)
ax = fig.add_subplot(111)
plt.title("Parallel force vs entropy gradient", fontsize=fsize)
ax.plot(z, tds_dz, 'r', label=r"$T\frac{ds}{dx}$")
ax.plot(z, f_p, 'b', label=r"$\frac{1}{\rho}F_p$")
plt.xlabel(r"x[m]", fontsize=fsize)
plt.ylabel(r"$T\frac{ds}{dx}[J kg^{-1} m^{-1}]$", fontsize=fsize)
ax.tick_params(axis='both', which='both', labelsize=ticksize)
ax.legend(fontsize=fsize, bbox_to_anchor=(1.04, 1), loc="upper left")
ax.grid()
fig.savefig('dS_dz.eps', bbox_inches='tight', format='eps')
err = []
Z = []
for i in range(len(z)):
    if z[i] > 0.01 and z[i] < 0.04:
        err.append(np.sqrt(((f_p[i] / tds_dz[i]) - 1)**2))
        Z.append(z[i])

fig = plt.figure(2)
plt.plot(Z, err)
plt.show()   
err = 100 * sum(err) / len(err)

print("Entropy error: "+str(err))
fig = plt.figure(3, figsize=figsize)
ax = fig.add_subplot(111)
plt.title("Parallel force vs pressure gradient", fontsize=fsize)
ax.plot(z, dp_dz, 'r', label=r"$\frac{dp}{dx}$")
ax.plot(z, -rho * f_p , 'b', label=r"$-F_p$")
#ax.plot(z, -rho * f_p -dmom_dz , 'k--', label=r"$-F_p - \frac{d(\rho u_x^2)}{dx}$")
plt.xlabel(r"x[m]", fontsize=fsize)
plt.ylabel(r"$\frac{dp}{dx} [kg m^{-2} s^{-2}]$", fontsize=fsize)
ax.tick_params(axis='both', which='both', labelsize=ticksize)
ax.legend(fontsize=fsize, bbox_to_anchor=(1.04, 1), loc="upper left")
ax.grid()
fig.savefig('dp_dz.eps', bbox_inches='tight', format='eps')

plt.show()

err = []
Z = []
for i in range(len(z)):
    if z[i] > 0.01 and z[i] < 0.04:
        err.append(np.sqrt((((-rho[i] * f_p[i] - dmom_dz[i]) / dp_dz[i]) - 1)**2))
        Z.append(z[i])
fig = plt.figure(3, figsize=figsize)
plt.plot(Z, err)
plt.show()  
err = 100 * sum(err) / len(err)
print("Pressure gradient error: "+str(err))
data = np.loadtxt("history.csv", delimiter=',', skiprows=1)
fig = plt.figure(4)
ax = fig.add_subplot(111)
#plt.title("Flat plate convergence history", fontsize=fsize)
ax.plot(data[:, 0], data[:, 13], 'k', label=r"$\rho$ residual")
plt.xlabel(r"Iteration[-]", fontsize=fsize)
plt.ylabel(r"Residual value[-]", fontsize=fsize)
ax.tick_params(axis='both', which='both', labelsize=ticksize)
ax.legend(fontsize=fsize, bbox_to_anchor=(1.04, 1), loc="upper left")
ax.grid()
fig.savefig('convergence.eps', bbox_inches='tight', format='eps')

plt.show()