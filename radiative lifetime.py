import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const
from matplotlib.pyplot import rcParams
from matplotlib.ticker import FixedLocator, FixedFormatter
import matplotlib.patches as mpatches

# rootgrp = nc.Dataset("./dipole/ndb.dipoles_fragment_2", "r", format="NETCDF4")
#
# print(rootgrp.dimensions)
# print(rootgrp.variables)
#
# #
# dipoles_P = np.array(rootgrp['DIP_P_k_0002_spin_0001'][:])
#
# print(dipoles_P[8,10,0,0])
# np.savetxt("dipoles_P",dipoles_P[:,:,2,0],fmt="%.8f")
#
mat = 'pbi'
num_excitons = 99
num_k = 332
val = 9

bohr_to_m = const.physical_constants['Bohr radius'][0]
m = const.physical_constants['electron mass'][0]
e = const.physical_constants['elementary charge'][0]
c = const.physical_constants['speed of light in vacuum'][0]
kB = const.physical_constants['Boltzmann constant'][0]
eps_0 = const.physical_constants['vacuum electric permittivity'][0]
h_to_j = const.physical_constants['hartree-joule relationship'][0]

square = 8.8 * 7.62102 * (bohr_to_m ** 2)
exciton_mass = 2.28 * m
GW_gap = 4
Nk = 60 * 60

speed_of_light_atom_unit = const.physical_constants['speed of light in vacuum'][0] / const.physical_constants['atomic unit of velocity'][0]
E = np.loadtxt('./o-2D_WR_WC.exc_qpt1_E_sorted') * e
weights = []
dipoles = []
for i in range(1,num_excitons+1):
    weights.append(np.loadtxt('./weights/o-2D_WR_WC.exc_qpt1_weights_at_'+str(i)))

for i in range(1,num_k):
    with nc.Dataset("./dipole/ndb.dipoles_fragment_%d"%i, "r", format="NETCDF4") as f:
        dipoles.append(np.array(f['DIP_P_k_%04d_spin_0001'%i][:]))
## dipole matrix elements in x and y direction

velocity_y_im = np.zeros(num_excitons)
velocity_y_re = np.zeros(num_excitons)
velocity_x_im = np.zeros(num_excitons)
velocity_x_re = np.zeros(num_excitons)
scattering = np.zeros((2,num_excitons))
lifetime = np.zeros((2,num_excitons))
for i in range(num_excitons):
    for j in weights[i]:
        integrand_x_re = j[-2] * dipoles[int(j[2])-1][int(j[0])-1,int(j[1])-val-1,0,0]
        integrand_x_im = j[-2] * dipoles[int(j[2])-1][int(j[0])-1,int(j[1])-val-1,0,1]
        velocity_x_re[i] = velocity_x_re[i] + integrand_x_re
        velocity_x_im[i] = velocity_x_im[i] + integrand_x_im
        integrand_y_re = j[-2] * dipoles[int(j[2])-1][int(j[0])-1,int(j[1])-val-1,1,0]
        integrand_y_im = j[-2] * dipoles[int(j[2])-1][int(j[0])-1,int(j[1])-val-1,1,1]
        velocity_y_re[i] = velocity_x_re[i] + integrand_x_re
        velocity_y_im[i] = velocity_x_im[i] + integrand_x_im

for i in range(99):
    scattering[0, i] = (e ** 2) * (velocity_x_re[i] ** 2 + velocity_x_im[i] ** 2) * h_to_j * m /((eps_0 * (m ** 2) * c * square * Nk * E[i,0]))
    lifetime[0, i] = (10 ** 12)  /scattering[0, i]

    scattering[1, i] = (e ** 2) * (velocity_x_re[i] ** 2 + velocity_x_im[i] ** 2) * h_to_j * m /((eps_0 * (m ** 2) * c * square * Nk * E[i,0]))
    lifetime[1, i] = (10 ** 12)  /scattering[1, i]

plt.rc('font',family='Times New Roman')
fig,ax=plt.subplots(figsize=[6,6])

ax.spines['bottom'].set_linewidth(2);###设置底部坐标轴的粗细
ax.spines['left'].set_linewidth(2);####设置左边坐标轴的粗细
ax.spines['right'].set_linewidth(2);###设置右边坐标轴的粗细
ax.spines['top'].set_linewidth(2);####设置上部坐标轴的粗细
ax.tick_params(axis="x", direction="out", width=2, length=4)
ax.tick_params(axis="y", direction="out", width=2, length=4)
ax.set_xlim(3.4,4.4)
ax.set_ylim(0,2.5)
ax.set_xticks([3.4,3.6,3.8,4.0,4.2,4.4])
ax.set_xticklabels([3.4,3.6,3.8,4.0,4.2,4.4],size=22,position=(0,-0.02))
ax.set_yticks([0,0.5,1.0,1.5,2.0,2.5])
ax.set_yticklabels([0,0.5,1.0,1.5,2.0,2.5],size=22,position=(-0.02,0))

ax.axvline(x=4.0,ls="-",c="black",alpha=0.6)

# ax.set_xlabel("Photon Energy (eV)",size=22)
ax.set_ylabel("%s (ps)"%chr(964),size=22,position=(0.2,0.5))

ax.axvline(x=GW_gap,ls="-",c="black",alpha=0.5)
a = plt.scatter(E[:99,0], lifetime[0,:], 100, edgecolors='dodgerblue', alpha=0.5, marker="o")
b = plt.scatter(E[:99,0], lifetime[1,:], 100, edgecolors='r', alpha=0.5, marker="o")
plt.tight_layout()
plt.savefig(mat+'.pdf',format='pdf', dpi=2000)

np.savetxt("E",E.T)
np.savetxt("lifetime",lifetime.T)