import matplotlib.pyplot as plt
import numpy as np

MAX_P = 20.
MAX_Z = 1.
N_P = 200
N_Z = 200
H_P = 2. * MAX_P / N_P
INV_H_P = 1. / H_P
H_Z = 2. * MAX_Z / N_Z
INV_H_Z = 1. / H_Z
TAU = 0.0004
COLL_NU = 100.
G = 10.

p_s = np.linspace(-MAX_P + .5*H_P, MAX_P - .5*H_P, N_P)
z_s = np.linspace(-MAX_Z + .5*H_Z, MAX_Z - .5*H_Z, N_Z)

fig_z_s, fig_p_s = np.mgrid[
    slice(-MAX_Z, MAX_Z + H_Z, H_Z),
    slice(-MAX_P, MAX_P + H_P, H_P)]

# problem_1
data = np.fromfile('output/1.bin', dtype=np.float64).reshape((N_Z, N_P))
plt.pcolormesh(fig_z_s, fig_p_s, data, cmap='viridis', vmin=0, vmax=0.13)
plt.tight_layout(1.08,None,None,(0.02,0.02,1.05,1))
plt.xlabel(r'$z$'); plt.ylabel(r'$p_z$')
plt.colorbar(label='$f$')
# plt.savefig(f'figures/1.eps'); plt.clf()
plt.show()

# problem_2
data = np.fromfile('output/2_H.bin', dtype=np.float64)
ts = np.arange(0, len(data)) * TAU
plt.plot(ts, data)
plt.tight_layout(1.08,None,None,(0.02,0.02,1,1))
plt.xlabel(r'$t$'); plt.ylabel(r'$H$')
# plt.savefig(f'figures/2_H.eps'); plt.clf()
plt.show()

for i in range(15):
    t = i * TAU * 10
    data = np.fromfile('output/2_phase.bin', dtype=np.float64, offset=8*i*N_Z*N_P, count=N_Z*N_P).reshape((N_Z, N_P))
    plt.plot(z_s, np.sum(data, axis=1), label=f'{t=:0<6.3}')
plt.tight_layout(1.08,None,None,(0.02,0.02,1,0.98))
plt.xlabel(r'$z$'); plt.ylabel(r'$f$')
plt.legend()
# plt.savefig(f'figures/2_f_z.eps'); plt.clf()
plt.show()

for i in range(15):
    t = i * TAU * 10
    data = np.fromfile('output/2_phase.bin', dtype=np.float64, offset=8*i*N_Z*N_P, count=N_Z*N_P).reshape((N_Z, N_P))
    plt.plot(p_s, np.sum(data, axis=0), label=f'{t=:0<6.3}')
plt.tight_layout(1.08,None,None,(0.02,0.02,1,0.98))
plt.xlabel(r'$p$'); plt.ylabel(r'$f$')
plt.legend()
# plt.savefig(f'figures/2_f_p.eps'); plt.clf()
plt.show()

for i in range(8):
    t = i * TAU * 10
    data = np.fromfile('output/2_phase.bin', dtype=np.float64, offset=8*i*N_Z*N_P, count=N_Z*N_P).reshape((N_Z, N_P))
    plt.pcolormesh(fig_z_s, fig_p_s, data, cmap='viridis', vmin=0)
    plt.tight_layout(1.08,None,None,(0.02,0.02,1.05,0.98))
    plt.xlabel(r'$z$'); plt.ylabel(r'$p_z$')
    plt.title(f'$t = {t:0<6.3}$')
    plt.colorbar(label='$f$')
    # plt.savefig(f'figures/2_{t:0<6.3}.eps'); plt.clf()
    plt.show()

# problem_3
data = np.fromfile('output/3_H.bin', dtype=np.float64)
ts = np.arange(0, len(data)) * TAU
plt.plot(ts, data)
plt.tight_layout(1.08,None,None,(0.02,0.02,1,1))
plt.xlabel(r'$t$'); plt.ylabel(r'$H$')
# plt.savefig(f'figures/3_H.eps'); plt.clf()
plt.show()

for i in range(15):
    t = i * TAU * 10
    data = np.fromfile('output/3_phase.bin', dtype=np.float64, offset=8*i*N_Z*N_P, count=N_Z*N_P).reshape((N_Z, N_P))
    plt.plot(z_s, np.sum(data, axis=1), label=f'{t=:0<6.3}')
plt.tight_layout(1.08,None,None,(0.02,0.02,1,0.98))
plt.xlabel(r'$z$'); plt.ylabel(r'$f$')
plt.legend()
# plt.savefig(f'figures/3_f_z.eps'); plt.clf()
plt.show()

for i in range(15):
    t = i * TAU * 10
    data = np.fromfile('output/3_phase.bin', dtype=np.float64, offset=8*i*N_Z*N_P, count=N_Z*N_P).reshape((N_Z, N_P))
    plt.plot(p_s, np.sum(data, axis=0), label=f'{t=:0<6.3}')
plt.tight_layout(1.08,None,None,(0.02,0.02,1,0.98))
plt.xlabel(r'$p$'); plt.ylabel(r'$f$')
plt.legend()
# plt.savefig(f'figures/3_f_p.eps'); plt.clf()
plt.show()


for i in range(8):
    t = i * TAU * 10
    data = np.fromfile('output/3_phase.bin', dtype=np.float64, offset=8*i*N_Z*N_P, count=N_Z*N_P).reshape((N_Z, N_P))
    plt.pcolormesh(fig_z_s, fig_p_s, data, cmap='viridis', vmin=0)
    plt.tight_layout(1.08,None,None,(0.02,0.02,1.05,0.98))
    plt.xlabel(r'$z$'); plt.ylabel(r'$p_z$')
    plt.title(f'$t = {t:0<6.3}$')
    plt.colorbar(label='$f$')
    # plt.savefig(f'figures/3_{t:0<6.3}.eps'); plt.clf()
    plt.show()
