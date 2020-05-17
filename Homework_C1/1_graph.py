import matplotlib.pyplot as plt
import numpy as np

etas = np.linspace(-3,3,400)
for xi, A in zip(
    (0.5, 0.5, 0.8),
    (0.2, 0.4, 0.2)
):
    plt.plot(etas[etas < 0], xi + A / etas[etas < 0], 'C0')
    plt.plot(etas[etas > 0], xi + A / etas[etas > 0], 'C0')
    plt.plot(etas, 1 / np.sqrt(1 + etas**2), 'C1')
    plt.ylim(bottom=-0.5, top=1.5)
    plt.tight_layout(1.08,None,None,(0,0.02,1,1))
    plt.xlabel(r'$\bar{\eta}$')
    # plt.savefig(f'figures/{xi}_{A}.eps'); plt.clf()
    plt.show()

xi_s = np.linspace(-1, 1, 400)
plt.plot(xi_s, (1 - np.abs(xi_s)**(2/3))**1.5)
plt.plot(xi_s, -(1 - np.abs(xi_s)**(2/3))**1.5, 'C0')
plt.plot(xi_s, 0 * xi_s, 'C0')
plt.tight_layout(1.08,None,None,(0.04,0.04,1,1))
plt.xlabel(r'$\xi = \frac{x}{l}$')
plt.ylabel(r'$A = \frac{mg}{kl}$')
# plt.savefig(f'figures/x_g.eps'); plt.clf()
plt.show()

for N in (4, 16, 64, 256):
    plt.plot(np.fromfile(f"output/adiabatic_1.2_y{N}.bin", dtype=np.float64),
             np.fromfile(f"output/adiabatic_1.2_p{N}.bin", dtype=np.float64))
    plt.tight_layout(1.08,None,None,(0.05,0.04,1,1))
    plt.xlabel(r'$\eta = \frac{y}{l}$')
    plt.ylabel(r'$\mathcal{P} = \frac{p_y}{\sqrt{mk}l}$')
    # plt.savefig(f'figures/adiabatic_1.2_yp{N}.eps'); plt.clf()
    plt.show()
    plt.plot(np.fromfile(f"output/adiabatic_1.2_x{N}.bin", dtype=np.float64),
             np.fromfile(f"output/adiabatic_1.2_J{N}.bin", dtype=np.float64))
    plt.tight_layout(1.08,None,None,(0.05,0.04,1,1))
    plt.ylim(bottom=0)
    plt.xlabel(r'$\xi = \frac{x}{l} = \frac{2l-vt}{l}$')
    plt.ylabel(r'$\frac{J}{\sqrt{mk}l^2}$')
    # plt.savefig(f'figures/adiabatic_1.2_xJ{N}.eps'); plt.clf()
    plt.show()

for N in (4, 16, 64, 256):
    plt.plot(np.fromfile(f"output/adiabatic_1.3_y{N}.bin", dtype=np.float64),
             np.fromfile(f"output/adiabatic_1.3_p{N}.bin", dtype=np.float64))
    plt.tight_layout(1.08,None,None,(0.05,0.04,1,1))
    plt.xlabel(r'$\eta = \frac{y}{l}$')
    plt.ylabel(r'$\mathcal{P} = \frac{p_y}{\sqrt{mk}l}$')
    # plt.savefig(f'figures/adiabatic_1.3_yp{N}.eps'); plt.clf()
    plt.show()
    plt.plot(np.fromfile(f"output/adiabatic_1.3_g{N}.bin", dtype=np.float64),
             np.fromfile(f"output/adiabatic_1.3_J{N}.bin", dtype=np.float64))
    plt.tight_layout(1.08,None,None,(0.05,0.04,1,1))
    plt.ylim(bottom=0)
    plt.xlabel(r'$\frac{mg}{kl}$')
    plt.ylabel(r'$\frac{J}{\sqrt{mk}l^2}$')
    # plt.savefig(f'figures/adiabatic_1.3_gJ{N}.eps'); plt.clf()
    plt.show()

for i in range(1, 7):
    N = 256 * i
    plt.plot(np.fromfile(f"output/adiabatic_1.4_y{N}.bin", dtype=np.float64),
             np.fromfile(f"output/adiabatic_1.4_p{N}.bin", dtype=np.float64))
    plt.tight_layout(1.08,None,None,(0.05,0.04,1,1))
    plt.xlabel(r'$\eta = \frac{y}{l}$')
    plt.ylabel(r'$\mathcal{P} = \frac{p_y}{\sqrt{mk}l}$')
    # plt.savefig(f'figures/adiabatic_1.4_yp{N}.eps'); plt.clf()
    plt.show()
    plt.plot(np.fromfile(f"output/adiabatic_1.4_t{N}.bin", dtype=np.float64),
             np.fromfile(f"output/adiabatic_1.4_omega{N}.bin", dtype=np.float64))
    plt.tight_layout(1.08,None,None,(0.05,0.04,1,1))
    plt.xlabel(r'$\tau = \sqrt{\frac{k}{m}}t$')
    plt.ylabel(r'$\sqrt{\frac{m}{k}}\omega$')
    # plt.savefig(f'figures/adiabatic_1.4_t_omega{N}.eps'); plt.clf()
    plt.show()

for i in range(1, 7):
    N = 256 * i
    plt.plot(np.fromfile(f"output/adiabatic_1.5_y{N}.bin", dtype=np.float64),
             np.fromfile(f"output/adiabatic_1.5_p{N}.bin", dtype=np.float64))
    plt.tight_layout(1.08,None,None,(0.05,0.04,1,1))
    plt.xlabel(r'$\eta = \frac{y}{l}$')
    plt.ylabel(r'$\mathcal{P} = \frac{p_y}{\sqrt{mk}l}$')
    # plt.savefig(f'figures/adiabatic_1.5_yp{N}.eps'); plt.clf()
    plt.show()
    plt.plot(np.fromfile(f"output/adiabatic_1.5_t{N}.bin", dtype=np.float64),
             np.fromfile(f"output/adiabatic_1.5_omega{N}.bin", dtype=np.float64))
    plt.tight_layout(1.08,None,None,(0.05,0.04,1,1))
    plt.xlabel(r'$\tau = \sqrt{\frac{k}{m}}t$')
    plt.ylabel(r'$\sqrt{\frac{m}{k}}\omega$')
    # plt.savefig(f'figures/adiabatic_1.5_t_omega{N}.eps'); plt.clf()
    plt.show()
