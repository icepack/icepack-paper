import argparse
import matplotlib.pyplot as plt
import numpy as np

m = 3

def τ_schoof(u, u_c, τ_c):
    return τ_c * u**(1/m) / (u + u_c)**(1/m)

def τ_modified(u, u_c, τ_c):
    return τ_c * u**(1/m) / (u**(1/m + 1) + u_c**(1/m + 1))**(1/(m + 1))

τ_c = 1e2
u_c = 2.5e2

u = np.geomspace(u_c / 100, 10 * u_c, 201)

fig, ax = plt.subplots(figsize=(8, 3))

τ_s = τ_schoof(u, u_c, τ_c)
τ_m = τ_modified(u, u_c, τ_c)

ax.plot([u_c, u_c], [0.95 * τ_s[0], 1.05 * τ_c], 'k--', linewidth=0.75, label='critical speed')

ax.plot(u, τ_c / u_c**(1/m) * u**(1/m), 'C1', label='Weertman')
ax.plot(u, τ_s, 'C2', label='Schoof')
ax.plot(u, τ_m, 'C0', label='Modified Schoof')

ax.set_xscale('log')
ax.set_ylim(0.95 * τ_s[0], 1.05 * τ_c)
ax.set_xlabel('speed (m / year)')
ax.set_ylabel('stress (kPa)')
ax.legend()

parser = argparse.ArgumentParser()
parser.add_argument('--output')
args = parser.parse_args()

fig.savefig(args.output, dpi=300, bbox_inches='tight')
