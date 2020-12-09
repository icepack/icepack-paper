import os
import argparse
import firedrake
import icepack.plot

parser = argparse.ArgumentParser()
parser.add_argument('--input')
parser.add_argument('--level', type=int)
parser.add_argument('--output')

args = parser.parse_args()

Lx, Ly = 640e3, 80e3
ny = 20
nx = int(Lx/Ly) * ny
coarse_mesh = firedrake.RectangleMesh(nx, ny, Lx, Ly)
mesh_hierarchy = firedrake.MeshHierarchy(coarse_mesh, args.level)
mesh = mesh_hierarchy[args.level]

Q = firedrake.FunctionSpace(mesh, family='CG', degree=1)
V = firedrake.VectorFunctionSpace(mesh, family='CG', degree=1)

h = firedrake.Function(Q)
u = firedrake.Function(V)

input_name = os.path.splitext(args.input)[0]
with firedrake.DumbCheckpoint(input_name, mode=firedrake.FILE_READ) as chk:
    timesteps, indices = chk.get_timesteps()
    chk.set_timestep(timesteps[-1], idx=indices[-1])

    chk.load(h, name='h')
    chk.load(u, name='u')

fig, axes = icepack.plot.subplots(
    nrows=2, sharex=True, sharey=True, figsize=(6.4, 2.8)
)

axes[0].get_xaxis().set_visible(False)
for ax in axes:
    ax.set_xlim(0, 640e3)
    ax.set_ylim(0, 80e3)
    ax.get_yaxis().set_visible(False)

colors_h = icepack.plot.tripcolor(h, axes=axes[0])
fig.colorbar(colors_h, ax=axes[0], fraction=0.0075, pad=0.04, label='m')
axes[0].set_title('Thickness')

colors_u = icepack.plot.tripcolor(u, axes=axes[1])
fig.colorbar(colors_u, ax=axes[1], fraction=0.0075, pad=0.04, label='m/year')
axes[1].set_title('Velocity')

fig.savefig(args.output, dpi=300, bbox_inches='tight')
