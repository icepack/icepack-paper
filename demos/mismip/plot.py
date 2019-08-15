import os
import argparse
import firedrake
import icepack.plot

parser = argparse.ArgumentParser()
parser.add_argument('--input')
parser.add_argument('--output')

args = parser.parse_args()

Lx, Ly = 640e3, 80e3
ny = 40
nx = int(Lx/Ly) * ny
mesh = firedrake.RectangleMesh(nx, ny, Lx, Ly)

Q = firedrake.FunctionSpace(mesh, family='CG', degree=1)
V = firedrake.VectorFunctionSpace(mesh, family='CG', degree=1)

h = firedrake.Function(Q)
u = firedrake.Function(V)

input_name = os.path.splitext(args.input)[0]
chk = firedrake.DumbCheckpoint(input_name, mode=firedrake.FILE_READ)
chk.load(h, name='h')
chk.load(u, name='u')

fig, axes = icepack.plot.subplots(nrows=2, sharex=True, sharey=True,
                                  figsize=(6.4, 2.8))
contours_h = icepack.plot.tricontourf(h, 40, axes=axes[0])
fig.colorbar(contours_h, ax=axes[0], fraction=0.0075, pad=0.04, label='m')
axes[0].set_title('Thickness')
contours_u = icepack.plot.tricontourf(u, 40, axes=axes[1])
fig.colorbar(contours_u, ax=axes[1], fraction=0.0075, pad=0.04, label='m/year')
axes[1].set_title('Velocity')
fig.savefig(args.output, dpi=300, bbox_inches='tight')
