import os
import argparse
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import firedrake
import icepack.plot

parser = argparse.ArgumentParser()
parser.add_argument('--input')
parser.add_argument('--output')
args = parser.parse_args()

mesh = firedrake.UnitDiskMesh(5)
R = 250e3
mesh.coordinates.dat.data[:] *= R

Q = firedrake.FunctionSpace(mesh, family='CG', degree=2)
V = firedrake.VectorFunctionSpace(mesh, family='CG', degree=2)

b = firedrake.Function(Q)
h0 = firedrake.Function(Q)
h = firedrake.Function(Q)
u = firedrake.Function(V)

input_name = os.path.splitext(args.input)[0]
with firedrake.DumbCheckpoint(input_name, mode=firedrake.FILE_READ) as chk:
    chk.load(b, name='bed')
    chk.load(h, name='thickness')
    chk.load(h0, name='thickness-initial')
    chk.load(u, name='velocity')

δh = firedrake.interpolate(h - h0, Q)

fig, axes = icepack.plot.subplots(ncols=2, sharex=True, sharey=True)

axes[0].set_ylabel('distance (meters)')
for ax in axes:
    ax.get_xaxis().set_visible(False)

axes[0].set_title('Bed elevation')
colors = icepack.plot.tripcolor(b, cmap='magma', axes=axes[0])
fig.colorbar(colors, ax=axes[0], orientation='horizontal', pad=0.02)

axes[1].set_title('Thickness change')
colors = icepack.plot.tripcolor(
    δh, vmin=-200, vmax=200, axes=axes[1], cmap='twilight'
)
fig.colorbar(colors, ax=axes[1], orientation='horizontal', pad=0.02)

fig.savefig(args.output, bbox_inches='tight')
