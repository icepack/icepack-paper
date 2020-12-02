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

fig, axes = icepack.plot.subplots(ncols=2)

axes[0].set_ylabel('distance (meters)')
for ax in axes:
    ax.get_xaxis().set_visible(False)

xmin, xmax = -50e3, +50e3
ymin, ymax = -175e3, -75e3

axes[0].set_title('Bed elevation')
colors = icepack.plot.tripcolor(b, cmap='magma', axes=axes[0])
fig.colorbar(
    colors, label='meters', ax=axes[0], orientation='horizontal', pad=0.02
)
axes[0].plot(
    [xmin, xmax, xmax, xmin, xmin],
    [ymin, ymin, ymax, ymax, ymin],
    color='k'
)

axes[1].set_title('Velocity')
axes[1].set_xlim(xmin, xmax)
axes[1].set_ylim(ymin, ymax)
axes[1].yaxis.set_label_position('right')
axes[1].yaxis.tick_right()
firedrake.tricontour(b, 20, cmap='Greys', axes=axes[1])
streamlines = firedrake.streamplot(
    u, resolution=5e3, cmap='magma', seed=1, axes=axes[1]
)
fig.colorbar(
    streamlines,
    label='meters / year',
    ax=axes[1],
    orientation='horizontal',
    pad=0.02
)

fig.savefig(args.output, bbox_inches='tight')
