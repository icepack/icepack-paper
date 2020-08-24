import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import geojson
import rasterio
import firedrake
import icepack.datasets, icepack.plot

parser = argparse.ArgumentParser()
parser.add_argument('--mesh')
parser.add_argument('--input')
parser.add_argument('--output')
args = parser.parse_args()

outline_filename = icepack.datasets.fetch_larsen_outline()
with open(outline_filename, 'r') as outline_file:
    outline = geojson.load(outline_file)

xmin, ymin, xmax, ymax = np.inf, np.inf, -np.inf, -np.inf
δ = 50e3
for feature in outline['features']:
    for line_string in feature['geometry']['coordinates']:
        xs = np.array(line_string)
        x, y = xs[:, 0], xs[:, 1]
        xmin, ymin = min(xmin, x.min() - δ), min(ymin, y.min() - δ)
        xmax, ymax = max(xmax, x.max() + δ), max(ymax, y.max() + δ)

image_filename = icepack.datasets.fetch_mosaic_of_antarctica()
with rasterio.open(image_filename, 'r') as image_file:
    height, width = image_file.height, image_file.width
    transform = image_file.transform
    window = rasterio.windows.from_bounds(
        left=xmin, bottom=ymin, right=xmax, top=ymax,
        width=width, height=height, transform=transform
    )
    image = image_file.read(indexes=1, window=window, masked=True)

def subplots(*args, **kwargs):
    fig, axes = icepack.plot.subplots(*args, **kwargs)
    xmin, ymin, xmax, ymax = rasterio.windows.bounds(window, transform)

    try:
        axes.imshow(image, extent=(xmin, xmax, ymin, ymax),
                    cmap='Greys_r', vmin=12e3, vmax=16.38e3)
    except AttributeError:
        for ax in axes:
            ax.imshow(image, extent=(xmin, xmax, ymin, ymax),
                      cmap='Greys_r', vmin=12e3, vmax=16.38e3)

    return fig, axes

mesh = firedrake.Mesh(args.mesh)

Q = firedrake.FunctionSpace(mesh, family='CG', degree=2)
V = firedrake.VectorFunctionSpace(mesh, family='CG', degree=2)

θ = firedrake.Function(Q)
u = firedrake.Function(V)

input_name = os.path.splitext(args.input)[0]
with firedrake.DumbCheckpoint(input_name, mode=firedrake.FILE_READ) as chk:
    chk.load(θ, name='log_fluidity')
    chk.load(u, name='velocity')

fig, axes = subplots(ncols=2, sharex=True, sharey=True)

axes[0].set_xlabel('x (meters)')
levels = np.linspace(-10, +10, 21)
contours = icepack.plot.tricontourf(
    θ, axes=axes[0], levels=levels, extend='both', alpha=0.8
)
fig.colorbar(contours, ax=axes[0], fraction=0.046, pad=0.05)
axes[0].set_title('Logarithm of inferred ice fluidity', y=1.08)

streamlines = icepack.plot.streamplot(
    u, axes=axes[1], precision=1e3, density=2.5e3
)
fig.colorbar(
    streamlines, ax=axes[1], label='meters / year', fraction=0.046, pad=0.05
)
axes[1].set_title('Computed ice velocity', y=1.08)

fig.savefig(args.output, dpi=300, bbox_inches='tight')
