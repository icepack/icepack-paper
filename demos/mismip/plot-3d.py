import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import firedrake
from firedrake import inner, sqrt
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
mesh2d = mesh_hierarchy[args.level]
mesh = firedrake.ExtrudedMesh(mesh2d, 1)

Q = firedrake.FunctionSpace(mesh, 'CG', 1, vfamily='R', vdegree=0)
V = firedrake.VectorFunctionSpace(mesh, 'CG', 1, vfamily='GL', vdegree=2, dim=2)

h = firedrake.Function(Q)
u = firedrake.Function(V)

input_name = os.path.splitext(args.input)[0]
with firedrake.DumbCheckpoint(input_name, mode=firedrake.FILE_READ) as chk:
    timesteps, indices = chk.get_timesteps()
    chk.set_timestep(timesteps[-1], idx=indices[-1])

    chk.load(h, name='h')
    chk.load(u, name='u')

x, y, ζ = firedrake.SpatialCoordinate(mesh)

P_0 = 1
P_1 = np.sqrt(3) * (2 * ζ - 1)
P_2 = np.sqrt(5) * (6 * ζ**2 - 6 * ζ + 1)

weight = P_0 - np.sqrt(3) * P_1 + np.sqrt(5) * P_2
u_b = icepack.depth_average(u, weight=weight)

weight = P_0 + np.sqrt(3) * P_1 + np.sqrt(5) * P_2
u_s = icepack.depth_average(u, weight=weight)

U_b = sqrt(inner(u_b, u_b))
U_s = sqrt(inner(u_s, u_s))
Q2d = firedrake.FunctionSpace(mesh2d, 'CG', 1)
ratio = firedrake.project(U_b / U_s, Q2d)

fig, axes = icepack.plot.subplots()
axes.set_xlim(0, 640e3)
axes.set_ylim(0, 80e3)
axes.get_yaxis().set_visible(False)
colors = icepack.plot.tripcolor(ratio, vmin=0.8, vmax=1.0, axes=axes)
fig.colorbar(colors, ax=axes, fraction=0.0075, pad=0.04)
fig.savefig(args.output, dpi=300, bbox_inches='tight')
