import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import firedrake
from firedrake import dx
import icepack.plot

parser = argparse.ArgumentParser()
parser.add_argument('--output')
parser.add_argument('--retreated')
parser.add_argument('--readvanced')
args = parser.parse_args()

Lx, Ly = 640e3, 80e3
ny = 20
nx = int(Lx/Ly) * ny
coarse_mesh = firedrake.RectangleMesh(nx, ny, Lx, Ly)
level = 3
mesh_hierarchy = firedrake.MeshHierarchy(coarse_mesh, level)
mesh = mesh_hierarchy[level]

Q = firedrake.FunctionSpace(mesh, family='CG', degree=1)
h = firedrake.Function(Q)

def get_volumes(name):
    with firedrake.DumbCheckpoint(name, mode=firedrake.FILE_READ) as chk:
        timesteps, indices = chk.get_timesteps()
        volumes = np.zeros_like(timesteps)

        for timestep, index in zip(timesteps, indices.astype(int)):
            chk.set_timestep(timestep, idx=index)
            chk.load(h, name='h')
            volumes[index] = firedrake.assemble(h * dx)

    return volumes

retreating_volumes = get_volumes(os.path.splitext(args.retreated)[0])
readvancing_volumes = get_volumes(os.path.splitext(args.readvanced)[0])
volumes = np.concatenate((retreating_volumes, readvancing_volumes[1:]))

fig, axes = plt.subplots()
axes.plot(1/24 * np.array(list(range(len(volumes)))), volumes / 1e9 / 917.)
axes.set_xlabel('Simulation time (years)')
axes.set_ylabel('Ice mass (gigatons)')
fig.savefig(args.output, dpi=300, bbox_inches='tight')
