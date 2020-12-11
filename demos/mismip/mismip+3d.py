import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import tqdm
from mpi4py import MPI
import firedrake
from firedrake import (
    inner, dx, Constant, interpolate, as_vector, exp, sqrt, max_value
)
import icepack, icepack.plot

parser = argparse.ArgumentParser()
parser.add_argument('--input')
parser.add_argument('--output')
parser.add_argument('--level', type=int, default=0)
parser.add_argument('--time', type=float)
parser.add_argument('--timestep', type=float)
args = parser.parse_args()

# Create the 2D mesh and function spaces
Lx, Ly = 640e3, 80e3
ny = 20
nx = int(Lx/Ly) * ny
area = Lx * Ly

coarse_mesh = firedrake.RectangleMesh(nx, ny, Lx, Ly)
mesh_hierarchy = firedrake.MeshHierarchy(coarse_mesh, args.level)
mesh2d = mesh_hierarchy[args.level]
Q2d = firedrake.FunctionSpace(mesh2d, 'CG', 1)
V2d = firedrake.VectorFunctionSpace(mesh2d, 'CG', 1)

# Load in the results from the previous stage
h2d = firedrake.Function(Q2d)
u2d = firedrake.Function(V2d)
input_name = os.path.splitext(args.input)[0]
with firedrake.DumbCheckpoint(input_name, mode=firedrake.FILE_READ) as chk:
    timesteps, indices = chk.get_timesteps()
    chk.set_timestep(timesteps[-1], idx=indices[-1])

    chk.load(h2d, name='h')
    chk.load(u2d, name='u')

# Create the 3D mesh and function spaces and lift the previous stage results
# into 3D
mesh = firedrake.ExtrudedMesh(mesh2d, layers=1)
Q = firedrake.FunctionSpace(mesh, 'CG', 1, vfamily='R', vdegree=0)
V_0 = firedrake.VectorFunctionSpace(
    mesh, 'CG', 1, vfamily='R', vdegree=0, dim=2
)
V = firedrake.VectorFunctionSpace(
    mesh, 'CG', 1, vfamily='GL', vdegree=2, dim=2
)

h0 = icepack.lift3d(h2d, Q)
u0 = firedrake.interpolate(icepack.lift3d(u2d, V_0), V)

# Create the bed topography
x, y, ζ = firedrake.SpatialCoordinate(mesh)

x_c = 300e3
X = x / x_c

B_0 = -150
B_2 = -728.8
B_4 = 343.91
B_6 = -50.57
B_x = B_0 + B_2 * X**2 + B_4 * X**4 + B_6 * X**6

f_c = 4e3
d_c = 500
w_c = 24e3

B_y = d_c * (
    1 / (1 + exp(-2 * (y - Ly / 2 - w_c) / f_c)) +
    1 / (1 + exp(+2 * (y - Ly / 2 + w_c) / f_c))
)

z_deep = -720
z_b = interpolate(max_value(B_x + B_y, z_deep), Q)

# We use units of megapascals-meters-years, whereas in the paper the parameters
# are reported in units of pascals, so these values include some conversion
# factors.
A = Constant(20)
C = Constant(1e-2)

# Create the bed friction dissipation functional
from icepack.constants import (
    ice_density as ρ_I,
    water_density as ρ_W,
    gravity as g,
    weertman_sliding_law as m
)

def friction(**kwargs):
    keys = ('velocity', 'thickness', 'surface', 'friction')
    u, h, s, C = map(kwargs.get, keys)

    p_W = ρ_W * g * max_value(0, -(s - h))
    p_I = ρ_I * g * h
    N = p_I - p_W
    τ_c = N / 2

    u_c = (τ_c / C)**m
    u_b = sqrt(inner(u, u))

    return τ_c * ((u_c**(1 / m + 1) + u_b**(1 / m + 1))**(m / (m + 1)) - u_c)

# Create the model object and extra options
model = icepack.models.HybridModel(friction=friction)
opts = {
    'dirichlet_ids': [1],
    'side_wall_ids': [3, 4],
    'diagnostic_solver_type': 'petsc',
    'diagnostic_solver_parameters': {
        'snes_type': 'newtontr',
        'ksp_type': 'preonly',
        'pc_type': 'lu',
        'pc_factor_mat_solver_type': 'mumps'
    }
}
solver = icepack.solvers.FlowSolver(model, **opts)

# Create the initial data
h = h0.copy(deepcopy=True)
s = icepack.compute_surface(thickness=h, bed=z_b)
a = firedrake.interpolate(firedrake.Constant(0.3), Q)

u = solver.diagnostic_solve(
    velocity=u0,
    thickness=h,
    surface=s,
    fluidity=A,
    friction=C
)

# Run the simulation
δt = args.timestep
num_steps = int(args.time / δt)

output_name = os.path.splitext(args.output)[0]
with firedrake.DumbCheckpoint(output_name, mode=firedrake.FILE_CREATE) as chk:
    trange = tqdm.trange(num_steps)
    for step in trange:
        chk.set_timestep(step * δt)
        chk.store(h, name='h')
        chk.store(u, name='u')

        h = solver.prognostic_solve(
            δt, thickness=h, velocity=u, accumulation=a, thickness_inflow=h0
        )
        s = icepack.compute_surface(thickness=h, bed=z_b)

        u = solver.diagnostic_solve(
            velocity=u,
            thickness=h,
            surface=s,
            fluidity=A,
            friction=C
        )

        hmin = h.dat.data_ro.min()
        min_thickness = firedrake.COMM_WORLD.allreduce(hmin, MPI.MIN)
        avg_thickness = firedrake.assemble(h * dx) / area
        msg = f'avg/min thickness: ({avg_thickness:4.2f}, {min_thickness:4.2f})'
        trange.set_description(msg)

    chk.set_timestep(num_steps * δt)
    chk.store(h, name='h')
    chk.store(u, name='u')
