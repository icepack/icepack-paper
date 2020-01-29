import os
import argparse
import tqdm
import firedrake
from firedrake import as_vector, inner, ds
import icepack

parser = argparse.ArgumentParser()
parser.add_argument('--mesh')
parser.add_argument('--input')
parser.add_argument('--output')
parser.add_argument('--damage', action='store_true')
parser.add_argument('--final-time', type=float)
parser.add_argument('--num-steps', type=int)
args = parser.parse_args()

# Load the mesh and create some function spaces
mesh = firedrake.Mesh(args.mesh)
Q = firedrake.FunctionSpace(mesh, family='CG', degree=2)
V = firedrake.VectorFunctionSpace(mesh, family='CG', degree=2)
Δ = firedrake.FunctionSpace(mesh, family='DG', degree=1)

# Load the initial data from an HDF5 file if available
if args.input:
    h0 = firedrake.Function(Q)
    u0 = firedrake.Function(V)

    input_name = os.path.splitext(args.input)[0]
    with firedrake.DumbCheckpoint(input_name, mode=firedrake.FILE_READ) as chk:
        chk.load(h0, name='h')
        chk.load(u0, name='u')

# Otherwise, create some synthetic initial data symbolically
else:
    import numpy as np
    from numpy import pi as π

    inlet_angles = π * np.array([-3/4, -1/2, -1/3, -1/6])
    inlet_widths = π * np.array([1/8, 1/12, 1/24, 1/12])

    x = firedrake.SpatialCoordinate(mesh)

    R = 200e3
    u_in = 300
    h_in = 350
    hb = 100
    dh, du = 400, 250

    hs, us = [], []
    for θ, ϕ in zip(inlet_angles, inlet_widths):
        x0 = R * as_vector((np.cos(θ), np.sin(θ)))
        v = -as_vector((np.cos(θ), np.sin(θ)))
        L = inner(x - x0, v)
        W = x - x0 - L * v
        Rn = 2 * ϕ / π * R
        q = firedrake.max_value(1 - (W / Rn)**2, 0)
        hs.append(hb + q * ((h_in - hb) - dh * L /R))
        us.append(firedrake.exp(-4 * (W/R)**2) * (u_in + du * L / R) * v)

    h_expr = firedrake.Constant(hb)
    for h in hs:
        h_expr = firedrake.max_value(h, h_expr)

    u_expr = sum(us)

    # Interpolate the initial data to finite element spaces
    h0 = firedrake.interpolate(h_expr, Q)
    u0 = firedrake.interpolate(u_expr, V)


# Create a model object and propagate the initial state forward by 400 years
from icepack.constants import glen_flow_law as n
from icepack.models.viscosity import viscosity_depth_averaged
def viscosity_damaged(**kwargs):
    u, h, D, A = map(kwargs.get, ('velocity', 'thickness', 'damage', 'fluidity'))
    return viscosity_depth_averaged(
        velocity=u, thickness=h, fluidity=(1 - D)**(-n) * A
    )

viscosity = viscosity_damaged if args.damage else viscosity_depth_averaged
flow_model = icepack.models.IceShelf(viscosity=viscosity)
flow_solver = icepack.solvers.FlowSolver(flow_model, dirichlet_ids=[1])

# Create the fluidity and damage fields
T = firedrake.Constant(255.15)
A = firedrake.interpolate(icepack.rate_factor(T), Q)

D = firedrake.Function(Δ)
D_inflow = firedrake.Constant(0.)
damage_model = icepack.models.DamageTransport()
damage_solver = icepack.solvers.DamageSolver(damage_model)

h = h0.copy(deepcopy=True)
u = flow_solver.diagnostic_solve(
    velocity=u0, thickness=h, fluidity=A, damage=D
)

final_time = args.final_time
num_timesteps = args.num_steps
dt = final_time / num_timesteps
a = firedrake.Constant(0.0)

for step in tqdm.trange(num_timesteps):
    h = flow_solver.prognostic_solve(
        dt, thickness=h, accumulation=a, velocity=u, thickness_inflow=h0
    )
    if args.damage:
        D = damage_solver.solve(
            dt, damage=D, velocity=u, fluidity=A, damage_inflow=D_inflow
        )

    u = flow_solver.diagnostic_solve(
        velocity=u, thickness=h, fluidity=A, damage=D
    )

output_name = os.path.splitext(args.output)[0]
with firedrake.DumbCheckpoint(output_name, mode=firedrake.FILE_CREATE) as chk:
    chk.store(D, name='D')
    chk.store(h, name='h')
    chk.store(u, name='u')
