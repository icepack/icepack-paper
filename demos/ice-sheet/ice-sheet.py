import os
import argparse
import tqdm
import numpy as np
import firedrake
from firedrake import (
    max_value,
    min_value,
    Constant,
    interpolate,
    sqrt,
    inner,
    exp
)
import icepack
from icepack.constants import (
    ice_density as ρ_I,
    glen_flow_law as n,
    gravity as g
)

parser = argparse.ArgumentParser()
parser.add_argument('--timestep', type=float)
parser.add_argument('--num-steps', type=int)
parser.add_argument('--output')
args = parser.parse_args()

# Create the mesh and function spaces
mesh = firedrake.UnitDiskMesh(5)
R = 250e3
mesh.coordinates.dat.data[:] *= R

Q = firedrake.FunctionSpace(mesh, family='CG', degree=2)
V = firedrake.VectorFunctionSpace(mesh, family='CG', degree=2)

# Create the bed elevation
x, y = firedrake.SpatialCoordinate(mesh)
r = sqrt(x**2 + y**2)

b_base = firedrake.Constant(400)
b_max = firedrake.Constant(1400)

# Radius of the plateau interior and ridge
ro = 125e3
Ro = firedrake.Constant(200e3)

def tanh(z):
    return (exp(z) - exp(-z)) / (exp(z) + exp(-z))

def θ(z):
    return (tanh(z) + 1) / 2

def sech(z):
    return 2 / (exp(z) + exp(-z))

a = firedrake.Constant(50e3)
ξ = (sqrt(x**2 + y**2) - ro) / a

b_expr_plateau = b_base * (1 - θ(3 * ξ))

ζ = (r - Ro) / Ro
b_expr_ridge = (b_max - b_base) * sech(3 * ξ)

ρ1 = firedrake.Constant(1/4)
μ1 = 1 - ρ1 * θ(3 * (x - ro/4) / a) * sech(2 * y / a)

ρ2 = firedrake.Constant(3/8)
μ2 = 1 - ρ2 * θ(3 * (y - ro/4) / a) * sech(2 * x / a)

ρ3 = firedrake.Constant(1/2)
μ3 = 1 - ρ3 * θ(3 * (-x + ro/4) / a) * sech(2 * y / a)

ρ4 = firedrake.Constant(5/8)
μ4 = 1 - ρ4 * θ(3 * (-y + ro/4) / a) * sech(2 * x / a)

μ = μ1 * μ2 * μ3 * μ4
S = 480 / (1 - Ro / R)

b_expr_valleys = (b_max - b_base) * sech(3 * ξ) * μ - θ(5 * ζ) * S * ζ

b_expr = b_expr_plateau + b_expr_valleys
b = interpolate(b_expr, Q)

# Create the surface elevation
max_radius = 195e3
dome_height = 2.4e3
dome = dome_height * max_value(1 - (x**2 + y**2) / max_radius**2, 0)
s0 = interpolate(dome, Q)
h0 = interpolate(max_value(s0 - b, 0), Q)

# Create the model and solver objects and initialize the velocity
model = icepack.models.ShallowIce()
solver = icepack.solvers.FlowSolver(model)

T = firedrake.Constant(273.15 - 5)
A = interpolate(icepack.rate_factor(T), Q)

h = h0.copy(deepcopy=True)
s = s0.copy(deepcopy=True)
u = solver.diagnostic_solve(
    velocity=firedrake.Function(V),
    thickness=h,
    surface=s0,
    fluidity=A
)

# Create the mass balance
ela = 300.
max_a = 0.  # alternative: 0.5
da_ds = 0.  # alternative: 0.5 / 1000

def mass_balance(s, max_a, da_ds, ela):
    return min_value((s - ela) * da_ds, max_a)

# Pick a timestep, a duration, and run the simulation
dt = args.timestep
num_timesteps = args.num_steps

for step in tqdm.trange(num_timesteps):
    a = interpolate(mass_balance(s, ela, max_a, da_ds), Q)
    h = solver.prognostic_solve(
        dt,
        thickness=h,
        accumulation=a,
        velocity=u
    )

    h.interpolate(max_value(h, 0))
    s = interpolate(h + b, Q)

    u = solver.diagnostic_solve(
        velocity=u,
        thickness=h,
        surface=s,
        fluidity=A
    )

# Write out the results to a file
output_name = os.path.splitext(args.output)[0]
with firedrake.DumbCheckpoint(output_name, mode=firedrake.FILE_CREATE) as chk:
    chk.store(b, name='bed')
    chk.store(h, name='thickness')
    chk.store(h0, name='thickness-initial')
    chk.store(u, name='velocity')
