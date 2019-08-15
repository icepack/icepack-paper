import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import firedrake
from firedrake import (inner, dx, Constant, interpolate, as_vector, exp, sqrt,
                       max_value)
import icepack, icepack.models, icepack.plot

parser = argparse.ArgumentParser()
parser.add_argument('--verbose', action='store_true')
parser.add_argument('--input')
parser.add_argument('--output', default='result')
parser.add_argument('--melt', action='store_true')
parser.add_argument('--time', type=float, default=4e3)
parser.add_argument('--timestep', type=float, default=1.0)

args = parser.parse_args()

# Create the mesh and function spaces
Lx, Ly = 640e3, 80e3
ny = 40
nx = int(Lx/Ly) * ny
mesh = firedrake.RectangleMesh(nx, ny, Lx, Ly)
area = firedrake.assemble(Constant(1) * dx(mesh))

Q = firedrake.FunctionSpace(mesh, family='CG', degree=1)
V = firedrake.VectorFunctionSpace(mesh, family='CG', degree=1)

# Create the bed topography
x, y = firedrake.SpatialCoordinate(mesh)

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

B_y = d_c * (1 / (1 + exp(-2 * (y - Ly/2 - w_c) / f_c)) +
             1 / (1 + exp(2 * (y - Ly / 2 + w_c) / f_c)))

z_deep = -720
z_b = interpolate(max_value(B_x + B_y, z_deep), Q)

# We use units of megapascals-meters-years, whereas in the paper the parameters
# are reported in units of pascals, so these values include some conversion
# factors.
A = Constant(20)
C = Constant(1e-2)

# Create the bed friction dissipation functional
from icepack.constants import (ice_density as ρ_I, water_density as ρ_W,
                               gravity as g, weertman_sliding_law as m)
def friction(u, h, s, C):
    u_b = sqrt(inner(u, u))
    p_W = ρ_W * g * max_value(0, -(s - h))
    p_I = ρ_I * g * h
    N = p_I - p_W
    ϕ = N / p_I
    return m / (m + 1) * C * ϕ * u_b**(1/m + 1) * dx

model = icepack.models.IceStream(friction=friction)
opts = {'dirichlet_ids': [1], 'side_wall_ids': [3, 4], 'tol': 1e-8}

# Read in the input data if there is any
if args.input:
    h0 = firedrake.Function(Q)
    u0 = firedrake.Function(V)
    input_name = os.path.splitext(args.input)[0]
    chk = firedrake.DumbCheckpoint(input_name, mode=firedrake.FILE_READ)
    chk.load(h0, name='h')
    chk.load(u0, name='u')
# Otherwise create some rough initial data
else:
    h0 = interpolate(Constant(100), Q)
    s0 = model.compute_surface(h=h0, b=z_b)
    u0 = model.diagnostic_solve(u0=interpolate(as_vector((35 * x / Lx, 0)), V),
                                h=h0, s=s0, A=A, C=C, **opts)

h = h0.copy(deepcopy=True)
s = model.compute_surface(h=h, b=z_b)
u = u0.copy(deepcopy=True)

# Accumulation and melt rate
accumulation = Constant(0.3)

def tanh(z):
    return (exp(z) - exp(-z)) / (exp(z) + exp(-z))

Ω = Constant(0.2 if args.melt else 0.0)
z_0 = Constant(-100)    # cut-off base elevation
h_c0 = Constant(75.0)   # cut-off water column thickness

# Run the simulation
δt = args.timestep
num_steps = int(args.time / δt)
for step in range(num_steps):
    z_d = s - h             # elevation of the ice base
    h_c = z_d - z_b         # water column thickness
    melt = Ω * tanh(h_c / h_c0) * max_value(z_0 - z_d, 0)
    a = accumulation - melt

    h = model.prognostic_solve(δt, h0=h, u=u, a=a, h_inflow=h0)
    s = model.compute_surface(h=h, b=z_b)
    u = model.diagnostic_solve(u0=u, h=h, s=s, A=A, C=C, **opts)

    if args.verbose:
        time = step * δt
        avg_thickness = firedrake.assemble(h * dx) / area
        min_thickness = h.dat.data_ro.min()
        max_thickness = h.dat.data_ro.max()
        print(time, avg_thickness, min_thickness, max_thickness, flush=True)

# Write out the results
output_name = os.path.splitext(args.output)[0]
chk = firedrake.DumbCheckpoint(output_name, mode=firedrake.FILE_CREATE)
chk.store(h, name='h')
chk.store(u, name='u')
