import os
import argparse
import subprocess
import numpy as np
import rasterio
import geojson
import firedrake
import icepack, icepack.plot

parser = argparse.ArgumentParser()
parser.add_argument('--geometry', default='larsen.geo')
parser.add_argument('--mesh', default='larsen.msh')
parser.add_argument('--output', default='results.h5')
parser.add_argument('--thickness-smoothing', type=float, default=2e3)
parser.add_argument('--regularization', type=float, default=7.5e3)
args = parser.parse_args()

# Read in the geometry and make a mesh
outline_filename = icepack.datasets.fetch_larsen_outline()
with open(outline_filename, 'r') as outline_file:
    outline = geojson.load(outline_file)
print(f'Loaded shapefile {outline_filename}')

geometry = icepack.meshing.collection_to_geo(outline)
with open(args.geometry, 'w') as geo_file:
    geo_file.write(geometry.get_code())
print(f'Generated geometry file {args.geometry}')

subprocess.run(
    f'gmsh -2 -format msh2 -v 2 -o {args.mesh} {args.geometry}'.split()
)
mesh = firedrake.Mesh(args.mesh)
print(f'Generated and loaded mesh file {args.mesh}')

# Read in the thickness data and smooth it a bit for stability
thickness_filename = icepack.datasets.fetch_bedmachine_antarctica()
thickness = rasterio.open(f'netcdf:{thickness_filename}:thickness', 'r')

Q = firedrake.FunctionSpace(mesh, family='CG', degree=2)
h_obs = icepack.interpolate(thickness, Q)

from firedrake import inner, grad, dx
h = h_obs.copy(deepcopy=True)
α_h = firedrake.Constant(args.thickness_smoothing)
J = 0.5 * ((h - h_obs)**2 + α_h**2 * inner(grad(h), grad(h))) * dx
F = firedrake.derivative(J, h)
parameters = {
    'solver_parameters': {
        'ksp_type': 'preonly',
        'pc_type': 'lu',
        'pc_factor_mat_solve_type': 'mumps'
    }
}
firedrake.solve(F == 0, h, **parameters)
print(f'Loaded and interpolated thickness measurements {thickness_filename}')

# Read in the velocity data and the measurement errors
velocity_filename = icepack.datasets.fetch_measures_antarctica()
vx = rasterio.open(f'netcdf:{velocity_filename}:VX', 'r')
vy = rasterio.open(f'netcdf:{velocity_filename}:VY', 'r')
stdx = rasterio.open(f'netcdf:{velocity_filename}:ERRX', 'r')
stdy = rasterio.open(f'netcdf:{velocity_filename}:ERRY', 'r')

V = firedrake.VectorFunctionSpace(mesh, family='CG', degree=2)
u_obs = icepack.interpolate((vx, vy), V)
σx = icepack.interpolate(stdx, Q)
σy = icepack.interpolate(stdy, Q)
print(f'Loaded and interpolated velocity measurements {velocity_filename}')

# Make a model and initialize things
T = firedrake.Constant(260)
A0 = icepack.rate_factor(T)
from icepack.constants import glen_flow_law as n
def viscosity(**kwargs):
    u = kwargs['velocity']
    h = kwargs['thickness']
    θ = kwargs['log_fluidity']

    A = A0 * firedrake.exp(θ)
    return icepack.models.viscosity.viscosity_depth_averaged(
        velocity=u, thickness=h, fluidity=A
    )

θ = firedrake.Function(Q)

model = icepack.models.IceShelf(viscosity=viscosity)
opts = {'dirichlet_ids': [2, 4, 5, 6, 7, 8, 9]}
solver = icepack.solvers.FlowSolver(model, **opts)

u = solver.diagnostic_solve(
    velocity=u_obs,
    thickness=h,
    log_fluidity=θ
)
print('Initialized log fluidity and modeled velocity')

# Make the misfit and regularization functionals
import icepack.inverse

def objective(u):
    δu = u - u_obs
    return 0.5 * ((δu[0] / σx)**2 + (δu[1] / σy)**2) * dx

Θ = firedrake.Constant(1.)
L = firedrake.Constant(7.5e3)
def regularization(θ):
    return 0.5 * (L / Θ)**2 * inner(grad(θ), grad(θ)) * dx

# Create an inverse problem and a solver and run the solver
problem = icepack.inverse.InverseProblem(
    model=model,
    objective=objective,
    regularization=regularization,
    state_name='velocity',
    state=u,
    parameter_name='log_fluidity',
    parameter=θ,
    solver_kwargs=opts,
    diagnostic_solve_kwargs={'thickness': h}
)

qdegree = model.quadrature_degree(velocity=u, thickness=h, log_fluidity=θ)
fc_params = {'form_compiler_parameters': {'quadrature_degree': qdegree}}
area = firedrake.assemble(firedrake.Constant(1) * dx(mesh))
def callback(solver):
    q = solver.search_direction
    dJ = solver.gradient
    dJ_dq = firedrake.action(dJ, q)
    decrement = firedrake.assemble(dJ_dq, **fc_params) / area

    error = firedrake.assemble(solver.objective) / area
    penalty = firedrake.assemble(solver.regularization) / area
    print(f'{error:10.4g} | {penalty:10.4g} | {decrement:10.4g}')

print('Running inverse solver...\n')
print(' error     | penalty    | decrement')
print('-----------|------------|-----------')
solver = icepack.inverse.GaussNewtonSolver(
    problem, callback, search_max_iterations=200
)

iterations = solver.solve(
    rtol=5e-3,
    etol=1e-6,
    atol=0.0,
    max_iterations=30
)

# Save the output to an HDF5 file
θ = solver.parameter
output_name = os.path.splitext(args.output)[0]
with firedrake.DumbCheckpoint(output_name, mode=firedrake.FILE_CREATE) as chk:
    chk.store(solver.parameter, name='log_fluidity')
    chk.store(solver.state, name='velocity')
