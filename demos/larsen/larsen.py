import os
import argparse
import subprocess
import numpy as np
import geojson
import rasterio
import firedrake
from firedrake import inner, as_vector, grad, sqrt, dx
import icepack, icepack.models, icepack.datasets, icepack.meshing

parser = argparse.ArgumentParser()
parser.add_argument('--input')
parser.add_argument('--output')
parser.add_argument('--geometry', default='larsen.geo')
parser.add_argument('--mesh', default='larsen.msh')
parser.add_argument('--level', type=int, default=0)
parser.add_argument('--regularization', type=float, default=2.5e3)
parser.add_argument('--thickness-smoothing', type=float, default=8e3)
parser.add_argument('--velocity-smoothing', type=float, default=2.5e3)
parser.add_argument('--gamma', type=float, default=0.1)
args = parser.parse_args()

# Read in the geometry and make a mesh
outline_filename = icepack.datasets.fetch_larsen_outline()
with open(outline_filename, 'r') as outline_file:
    outline = geojson.load(outline_file)

geometry = icepack.meshing.collection_to_geo(outline)
with open(args.geometry, 'w') as geo_file:
    geo_file.write(geometry.get_code())

subprocess.call('gmsh -v 0 -2 -format msh2 -o {:s} {:s}'
                .format(args.mesh, args.geometry).split(' '))
coarse_mesh = firedrake.Mesh(args.mesh)
mesh_hierarchy = firedrake.MeshHierarchy(coarse_mesh, args.level)
mesh = mesh_hierarchy[args.level]

# Read in the observational data
thickness_filename = icepack.datasets.fetch_bedmachine_antarctica()
thickness = rasterio.open('netcdf:' + thickness_filename + ':thickness', 'r')

velocity_filename = icepack.datasets.fetch_measures_antarctica()
vx = rasterio.open('netcdf:' + velocity_filename + ':VX', 'r')
vy = rasterio.open('netcdf:' + velocity_filename + ':VY', 'r')
stdx = rasterio.open('netcdf:' + velocity_filename + ':ERRX', 'r')
stdy = rasterio.open('netcdf:' + velocity_filename + ':ERRX', 'r')

# Make function spaces and interpolate the gridded data
Q = firedrake.FunctionSpace(mesh, family='CG', degree=1)
V = firedrake.VectorFunctionSpace(mesh, family='CG', degree=1)

# The thickness requires some smoothing to avoid pathological values
params = {
    'solver_parameters': {
        'mat_type': 'aij',
        'ksp_type': 'preonly',
        'pc_type': 'lu',
        'pc_factor_mat_solver_type': 'mumps'
    }
}
h_obs = icepack.interpolate(thickness, Q)
h = h_obs.copy(deepcopy=True)
α = firedrake.Constant(args.thickness_smoothing)
J = 0.5 * ((h - h_obs)**2 + α**2 * inner(grad(h), grad(h))) * dx
firedrake.solve(firedrake.derivative(J, h) == 0, h, **params)

u_obs = icepack.interpolate((vx, vy), V)
σx = icepack.interpolate(stdx, Q)
σy = icepack.interpolate(stdy, Q)
v = u_obs.copy(deepcopy=True)
β = firedrake.Constant(args.velocity_smoothing)
J = 0.5 * (inner(v - u_obs, v - u_obs) + β**2 * inner(grad(v), grad(v))) * dx
firedrake.solve(firedrake.derivative(J, v) == 0, v, **params)

# Make a model and initialize things
T0 = 260
A0 = icepack.rate_factor(T0)
from icepack.constants import glen_flow_law as n
def viscosity(u, h, θ):
    A = A0 * firedrake.exp(θ)
    return icepack.models.viscosity.viscosity_depth_averaged(u, h, A)

θ = firedrake.Function(Q)
if args.input:
    with firedrake.DumbCheckpoint(input_name, mode=firedrake.FILE_READ) as chk:
        chk.load(θ, name='θ')

ice_shelf = icepack.models.IceShelf(viscosity=viscosity)
opts = {'dirichlet_ids': [2, 4, 5, 6, 7, 8, 9], 'tol': 1e-6}
u = ice_shelf.diagnostic_solve(u0=v, h=h, θ=θ, **opts)

# Make the misfit and regularization functionals
γ = firedrake.Constant(0.1)
def objective(u):
    δu = as_vector(((u[0] - u_obs[0]) / σx, ((u[1] - u_obs[1]) / σy)))
    return (sqrt(inner(δu, δu) + γ**2) - γ) * dx

L = firedrake.Constant(args.regularization)
def regularization(θ):
    return 0.5 * L**2 * inner(grad(θ), grad(θ)) * dx

# Create an inverse problem and a solver and run the solver
import icepack.inverse
problem = icepack.inverse.InverseProblem(
    model=ice_shelf,
    method=icepack.models.IceShelf.diagnostic_solve,
    objective=objective,
    regularization=regularization,
    state_name='u',
    state=u,
    parameter_name='θ',
    parameter=θ,
    model_args={'h': h, 'u0': u, 'tol': 1e-6},
    dirichlet_ids=opts['dirichlet_ids']
)

area = firedrake.assemble(firedrake.Constant(1) * dx(mesh))
def callback(solver):
    E = firedrake.assemble(solver.objective)
    R = firedrake.assemble(solver.regularization)
    print('{:g}, {:g}'.format(E/area, R/area))

solver = icepack.inverse.GaussNewtonSolver(problem, callback)
iterations = solver.solve(
    rtol=5e-3,
    atol=0.0,
    max_iterations=30
)

# Save the output to an HDF5 file
output_name = os.path.splitext(args.output)[0]
with firedrake.DumbCheckpoint(output_name, mode=firedrake.FILE_CREATE) as chk:
    chk.store(solver.parameter, name='θ')
    chk.store(solver.state, name='u')
