import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
import firedrake
from firedrake import sqrt, inner
import icepack.plot

def colorbar(figure, axes, mappable, *args, **kwargs):
    divider = make_axes_locatable(axes)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    return figure.colorbar(mappable, *args, cax=cax, **kwargs)

mesh = firedrake.Mesh('ice-shelf.msh')
Q = firedrake.FunctionSpace(mesh, family='CG', degree=2)
V = firedrake.VectorFunctionSpace(mesh, family='CG', degree=2)
Δ = firedrake.FunctionSpace(mesh, family='DG', degree=1)

filename = 'steady-state-undamaged'
with firedrake.DumbCheckpoint(filename, mode=firedrake.FILE_READ) as chk:
    h_undamaged = firedrake.Function(Q)
    u_undamaged = firedrake.Function(V)
    chk.load(h_undamaged, name='h')
    chk.load(u_undamaged, name='u')

filename = 'steady-state-damaged'
with firedrake.DumbCheckpoint(filename, mode=firedrake.FILE_READ) as chk:
    h_damaged = firedrake.Function(Q)
    u_damaged = firedrake.Function(V)
    D = firedrake.Function(Δ)
    chk.load(h_damaged, name='h')
    chk.load(u_damaged, name='u')
    chk.load(D, name='D')

δh = firedrake.interpolate(h_damaged - h_undamaged, Q)

speed = sqrt(inner(u_undamaged, u_undamaged))
v = u_undamaged / speed
δu = firedrake.interpolate(inner(u_damaged - u_undamaged, v), Q)

fig, axes = icepack.plot.subplots(nrows=1, ncols=2, sharex=True, sharey=True)
axes[0].set_ylabel('distance (m)')

umax = icepack.norm(u_undamaged, norm_type='Linfty')
umax = (int(umax // 50) + 1) * 50
streamlines = icepack.plot.streamplot(u_undamaged, precision=1e3, density=5e3, axes=axes[0])
colorbar(fig, axes[0], streamlines, label='velocity (m/a)')

dumax = icepack.norm(δu, norm_type='Linfty')
dumax = int(dumax) + 1
levels = np.linspace(-dumax, +dumax, 2 * dumax + 1)
contours = icepack.plot.tricontourf(δu, levels=levels, cmap='seismic', axes=axes[1])
colorbar(fig, axes[1], contours, label='$\Delta$ velocity (m/a)')

fig.tight_layout(h_pad=1)
fig.savefig('velocity.png', dpi=300, bbox_inches='tight')


fig, axes = icepack.plot.subplots(nrows=1, ncols=2, sharex=True, sharey=True)
axes[0].set_ylabel('distance (m)')

hmax = icepack.norm(h_damaged, norm_type='Linfty')
hmax = (int(hmax) // 50) * 50
levels = np.linspace(0, hmax, 11)
contours = icepack.plot.tricontourf(h_damaged, levels=levels, axes=axes[0])
colorbar(fig, axes[0], contours, label='thickness (m)')

Dmax = icepack.norm(D, norm_type='Linfty')
triangles = icepack.plot.tripcolor(D, axes=axes[1])
contours = icepack.plot.tricontour(h_damaged, 5, colors='black', alpha=.5, axes=axes[1])
colorbar(fig, axes[1], triangles, label='damage')

fig.tight_layout(h_pad=1)
fig.savefig('thickness.png', dpi=300, bbox_inches='tight')
