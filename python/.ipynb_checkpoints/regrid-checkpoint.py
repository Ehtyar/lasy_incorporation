import numpy as np
from scipy.interpolate import interp1d

from lasy.utils.grid import Grid

try:
    from tqdm.auto import tqdm
    tqdm_available = True
    bar_format='{l_bar}{bar}| {elapsed}<{remaining} [{rate_fmt}{postfix}]'
except Exception:
    tqdm_available = False


def regrid_laser(laser, Nx=None, Ny=None, Nt=None):
    """rescales the laser to have (Nx,Ny,Nt) points and returns the new laser.
    Only works for xyt laser"""

    grid = laser.grid
    
    dirs = []
    if Nx is not None:
        dirs.append(0)
    else:
        Nx = grid.npoints[0]
    if Ny is not None:
        dirs.append(1)
    else:
        Ny = grid.npoints[1]
    if Nt is not None:
        dirs.append(2)
    else:
        Nt = grid.npoints[2]
    if len(dirs) == 0:
        return laser # nothing to do here
    
    field = grid.get_temporal_field()

    field_new = np.zeros((Nx, Ny, Nt), dtype=field.dtype)
    
    if len(dirs) == 1:
        d = dirs[0]
        Ntr = [Nx, Ny, Nt]
        Nd = Ntr.pop(d)
        axis_old = grid.axes[d]
        axis_new = np.linspace(grid.lo[d], grid.hi[d], Nd)

        if tqdm_available:
            pbar = tqdm(total=Ntr[0], bar_format=bar_format)
        for i0 in range(Ntr[0]):
            for i1 in range(Ntr[1]):
                if d == 0:
                    sl = field[:,i0,i1]
                if d == 1:
                    sl = field[i0,:,i1]
                if d == 2:
                    sl = field[i0,i1,:]
                # interpolate
                interp_fu_abs = interp1d(axis_old, np.abs(sl))
                slice_abs = interp_fu_abs(axis_new)
                interp_fu_angl = interp1d(axis_old, np.unwrap(np.angle(sl)))
                slice_angl = interp_fu_angl(axis_new)
                if d == 0:
                    field_new[:,i0,i1] = slice_abs * np.exp(1j * slice_angl)
                if d == 1:
                    field_new[i0,:,i1] = slice_abs * np.exp(1j * slice_angl)
                if d == 2:
                    field_new[i0,i1,:] = slice_abs * np.exp(1j * slice_angl)
            
            if tqdm_available:
                pbar.update(1)
            else:
                if i0%20 == 19:
                    print(i0+1, "out of", Ntr[0])
        if tqdm_available:
            pbar.close()
    if len(dirs) == 2:
        pass # interpolate
    if len(dirs) == 3:
        raise NotImplementedError("Cannot interpolate all 3 dimensions (yet)")

    grid_new = Grid("xyt", grid.lo, grid.hi, (Nx, Ny, Nt))
    grid_new.set_temporal_field(field_new)
    laser.grid = grid_new
    return laser
        
def cutgrid_laser(laser, Nx=None, Ny=None, Nt=None):
    """Cuts the middle Nx, Ny and Nt from the laser grid and returns the new laser.
    only works for laser in xyt"""
    grid = laser.grid
    field = grid.get_temporal_field()

    if Nx is not None:
        assert Nx <= grid.npoints[0]
        xmin = grid.npoints[0]//2 - Nx//2
        xmax = grid.npoints[0]//2 + Nx//2
        xmid = (grid.hi[0]+grid.lo[0])/2.
        xran = (grid.hi[0]-grid.lo[0])/2.
        xlo = xmid - xran*(Nx / grid.npoints[0])
        xhi = xmid + xran*(Nx / grid.npoints[0])
        field = field[xmin:xmax,:,:]
    else:
        xlo = grid.lo[0]
        xhi = grid.hi[0]
        Nx = grid.npoints[0]
        
    if Ny is not None:
        assert Ny <= grid.npoints[1]
        ymin = grid.npoints[1]//2 - Ny//2
        ymax = grid.npoints[1]//2 + Ny//2
        ymid = (grid.hi[1]+grid.lo[1])/2.
        yran = (grid.hi[1]-grid.lo[1])/2.
        ylo = ymid - yran*(Ny / grid.npoints[1])
        yhi = ymid + yran*(Ny / grid.npoints[1])
        field = field[:,ymin:ymax,:]
    else:
        ylo = grid.lo[1]
        yhi = grid.hi[1]
        Ny = grid.npoints[1]
        
    if Nt is not None:
        assert Nt <= grid.npoints[2]
        tmin = grid.npoints[2]//2 - Nt//2
        tmax = grid.npoints[2]//2 + Nt//2
        tmid = (grid.hi[2]+grid.lo[2])/2.
        tran = (grid.hi[2]-grid.lo[2])/2.
        tlo = tmid - tran*(Nt / grid.npoints[2])
        thi = tmid + tran*(Nt / grid.npoints[2])
        field = field[:,:,tmin:tmax]
    else:
        tlo = grid.lo[2]
        thi = grid.hi[2]
        Nt = grid.npoints[2]

    lo = (xlo, ylo, tlo)
    hi = (xhi, yhi, thi)
    grid_new = Grid(laser.dim, lo, hi, (Nx, Ny, Nt))
    grid_new.set_temporal_field(field)
    laser.grid = grid_new
    return laser