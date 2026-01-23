"""read openPMD-files and display their metadata and display their data on a symlog-plot."""

import matplotlib.pyplot as plt
import matplotlib.colors as clr
import numpy as np
import openpmd_api as io
from .ptime import (start_clock, set_clock, read_clock, get_clocks)
from scipy.constants import c, mu_0, epsilon_0, e, m_e
from scipy.optimize import curve_fit

try:
    from tqdm.auto import tqdm
    tqdm_available = True
    bar_format='{l_bar}{bar}| {elapsed}<{remaining} [{rate_fmt}{postfix}]'
except Exception:
    tqdm_available = False

from .full_field import show


def _gauss(pos, mux, muy, w, A):
    """for fitting w"""
    x, y = pos
    xx, yy = np.meshgrid(x, y)
    g = A * np.exp(-((xx-mux)**2 + (yy-muy)**2) / w**2)
    return g.ravel()

def _get_density(filename=None, series=None, iteration=None):
    """calculating the energy desity (if possible) or an approximation for it form the given series. Also returns the spacing.
    """
    close = False
    if series is None:
        assert filename is not None, "filename or series mus be given."
        close = True
        series = io.Series(filename, io.Access.read_only)
    else:
        assert filename is None, "can not work with both series and filename being given."
    if iteration is None:
        iteration is min(series.iterations)
    # get all the data
    if "E" in series.iterations[iteration].meshes and "B" in series.iterations[iteration].meshes:
        E_x = series.iterations[iteration].meshes["E"]["x"][:,:,:]
        E_y = series.iterations[iteration].meshes["E"]["y"][:,:,:]
        E_z = series.iterations[iteration].meshes["E"]["z"][:,:,:]
        series.flush()
        E2 = E_x**2 + E_y**2 + E_z**2
        E2 *= series.iterations[iteration].meshes["E"]["x"].get_attribute("unitSI")**2

        B_x = series.iterations[iteration].meshes["B"]["x"][:,:,:]
        B_y = series.iterations[iteration].meshes["B"]["y"][:,:,:]
        B_z = series.iterations[iteration].meshes["B"]["z"][:,:,:]
        series.flush()
        B2 = B_x**2 + B_y**2 + B_z**2
        B2 *= series.iterations[iteration].meshes["B"]["x"].get_attribute("unitSI")**2

        # actual energy density of the electromagnetic field
        density = E2 * epsilon_0/2 + B2 * 1/(2*mu_0)

        spacing = [s * series.iterations[iteration].meshes["E"].get_attribute("gridUnitSI")
                    for s in series.iterations[iteration].meshes["E"].get_attribute("gridSpacing")]
    elif "E" in series.iterations[iteration].meshes:
        print("Warning: without the B-field the w-caculation may be inaccurate")
        E_x = series.iterations[iteration].meshes["E"]["x"][:,:,:]
        E_y = series.iterations[iteration].meshes["E"]["y"][:,:,:]
        E_z = series.iterations[iteration].meshes["E"]["z"][:,:,:]
        series.flush()

        # standin for the energy density of the electromagnetic field
        density = E_x**2 + E_y**2 + E_z**2
        density *= series.iterations[iteration].meshes["E"]["x"].get_attribute("unitSI")**2

        spacing = [s * series.iterations[iteration].meshes["E"].get_attribute("gridUnitSI")
                    for s in series.iterations[iteration].meshes["E"].get_attribute("gridSpacing")]
    elif "laserEnvelope" in series.iterations[iteration].meshes:
        env = series.iterations[iteration].meshes["laserEnvelope"][:,:,:]
        series.flush()

        # standin for the energy density of the electromagnetic field
        density = np.abs(env)**2

        spacing = [s * series.iterations[iteration].meshes["laserEnvelope"].get_attribute("gridUnitSI")
                    for s in series.iterations[iteration].meshes["laserEnvelope"].get_attribute("gridSpacing")]
    else:
        raise ValueError("file needs to contain either E and B fields or the laserEnvelope to calculate w")
    if close:
        series.close()
    return density, spacing

def show_w(filename=None, series=None, iteration=None, forward=0, density=None, spacing=None, method="both", p=True):
    """Calculates the value of w for the iteration in the series in two ways: first by fittig a 
    gauss function and then by statistical standard deviation. Prints the values and returns them.

    Parameters:
    filename : str
        The full path at which the file can be found.
        Either filename, series or density must be given.
        
    series : openpmd_api Series object
        The series for which w should be calculated
        Either filename, series or density must be given.

    density : 3D array of float
        The energy density of the laser field. requires spacing.
        Either filename, series or density must be given.

    spacing : tuple or list of 3 float
        Only relevant if density is given. Describes the spacing of the points in density.

    iteration : int (optional)
        the iteration in the series for which w should be calculated.
        if None the first iteration in the series will be taken.

    forward : int in {0,1,2} (optional)
        the forward direction in the mesh(es)

    method : str in {"fit", "stat", "both"} (optional)
        defines the method by which w should be calculated. "fit" fits a gauss function to the data,
        "stat" calculates the statistical deviation and "both" does both methods in this order.

    p : bool (optional)
        Wether to print the results

    Returns:
    w1 : float
        calculated by a gauss fit

    w2 : float
        calculated via standard deviation
    """
    assert method in ["stat", "fit", "both"], "the given method is not defined."
    if density is None:
        if spacing is not None:
            print("ignoring given spacing as density must be taken from file")
        density, spacing = _get_density(filename=filename, series=series, iteration=iteration)
    else:
        assert spacing is not None, "if density is given spacing is required."
        assert len(spacing) == len(density.shape)
        assert filename is None, "can not work with both density and filename being given."
        assert series is None, "can not work with both density and series being given."
    fl = np.sum(density, axis=forward)
    if np.max(fl) <= 0:
        if p:
            print("Cant calculate w, all the fields are 0")
        if method == "both":
            return 0, 0
        else:
            return 0
    dirs = [0, 1, 2]
    dirs.remove(forward)
    
    if method in ["fit", "both"]:
        if p:
            print("fitting a gauss-function:")
        p0 = [0, 0, spacing[dirs[0]]*fl.shape[0]/4, np.max(np.abs(fl))]
        x1 = np.linspace(-spacing[dirs[0]]*fl.shape[0]/2, spacing[dirs[0]]*fl.shape[0]/2, fl.shape[0])
        x2 = np.linspace(-spacing[dirs[1]]*fl.shape[1]/2, spacing[dirs[1]]*fl.shape[1]/2, fl.shape[1])
        coeff, _ = curve_fit(_gauss, (x1, x2), fl.ravel(), p0=p0)
        if p:
            print("w1 =", np.sqrt(2) * np.abs(coeff[2]))
            print("center at", coeff[0], coeff[1])
        if method == "fit":
            return np.sqrt(2) * np.abs(coeff[2])
            
    if method in ["stat", "both"]:
        if p:
            print("calculating statistical deviation:")
        r2 = np.zeros(fl.shape[0]*fl.shape[1])
        fl_r = np.zeros(fl.shape[0]*fl.shape[1])
        for n1 in range(fl.shape[0]):
            for n2 in range(fl.shape[1]):
                r2[n1*fl.shape[1]+n2] = (n1-fl.shape[0]/2)**2*spacing[dirs[0]]**2+(n2-fl.shape[1]/2)**2*spacing[dirs[1]]**2
                fl_r[n1*fl.shape[1]+n2] = fl[n1, n2]
        w = np.sqrt(2 * np.average(r2, weights=fl_r))
        if p:
            print("w2 =", w)
        if method == "stat":
            return w
            
    return np.sqrt(2) * np.abs(coeff[2]), w

def plot_w(filename=None, series=None, forward=0, ret=False, startit=0, maxtime=np.inf, apply_maxtime_to_measurement=True, comp_func=None, title="", unit=1, unitname="m", ret_ax=False):
    """calculates the values of w for all the iterations in the series and plots them against time.
    
    Parameters:
    filename : str
        The full path at which the file can be found.
        Either filename or series must be given.
        
    series : openpmd_api Series object
        The series for which w should be calculated
        Either filename or series must be given.

    forward : int or str (optional)
        The forward direction of the pulse.
        Can either be the index of the field as it is in the file or a direction from x, y and z.

    ret : bool (optional)
        Wether to return the found arrays after plotting

    startit : int (optional)
        At which iteration to start the calculation.
        Meant for PIC-simulations where in the first iterations the pulse is still entering the moving window.

    maxtime : float (optional)
        A maximum time of the iterations beyond which the values are neither calculated nor plotted.
        Assumes, that the iterations are in time order in the file.

    apply_maxtime_to_measurement : bool (optional)
        Whether to apply maxtime to the measured data as well as the theoretical values.

    comp_func : function or list of functions (optional)
        if not None it will plot both the w values and this function/these functions.
        The functions docstring will be used as label in the plot.

    title : str (optional)
        The title of the plot

    unit : float (optional)
        plots the values of w divided by this number. Set unitname for axis to be correct.
        If ret is True the returned values are still in SI units (i.e. as if unit=1).

    unitname : str(optional)
        The name that should be displayed as the unit of w.

    ret_ax : bool (optional)
        Whether to return the axes oject and the fig object instead of showing it.

    Returns:
    fig : pyplot figure object (optional)
        only returned if ret_ax is True.

    ax : pyplot axes object (optional)
        only returned if ret_ax is True.

    ts : array of float (optional)
        The times of the iterations. Only returned if ret is True.
        
    ws : array of float (optional)
        The values of w in the iterations. Only returned if ret is True.
    """
    close = False
    if series is None:
        assert filename is not None, "filename or series must be given."
        close = True
        series = io.Series(filename, io.Access.read_only)
    else:
        assert filename is None, "can not work with both series and filename being given."

    if forward in ["x", "y", "z"]:
        n_1 = min(series.iterations)
        if "E" in series.iterations[n_1].meshes:
            labels = series.iterations[n_1].meshes["E"].get_attribute("axisLabels")
        elif "laserEnvelope" in series.iterations[n_1].meshes:
            labels = series.iterations[n_1].meshes["laserEnvelope"].get_attribute("axisLabels")
        else:
            raise ValueError("cannot calculate w: no mesh E or laserEnvelope found")
        forward = labels.index(forward)
    else:
        assert forward in [0, 1, 2], "forward must be a direction or valid index."
        
    print("calculating w")
    N = len(series.iterations)
    ws = np.zeros(N)
    ts = np.zeros(N)
    Nmax = N
    Nmin = 0
    if tqdm_available:
        pbar = tqdm(total=N, bar_format=bar_format)
    for n, n_i in enumerate(series.iterations):
        if n_i < startit:
            Nmin = n+1
        else:
            ts[n] = series.iterations[n_i].get_attribute("time") * series.iterations[n_i].get_attribute("timeUnitSI")
            if (ts[n] >= maxtime) and (Nmax == N):
                Nmax = n
                if apply_maxtime_to_measurement:
                    print("maxtime reached")
                    break
            ws[n] = show_w(series=series, iteration=n_i, forward=forward, method="stat", p=False)
        if tqdm_available:
            pbar.update(1)
        else:
            print(n+1, "out of", N, "complete")

    if tqdm_available:
        pbar.close()
    if close:
        series.close()
    
    print("plotting w")
    fig = plt.figure()
    ax = fig.add_subplot()
    if comp_func is not None:
        try:
            comp_func[0]
            mul=True
        except (TypeError, IndexError):
            mul=False
        if mul:
            for func in comp_func:
                ax.plot(ts[Nmin:Nmax], func(ts[Nmin:Nmax])/unit, label=func.__doc__)
        else:
            ax.plot(ts[Nmin:Nmax], comp_func(ts[Nmin:Nmax])/unit, label=comp_func.__doc__)
    ax.set_xlabel("$t$/s")
    ax.set_ylabel("$w$/"+unitname)
    if apply_maxtime_to_measurement:
        ax.plot(ts[Nmin:Nmax], ws[Nmin:Nmax]/unit, ".", label="measured w")
    else:
        ax.plot(ts[Nmin:], ws[Nmin:]/unit, ".", label="measured w")
    ax.set_title(title)
    if comp_func is not None:
        ax.legend()

    if ret_ax:
        if ret:
            return fig, ax, ts, ws
        else:
            return fig, ax
            
    fig.show()
    plt.show()
    
    if ret:
        return ts, ws

def show_lpeak(filename=None, series=None, iteration=None, forward=0, density=None, spacing=None, method="stat", p=True):
    """Calculates the distance form the center of the peak of the pulse for the iteration in the series.

    Parameters:
    filename : str
        The full path at which the file can be found.
        Either filename, series or density must be given.
        
    series : openpmd_api Series object
        The series for which l_peak should be calculated
        Either filename, series or density must be given.

    density : 3D array of float
        The energy density of the laser field. requires spacing.
        Either filename, series or density must be given.

    spacing : tuple or list of 3 float
        Only relevant if density is given. Describes the spacing of the points in density.

    iteration : int (optional)
        the iteration in the series for which l_peak should be calculated.
        if None the first iteration in the series will be taken.

    forward : int in {0,1,2} (optional)
        the forward direction in the mesh(es)

    method : str (optional)
        defines the method by which w should be calculated.
        "stat": calculates the statistical average of the intensity on the axis.
        "avstat": calculates the statistical average of the intensity.
        "fit": fits a gaussian onto the intensity on axis.
        "avfit": fits a gaussian onto the intensity.
        "max": finds the maximum value in the laser field and returns its time.
        "avmax": like "max" but averages the transversal directionas first.

    p : bool (optional)
        Wether to print the results

    Returns:
    l_peak : float
        the distance from the center of the data the peak is at according to the given method
    """
    assert method in ["stat", "avstat", "fit", "avfit", "max", "avmax"], "the given method does not exist."
    if density is None:
        if spacing is not None:
            print("ignoring given spacing as density must be taken from file")
        density, spacing = _get_density(filename=filename, series=series, iteration=iteration)
    else:
        assert spacing is not None, "if density is given spacing is required."
        assert len(spacing) == len(density.shape)
        assert filename is None, "can not work with both density and filename being given."
        assert series is None, "can not work with both density and series being given."
    if np.max(density) <= 0:
        if p:
            print("Cant calculate l_peak, all the fields are 0")
        return 0
    dirs = [0, 1, 2]
    dirs.remove(forward)

    axis = np.linspace(-density.shape[forward]/2 * spacing[forward], density.shape[forward]/2 * spacing[forward], density.shape[forward])
    s_int = np.sum(density, axis=tuple(dirs))
    Nx, Ny = density.shape[dirs[0]]//2, density.shape[dirs[1]]//2
    match forward:
        case 0:
            lineout = density[:,Nx,Ny]
        case 1:
            lineout = density[Nx,:,Ny]
        case 2:
            lineout = density[Nx,Ny,:]
        case _:
            raise ValueError("forward needs to be either 0, 1 or 2")
    match method:
        case "stat":
            l_peak = np.average(axis, weights=lineout)
        case "avstat":
            l_peak = np.average(axis, weights=s_int)
        case "fit":
            p0 = (0,(axis[-1] - axis[0])/4,np.max(lineout))
            coeff, _ = curve_fit(_gauss, axis, lineout, p0=p0)
            l_peak = coeff[0]
        case "avfit":
            p0 = (0,(axis[-1] - axis[0])/4,np.max(s_int))
            coeff, _ = curve_fit(_gauss, axis, s_int, p0=p0)
            l_peak = coeff[0]
        case "max":
            l_peak = axis[np.argmax(lineout)]
        case "avmax":
            l_peak = axis[np.argmax(s_int)]
        case _:
            print("?")
            raise ValueError
    if p:
        print("Peak at", l_peak, "m")
    return l_peak

def plot_lpeak(filename=None, series=None, forward=0, ret=False, startit=0, maxtime=np.inf, apply_maxtime_to_measurement=True, comp_func=None, title="", ret_ax=False):
    """calculates the values of lpeak for all the iterations in the series and plots them against time.
    
    Parameters:
    filename : str
        The full path at which the file can be found.
        Either filename or series must be given.
        
    series : openpmd_api Series object
        The series for which w should be calculated
        Either filename or series must be given.

    forward : int or str (optional)
        The forward direction of the pulse.
        Can either be the index of the field as it is in the file or a direction from x, y and z.

    ret : bool (optional)
        Wether to return the found arrays after plotting

    startit : int (optional)
        The iteration at which the calculation should be started.
        Meant for PIC-simulations where the pulse needs to move into the moving window first.

    maxtime : float (optional)
        A maximum time of the iterations beyond which the values are neither calculated nor plotted.
        Assumes, that the iterations are in time order in the file.

    apply_maxtime_to_measurement : bool (optional)
        Whether to apply maxtime to the measured data as well as the theoretical values.

    comp_func : function (optional)
        if not None it will plot both the lpeak values and this function.
        The functions docstring will be used as label in the plot.

    title : str (optional)
        The title of the plot

    ret_ax : bool (optional)
        Whether to return the axes oject and the fig object instead of showing it.

    Returns:
    fig : pyplot figure object (optional)
        only returned if ret_ax is True.

    ax : pyplot axes object (optional)
        only returned if ret_ax is True.

    ts : array of float (optional)
        The times of the iterations. Only returned if ret is True.
        
    lpeakss : array of float (optional)
        The values of lpeak in the iterations. Only returned if ret is True.
    """
    close = False
    if series is None:
        assert filename is not None, "filename or series must be given."
        close = True
        series = io.Series(filename, io.Access.read_only)
    else:
        assert filename is None, "can not work with both series and filename being given."

    if forward in ["x", "y", "z"]:
        n_1 = min(series.iterations)
        if "E" in series.iterations[n_1].meshes:
            labels = series.iterations[n_1].meshes["E"].get_attribute("axisLabels")
        elif "laserEnvelope" in series.iterations[n_1].meshes:
            labels = series.iterations[n_1].meshes["laserEnvelope"].get_attribute("axisLabels")
        else:
            raise ValueError("cannot calculate lpeak: no mesh E or laserEnvelope found")
        forward = labels.index(forward)
    else:
        assert forward in [0, 1, 2], "forward must be a direction or valid index."
        
    print("calculating lpeak")
    N = len(series.iterations)
    lpeaks = np.zeros(N)
    ts = np.zeros(N)
    if tqdm_available:
        pbar = tqdm(total=N, bar_format=bar_format)
    lps = 0
    Nmax = N
    Nmin = 0
    lps = None
    for n, n_i in enumerate(series.iterations):
        if n_i < startit:
            Nmin = n+1
        else:
            ts[n] = series.iterations[n_i].get_attribute("time") * series.iterations[n_i].get_attribute("timeUnitSI")
            if (ts[n] >= maxtime) and (Nmax == N):
                Nmax = n
                if apply_maxtime_to_measurement:
                    print("maxtime reached")
                    break
            if lps is None:
                lps = show_lpeak(series=series, iteration=n_i, forward=forward, method="stat", p=False)
            else:
                lpeaks[n] = show_lpeak(series=series, iteration=n_i, forward=forward, method="stat", p=False) - lps
        if tqdm_available:
            pbar.update(1)
        else:
            print(n+1, "out of", N, "complete")

    if tqdm_available:
        pbar.close()
    if close:
        series.close()
    
    print("plotting lpeak")
    fig = plt.figure()
    ax = fig.add_subplot()
    if comp_func is not None:
        ax.plot(ts[Nmin:Nmax], comp_func(ts[Nmin:Nmax])-comp_func(ts[Nmin]), label=comp_func.__doc__)
    if apply_maxtime_to_measurement:
        ax.plot(ts[Nmin:Nmax], lpeaks[Nmin:Nmax], ".", label="measured lpeak")
    else:
        ax.plot(ts[Nmin:], lpeaks[Nmin:], ".", label="measured lpeak")
    ax.set_xlabel("$t$/s")
    ax.set_ylabel("$l_{peak}$/m")
    ax.set_title(title)
    if comp_func is not None:
        ax.legend()

    if ret_ax:
        if ret:
            return fig, ax, ts, ws
        else:
            return fig, ax
            
    fig.show()
    plt.show()
    
    if ret:
        return ts, lpeaks

_w = show_w
_lpeak = show_lpeak

def _prepcut(N, center_N, shape, spacing, direction, show_SI):
    """little helper function for _showit"""
    if show_SI:
        extent = [-shape/2*spacing, shape/2*spacing]
    else:
        extent = [0, shape*spacing]
    Nlo = 0
    Nhi = shape
    if N is not None:
        if center_N is not None:
            Nlo = center_N - N//2
            Nhi = center_N + N//2
            assert Nlo < Nhi
            assert Nlo >= 0, "N"+direction+" does not fit at this center."
            assert Nhi <= shape, "N"+direction+" does not fit at this center."
            if show_SI:
                center = spacing * (center_N - shape/2)
                extent[0] = center - spacing * (N//2)
                extent[1] = center + spacing * (N//2)
            else:
                extent[0] = Nlo
                extent[1] = Nhi
        else:
            assert N < shape, "N"+direction+" does not fit."
            f = 1 - (N / shape)
            Nlo = int(f/2*shape)
            Nhi = shape - int(f/2*shape)
            if show_SI:
                extent[0] *= 1-f
                extent[1] *= 1-f
            else:
                extent[0] = Nlo
                extent[1] = Nhi
    return Nlo, Nhi, extent
    

def _showit(n_i, series, direction="x", cutdim="x", cutfrac=0.5, Nx=None, Nz=None, Ny=None, center_Nx=None,
            center_Ny=None, center_Nz=None, forward="z", linthresh_frac=0.01, 
            title="", titleit=False, titleit_mul=1, titleit_after="", show_SI=True, figsize=(16, 12), mesh="E",
            show_w=False, show_lpeak=False):
    """internal function for showing a single iteration"""

    print("iteration:", n_i)
    
    assert direction in ["x", "y", "z"]
    assert cutdim in ["x", "y", "z"]
    # figuring out, which direction is neither forward nor cut
    trans = ["x", "y", "z"]
    trans.remove(cutdim)
    if forward[-1:] != cutdim:
        trans.remove(forward[-1:])
    else:
        print("both directions transversal")

    spacing = [s * series.iterations[n_i].meshes[mesh].get_attribute("gridUnitSI")
                for s in series.iterations[n_i].meshes[mesh].get_attribute("gridSpacing")]
    l = dict(zip(series.iterations[n_i].meshes[mesh].get_attribute("axisLabels"), [0,1,2]))
    if len(series.iterations[n_i].meshes[mesh]) == 3:
        shape = series.iterations[n_i].meshes[mesh][direction].shape
    elif len(series.iterations[n_i].meshes[mesh]) == 1:
        shape = series.iterations[n_i].meshes[mesh].shape
    extent = [[-sh/2*sp, sh/2*sp] for sh,sp in zip(shape, spacing)]
    Nlo = [0,0,0]
    Nhi = shape.copy()
    Nd = [Nx, Ny, Nz]
    center_Nd = [center_Nx, center_Ny, center_Nz]
    for nd, d in enumerate("xyz"):
        Nlo[l[d]], Nhi[l[d]], extent[l[d]] = _prepcut(Nd[nd], center_Nd[nd], shape[l[d]], spacing[l[d]], d, show_SI)
    Nlo[l[cutdim]] = int(shape[l[cutdim]] * cutfrac)
    Nhi[l[cutdim]] = int(shape[l[cutdim]] * cutfrac) + 1
    extent.pop(l[cutdim])

    alabel = mesh
    if len(series.iterations[n_i].meshes[mesh]) == 3:
        alabel += "_" + direction
        field = series.iterations[n_i].meshes[mesh][direction][Nlo[0]:Nhi[0],Nlo[1]:Nhi[1],Nlo[2]:Nhi[2]]
    elif len(series.iterations[n_i].meshes[mesh]) == 1:
        field = series.iterations[n_i].meshes[mesh][Nlo[0]:Nhi[0],Nlo[1]:Nhi[1],Nlo[2]:Nhi[2]]
    series.flush()
    
    if show_SI:
        if len(series.iterations[n_i].meshes[mesh]) == 3:
            field *= series.iterations[n_i].meshes[mesh][direction].get_attribute("unitSI")
        elif len(series.iterations[n_i].meshes[mesh]) == 1:
            field *= series.iterations[n_i].meshes[mesh].get_attribute("unitSI")
        alabel += "/(V/m)"
        
    print("absolute maximum:",np.max(np.abs(field)))
    if show_w or show_lpeak:
        density, sp = _get_density(series=series, iteration=n_i)
        if show_w:
            _w(density=density, spacing=sp, forward=l[forward[-1:]])
        if show_lpeak:
            _lpeak(density=density, spacing=sp, forward=l[forward[-1:]])

    assert len(trans) in [1,2], "the dimensions dont add up here"
    if len(trans) == 1:
        xlabel = forward[-1:]
        ylabel = trans[0]
        transpose = (l[trans[0]] > l[forward[-1:]])
    else:
        if l[trans[0]] < l[trans[1]]:
            xlabel = trans[0]
            ylabel = trans[1]
        else:
            xlabel = trans[1]
            ylabel = trans[0]
        transpose = False
    if show_SI:
        xlabel += "/m"
        ylabel += "/m"
    if titleit:
        title += str(n_i * titleit_mul) + titleit_after
    if "-" in forward:
        flipx = True
    else:
        flipx = False
    
    show(np.squeeze(field), extent, cutdim=l[cutdim], transpose=transpose, linthresh_frac=linthresh_frac,
            title=title, xlabel=xlabel, ylabel=ylabel, alabel=alabel, flipx=flipx, figsize=figsize)

def show_metadata(filename=None, iteration=None, series=None, s_ser=True, s_it=True, s_mesh=True, s_meshes=True,
                  s_rec=True, s_rec_c=True, s_part=True, s_parts=True, h_meshes=[], h_particles=[]):
    """Print the metadata of the openPMD-file at filename or the openPMD-series series.

    Parameters:
    filename : str
        The complete path and name of the file to be shown. Only neccessary, if no series is given.
        If series is given this is ignored.

    series : openpmd-api series object
        The series object to show the metadata of. Only neccessary, if no filename is given.

    iteration : int (optional)
        Show a specific iteration. If not set the first iteration in the file is shown.

    s_[name] : bool (very optional)
        For turning off printing of specific parts of the metadata. Set to False to turn off.
        s_ser: Series
        s_it: Iteration
        s_mesh: Meshes object
        s_meshes: Meshes
        s_rec: Records
        s_rec_c: Record components
        s_part: Particles object
        s_parts: Particles

    h_meshes : list of str (optional)
        Hides the meshes with names in the given list.

    h_particles : list of str (optional)
        Hides the particles with names in the given list.
    """
    close = False
    if series is None:
        assert filename is not None, "filename or series must be given."
        close = True
        series = io.Series(filename, io.Access.read_only)
    if s_ser:
        print("The series object",series)
        print("Attributes:")
        for attr in series.attributes:
            print(attr, series.get_attribute(attr))
    if iteration is None:
        iteration = min(series.iterations)
    i = series.iterations[iteration]
    print("\nshowing iteration number", iteration)
    if s_it:
        print("\nThe iteration object",i)
        print("Attributes:")
        for attr in i.attributes:
            print(attr, i.get_attribute(attr))

    if s_mesh:
        print("\nThe meshes:",i.meshes)
        if len(i.meshes.attributes)>0:
            print("Attributes:")
            for attr in i.meshes.attributes:
                print(attr, i.meshes.get_attribute(attr))
    if len(i.meshes)>0 and s_meshes:
        print("\nMeshes:")
        for mesh in i.meshes:
            if not mesh in h_meshes:
                print("\nmesh",mesh)
                print(i.meshes[mesh])
                print("Attributes:")
                for attr in i.meshes[mesh].attributes:
                    print(attr, i.meshes[mesh].get_attribute(attr))
                if s_rec_c:
                    if len(i.meshes[mesh]) > 1:
                        print("\nRecord components:")
                        for l in i.meshes[mesh]:
                            print("\nrecord component", l)
                            print(i.meshes[mesh][l])
                            print("Attributes:")
                            for attr in i.meshes[mesh][l].attributes:
                                print(attr, i.meshes[mesh][l].get_attribute(attr))
                    else:
                        print("\nDtype", i.meshes[mesh].dtype)
                        print("shape", i.meshes[mesh].shape)

    if s_part:
        print("\nThe Particles:",i.particles)
        if len(i.particles.attributes)>0:
            print("Attributes:")
            for attr in i.particles.attributes:
                print(attr, i.particles.get_attribute(attr))
    if len(i.particles)>0 and s_parts:
        print("\nParticles:")
        for particle in i.particles:
            if not particle in h_particles:
                print("\nparticle", particle)
                print(i.particles[particle])
                print("Attributes:")
                for attr in i.particles[particle].attributes:
                    print(attr, i.particles[particle].get_attribute(attr))
                if s_rec:
                    print("\nRecords:")
                    for rec in i.particles[particle]:
                        print("\nRecord", rec)
                        print(i.particles[particle][rec])
                        print("Attributes:")
                        for attr in i.particles[particle][rec].attributes:
                            print(attr, i.particles[particle][rec].get_attribute(attr))
                        if s_rec_c:
                            if len(i.particles[particle][rec]) > 1:
                                print("\nRecord components:")
                                for l in i.particles[particle][rec]:
                                    print("\nRecord component", l)
                                    print(i.particles[particle][rec][l])
                                    print("Attributes:")
                                    for attr in i.particles[particle][rec][l].attributes:
                                        print(attr, i.particles[particle][rec][l].get_attribute(attr))
                            else:
                                print("\nDtype", i.particles[particle][rec].dtype)
                                print("shape", i.particles[particle][rec].shape)
    if close:
        series.close()

def _lwfa_data(series, iteration, n, xmax):
    """helper function for show_lwfa"""
    it = series.iterations[iteration]

    # Electric field
    shape = it.meshes["E"]["y"].shape
    Ey = it.meshes["E"]["y"][shape[0]//2,:,shape[2]//2]
    Ex = it.meshes["E"]["x"][:,:,shape[2]//2]
    series.flush()
    Ey *= it.meshes["E"]["y"].get_attribute("unitSI")
    Ex *= it.meshes["E"]["x"].get_attribute("unitSI")

    # electron positions
    z = it.particles["e"]["positionOffset"]["z"][:]
    y = it.particles["e"]["positionOffset"]["y"][:]
    x = it.particles["e"]["positionOffset"]["x"][:]
    series.flush()
    z = z * it.particles["e"]["positionOffset"]["z"].get_attribute("unitSI")
    y = y * it.particles["e"]["positionOffset"]["y"].get_attribute("unitSI")
    x = x * it.particles["e"]["positionOffset"]["x"].get_attribute("unitSI")

    # correcting electron positions
    dz = it.particles["e"]["position"]["z"][:]
    dy = it.particles["e"]["position"]["y"][:]
    dx = it.particles["e"]["position"]["x"][:]
    series.flush()
    z += dz * it.particles["e"]["position"]["z"].get_attribute("unitSI")
    y += dy * it.particles["e"]["position"]["y"].get_attribute("unitSI")
    x += dx * it.particles["e"]["position"]["x"].get_attribute("unitSI")
    
    z -= (np.max(z) + np.min(z))/2
    y -= (np.max(y) + np.min(y))/2

    # electron momentum
    px = it.particles["e"]["momentum"]["x"][:]
    py = it.particles["e"]["momentum"]["y"][:]
    pz = it.particles["e"]["momentum"]["z"][:]
    weight = it.particles["e"]["weighting"][:]
    series.flush()
    px = px * it.particles["e"]["momentum"]["x"].get_attribute("unitSI") / weight
    py = py * it.particles["e"]["momentum"]["y"].get_attribute("unitSI") / weight
    pz = pz * it.particles["e"]["momentum"]["z"].get_attribute("unitSI") / weight

    # calculating the electron energy
    E = np.sqrt((px**2+py**2+pz**2)*c**2 + m_e**2*c**4) / e

    # cutting out the electrons far from the center slice
    midx = (np.max(x) + np.min(x)) / 2
    outx = (np.max(x) - np.min(x)) / 2
    maxx = midx + outx * xmax
    minx = midx - outx * xmax
    
    return (Ey, Ex, np.extract(np.logical_and(x[::n] < maxx, x[::n] > minx), z[::n]), 
            np.extract(np.logical_and(x[::n] < maxx, x[::n] > minx), y[::n]), np.extract(np.logical_and(x[::n] < maxx, x[::n] > minx), E[::n]))

def show_lwfa(filename=None, series=None, iteration=None, n=1, xmax=0.01, s=0.1, xbound=(None, None), ybound=(None, None), title="", figsize=(20, 12), ret_ax=False):
    """Generates an LWFA image from the file at filename and the iteration iteration similar to the Bachelor's thesis.
    Only works on openPMD files generated by PIConGPU.
    
    Parameters:
    filename : str
        The full path at which the file can be found.
        Either filename or series must be given.
        
    series : openpmd_api Series object
        The series from which this should be plotted.
        Either filename or series must be given.

    iteration : int (optional)
        The iteration to use. If None it uses the first iteration inn the file.

    n : int (optional)
        only uses every nth electron

    xmax : float in [0,1] (optional)
        plot electrons only this far from the center

    s : float (optional)
        The size of the drawn electrons

    xbound : tuple of float (optional)
        Sets the bounds in the forward direction to which the plot will be drawn. x refers to the plot axis.

    ybound : tuple of float (optional)
        Sets the bounds in the transverse direction to which the plot will be drawn. y refers to the plot axis.

    title : str (optional)
        The title of the plot.

    figsize : tuple of float (optional)
        Sets the size of the matplotlib figure.

    ret_ax : bool (optional)
        Whether to return the axes objects instead of drawing them to the screen.

    Returns:
    fig : pyplot figure object (optional)
        only returned if ret_ax is True.

    ax : pyplot axes object (optional)
        only returned if ret_ax is True.
    """
    print("fetching data")
    close = False
    if series is None:
        assert filename is not None, "filename or series must be given."
        close = True
        series = io.Series(filename, io.Access.read_only)
    if iteration is None:
        iteration = min(series.iterations)
    it = series.iterations[iteration]

    # gathering data
    Eshape = it.meshes["E"]["y"].shape
    Espacing = [s * it.meshes["E"].get_attribute("gridUnitSI")
                for s in it.meshes["E"].get_attribute("gridSpacing")]
    Ey, Ex, z, y, E = _lwfa_data(series, iteration, n, xmax)
    mz = (np.arange(Eshape[0]+1) - Eshape[0]/2) * Espacing[0]
    my = (np.arange(Eshape[1]+1) - Eshape[1]/2) * Espacing[1]
    mEx = np.max(np.abs(Ex))
    mEy = np.max(np.abs(Ey))

    if close:
        series.close()
    
    # making the figure
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot()
    # first the laser
    pcolormesh = ax.pcolormesh(my, mz, Ex, label="laser pulse", cmap="seismic", vmax=mEx, vmin=-mEx)
    # then the electrons
    scatter = ax.scatter(y, z, s=0.1, c=E, label="electrons", cmap="inferno", alpha=0.7, norm=clr.LogNorm())
    # then the accelerating field on a different y-axis
    ax2 = ax.twinx()
    ax2.plot(my[:-1], 0*my[:-1], c="black", alpha=0.5)
    ax2.plot(my[:-1], Ey, label="accelerating field")
    # then all the labels and stuff
    ax.set_xlabel("$y$/m")
    ax.set_ylabel("$z$/m")
    ax.set_title(title)
    ax2.set_ylabel("accelerating field $E_y$/(V/m)")
    plt.colorbar(pcolormesh, ax=ax, label="laser field $E_x$/(V/m)")
    plt.colorbar(scatter, ax=ax, label="electron energy E/eV")
    ax.set_xbound(xbound)
    ax.set_ybound(ybound)
    if ret_ax:
        return fig, ax
    print("drawing")
    plt.show()

def _show_lwfa(n_i, series, title="", titleit=False, titleit_mul=1, titleit_after="",  **kwargs):
    """helper function for the show_lwfa option of the show_file function."""
    print("Iteration", n_i)
    if titleit:
        title += str(n_i * titleit_mul) + titleit_after
    return show_lwfa(series=series, iteration=n_i, title=title, **kwargs)

def show_file(filename, itstep=1, metadata=False, iteration=None, show_time=False, show_lwfa=False, **kwargs):
    """shows the E-field in the iterations of the openPMD-file.

    Parameters:
    filename : str
        Complete path to the file including the name

    itstep : int (optional)
        Only show every itstep iteration. Make sure the iterations you want to see are a multiple of this value.

    metadata : bool (optional)
        Show metadata before showing the contents.

    iteration : int or list of int (optional)
        Only show specific iteration(s). ignores itstep.

    show_time : bool (optional)
        prints the runtime when an iteration is done.

    show_lwfa : bool (optional)
        If false uses the show function to show the electric field on a symlog plot.
        If True uses the show_lwfa function to display a plot showing LWFA.
        In that case it is assumed, that the file is a standard PIConGPU output.
        

    Other Parameters:
    title : str (optional)
        The title for the plots

    titleit : bool (optional)
        Wether the iteration number should be added to the title of each plot.

    titleit_mul : float or int (optional)
        If titleit is set to True it will be multiplied by this number.

    titleit_after : str (optional)
        If titleit is set to True this text will be after the iteration number in the title.

    figsize : tuple of float (optional)
        Sets the size of the matplotlib figure.
        

    Parameters specific to show_lwfa=False:
    direction : "x", "y" or "z" (optional)
        The field direction which is to be shown

    cutdim : "x", "y" or "z" (optional)
        Defines the dimension of the array going out of the screen in the plot.

    cutfrac : float [0,1] (optional)
        Defines, where along the cutdim the plot should be in terms of fraction along it.

    Nx, Ny, Nz : int (optional)
        Only show the central N points in each direction respectively

    center_Nx, center_Nx, center_Nx : int (optional)
        Only relevant, if the respective N is given. This is the center point to be shown.

    linthresh_frac : float (optional)
        Fraction of the absolute maximum from where the normalistion should be linear.

    show_SI : bool (optional)
        If True the axes are written in SI-units, otherwise they are written in points.

    show_w : bool (optional)
        wether to calculate and print the beam waist at this point.

    show_lpeak : bool (optional)
        wether to calculate and print the peak position on the longitudinal axis.


    Parameters specific to show_lwfa=True:
    n : int (optional)
        only show every nth electron.
        
    xmax : float in [0, 1] (optional)
        Only plot electron this far from the center. 
        This is a fraction of the longest distance of an electron and the center. 
        
    s : float (optional)
        the size of the drawn electrons.
        
    xbound, ybound : tuple of float (optional)
        Draw the plot only wwithin these bounds. In SI units.
    """
    if show_time:
        if "time showing file" in get_clocks():
            set_clock("time showing file", 0)
        else:
            start_clock("time showing file")
    series = io.Series(filename, io.Access.read_only)
    
    if metadata:
        show_metadata(iteration=iteration, series=series)
        print("\nData:")

    func = _showit
    if show_lwfa:
        func = _show_lwfa
    
    if iteration is not None:
        try:
            n = iteration[0]

            mul = True
        except (TypeError, IndexError):
            mul = False
        
        if mul:
            # show the specified iterations
            for n_i in iterations:
                if n_i in series.iterations:
                    func(n_i, series, **kwargs)
                    if show_time:
                        print_clock("time showing file")
                else:
                    print("Warning: Iteration", n_i, "does not exist.")
        else:
            # show just one iteration
            assert iteration in series.iterations, "Specified iteration does not exist."
            func(iteration, series, **kwargs)
            if show_time:
                print_clock("time showing file")
    else:
        # show multiple iterations
        for n in range(int(max(series.iterations)/itstep)+1):
            if itstep * n in series.iterations:
                func(itstep * n, series, **kwargs)
                if show_time:
                    print_clock("time showing file")
    series.close()
    plt.show()

