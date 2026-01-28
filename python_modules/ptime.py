import time

_starts = {"__rtime__": time.time()}

def _printf(string, filename):
    file = open(filename, "r")
    lines = file.readlines()
    file.close()
    lines.append(string+"\n")
    file = open(filename, "w")
    file.writelines(lines)
    file.close()

def _format_time(t):
    """format the given time"""
    s = t % 60
    t_str = str(s)[:6] + " s"
    if t > 60:
        m = (int(t) // 60) % 60
        t_str = str(m) + " min " + t_str
    if t > 3600:
        h = (int(t) // 3600) % 24
        t_str = str(h) + " h " + t_str
    if t > 86400:
        d = int(t) // 86400
        t_str = str(d) + "d " + t_str
    return t_str

def rtime():
    """returns the time since this module was imported.
    """
    return time.time() - _starts["__rtime__"]

def ptime(title="runtime", filename=None):
    """prints the time since this module was imported with nice formatting.
    
    Parameters:
    title : str (optional)
        the text that should be printed before the actual time.

    filename : str (optional)
        If this is not None, the time will be appended to the specified fiel instead of printed.
    """
    rt = rtime()
    s = title+": "+_format_time(rt)
    if filename is None:
        print(s)
    else:
        _printf(s, filename)

def start_clock(name, value=0, start_paused=False):
    """start a clock that can later be read out. No process is running in the background for this.

    Parameters:
    name : str
        The name of the clock.

    value : float (optional)
        The value the clock should be set to.
    """
    assert name != "__rtime__", "clock name '__rtime__' cannot be used."
    if name in _starts:
        raise ValueError(f"clock '{name}' already exists")
    if start_paused:
        _starts[name] = (-1, value)
    else:
        _starts[name] = (time.time(), value)

def read_clock(name):
    """returns the time on the clock.

    Parameters:
    name : str
        The name of the clock that should be read.
    """
    assert name in _starts, f"clock {name} does not exist"
    assert name != "__rtime__", "clock name '__rtime__' cannot be used."
    if _starts[name][0] > 0:
        return time.time() - _starts[name][0] + _starts[name][1]
    else:
        return _starts[name][1]

def set_clock(name, value):
    """Sets the time on the clock to the specified value. If the clock is paused it stays that way.

    Parameters:
    name : str
        The name of the clock that should be set.

    value : float
        The value the clock should be set to.
    """
    assert name in _starts, f"clock {name} does not exist"
    assert name != "__rtime__", "clock name '__rtime__' cannot be used."
    if _starts[name][0] > 0:
        _starts[name] = (time.time(), value)
    else:
        _starts[name] = (-1, value)

def pause_clock(name):
    """pauses the specified clock. If it is already paused, nothing changes.

    Parameters:
    name : str
        The name of the clock that should be paused.
    """
    assert name in _starts, f"clock {name} does not exist"
    assert name != "__rtime__", "clock name '__rtime__' cannot be used."
    _starts[name] = (-1, read_clock(name))

def unpause_clock(name):
    """unpauses the specified clock. If it is already unpaused, nothing changes.

    Parameters:
    name : str
        The name of the clock that should be unpaused.
    """
    assert name in _starts, f"clock {name} does not exist"
    assert name != "__rtime__", "clock name '__rtime__' cannot be used."
    _starts[name] = (time.time(), read_clock(name))

def print_clock(name, filename=None):
    """prints the time on the clock.

    Parameters:
    name : str
        The name of the clock that should be read.

    filename : str (optional)
        If this is not None, the time will be appended to the specified fiel instead of printed.
    """
    assert name in _starts, f"clock {name} does not exist"
    assert name != "__rtime__", "clock name '__rtime__' cannot be used."
    s = name+": "+_format_time(read_clock(name))
    if filename is None:
        print(s)
    else:
        _printf(s, filename)

def get_clocks():
    """returns a list of all active clocks
    """
    l = list(_starts.keys())
    l.remove("__rtime__")
    return l