"""
  Handle kTP rate objects

ktps:
    # TP DCT
    {P1: (T1, T2, T3, ...),
     P2: (T1, T2, T3, ...),
     P3: (T1, T2, T3, ...)},
    # kTP DCT
    {P1: (kTP1, kTP2, kTP3, ...)
     P2: (kTP1, kTP2, kTP3, ...)
     P3: (kTP1, kTP2, kTP3, ...)}

params
    # TP DCT
    {P1: (T1, T2, T3, ...),
     P2: (T1, T2, T3, ...),
     P3: (T1, T2, T3, ...)},
    # PARAM DCT
    {P1: params
     P2: params
     P3: params}
    # FIT DCT
    {P1: fit
     P2: fit
     P3: fit}

err
    # TP DCT
    {P1: (T1, T2, T3, ...),
     P2: (T1, T2, T3, ...),
     P3: (T1, T2, T3, ...)},
    # PARAM DCT
    {P1: (max_err, mean_err)
     P2: (max_err, mean_err)
     P3: (max_err, mean_err)}

"""

# Constructor
def from_data(temps, pressures, ktps):
    """ Build
    """


def from_dcts(tp_dct, ktp_dct):
    """ Build
    """
    return (tp_dct, ktp_dct)


# getters
def tp_dct(ratek_obj):
    """ Obtain dct listing out what temps have calc'd rate constants
        at each pressure
    """
    dct, _ = ratek_obj
    return dct


def ktp_dct(ratek_obj):
    """ Obtain dct listing out what temps have calc'd rate constants
        at each pressure
    """
    _, dct = ratek_obj
    return dct


def pressures(ktp):
    """ grab the pressures
    """
    dct = tp_dct(ktp)
    return tuple(dct.keys())


def temperatures(ktp, pressure):
    """ grab the temps
    """
    dct = tp_dct(ktp)
    return dct[pressure]


def kts(ktp, pressure):
    """ grab the temps
    """
    dct = ktp_dct(ktp)
    return dct[pressure]


# converter
def ktp_matrix(ktp, inpressures=(), intemps=()):
    """ Build numpy matrix Obtain rate constants from structure
    """

    # Get the pressures
    if inpressures:
        assert set(inpressures) <= set(pressures(ktp)), (
            'Requested pressures for which no rateks exist')
        rateps = inpressures
    else:
        rateps = pressures(ktp)

    # Get the rate constant at each pressure
    for pressure in pressures:
        # Get the temperatures
        if intemps:
            # need assert
            alltemps = temperatures(ktp, pressure)
            idxs = (alltemps.index(temp) for temp in intemps)
        else:
            kts = kts(ktp, pressure)


# formatter
def _highp_str():
    if val == 1e20:
        out = 'High'
    else:
        out = val
    return out
