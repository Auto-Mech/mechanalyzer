"""
Wrappers for the new xarray system. Constructors, Getters, then Setters.
"""

import xarray

# Constructors
def from_data(temps, press, rates):
    """
    Construct a KTP DataArray from data
    """

    ktp = xarray.DataArray(rates, [("pres", press), ("temp", temps)])

    return ktp



# Getters
def get_pressures(ktp):
    """
    Gets the temperature values
    """

    return ktp.pres.data


def get_temperatures(ktp):
    """
    Gets the pressure values
    """

    return ktp.temp.data


def get_values(rates):
    """
    WORK IN PROGRESS
    Gets the KTP values
    """

    return rates


def get_pslice(ktp, ip):
    """
    Get a slice at a selected pressure value
    """

    return ktp.sel(pres=ip)


def get_tslice(ktp, it):
    """
    Get a slice at a selected temperature value
    """

    return ktp.sel(temp=it)


def get_spec_vals(ktp, it, ip):
    """
    Get a specific value at a selected temperature and pressure value
    """

    return ktp.sel(temp=it, pres=ip)


def get_ipslice(ktp, ip):
    """
    Get a slice at a selected pressure index
    """

    return ktp.isel(pres=ip)


def get_itslice(ktp, it):
    """
    Get a slice at a selected temperature index
    """

    return ktp.isel(temp=it)



# Setters
def set_rates(ktp, rates):
    """
    DOES NOT WORK YET. Still fixing!
    Sets the KTP values
    """

    #ktp = ktp.loc[rates]
    #return ktp
    pass
