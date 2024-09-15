import xarray
import numpy

# Constructors
def from_data(temps, press, ktp_vals):
    # WORK IN PROGRESS - SHOULD NOT WORK YET JUST TESTING THINGS OUT
    # Construct a KTP DataArray from data
    ktp = xarray.DataArray(ktp_vals, [("pres", press), ("temp", temps)])

    return ktp

# Getters
def get_temperatures(ktp):
    # WORK IN PROGRESS
    # Gets the temperature values

    temps = ktp.get("temps")

    return(temps)


def get_pressures(ktp):
    # WORK IN PROGRESS
    # Gets the pressure values

    press = ktp["press"]

    return(press)


def get_values(ktp):
    # WORK IN PROGRESS
    # Gets the KTP values

    vals = ktp(ktp_vals)

    return(vals)

def get_pslice(ktp):
    # WORK IN PROGRESS
    # Get a slice at a selected pressure value

    pslice = ktp.sel(pres = numpy.inf)

    return pslice

def get_tslice(ktp, t):
    # WORK IN PROGRESS
    # Get a slice at a selected temperature value

    tslice = ktp.sel(temp = t)

    return tslice

def get_spec_vals(ktp, t, p):
    # WORK IN PROGRESS
    # Get a specific value at a selected temperature and pressure value

    spec_vals = ktp.sel(temp = t, pres = p)

    return spec_vals

def get_ipslice(ktp, i):
    # WORK IN PROGRESS
    # Get a slice at a selected pressure index

    ipslice = ktp.isel(pres = i)

    return ipslice

def get_itslice(ktp, i):
    # WORK IN PROGRESS
    # Get a slice at a selected temperature index

    itslice = ktp.isel(temp = i)

    return itslice


# Setters
def set_values(ktp_obj, ktp_vals):
    # WORK IN PROGRESSS
    # Sets the KTP values

    return(ktp_vals)
