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

    return (list(ktp["temp"]))


def get_pressures(ktp):
    # WORK IN PROGRESS
    # Gets the pressure values

    return (list(ktp["pres"]))


def get_values(ktp_vals):
    # WORK IN PROGRESS
    # Gets the KTP values

    return (ktp_vals)

def get_pslice(ktp):
    # WORK IN PROGRESS
    # Get a slice at a selected pressure value

    return (ktp.sel(pres = numpy.inf))

def get_tslice(ktp, t):
    # WORK IN PROGRESS
    # Get a slice at a selected temperature value

    return(ktp.sel(temp = t))

def get_spec_vals(ktp, t, p):
    # WORK IN PROGRESS
    # Get a specific value at a selected temperature and pressure value

    return (ktp.sel(temp = t, pres = p))

def get_ipslice(ktp, i):
    # WORK IN PROGRESS
    # Get a slice at a selected pressure index

    return(ktp.isel(pres = i))

def get_itslice(ktp, i):
    # WORK IN PROGRESS
    # Get a slice at a selected temperature index

    return (ktp.isel(temp = i))

def get_dims(ktp):
    # WORK IN PROGRESS
    # Get all dimensions associated with an array

    return(ktp.dims)

def get_t_coords(ktp):
    # WORK IN PROGRESS
    # Get co-ordinates associated with a dimension for temperature
    
    return(ktp.temp.data)

def get_p_coords(ktp):
    # WORK IN PROGRESS
    # Get co-ordinates associated with a dimension for pressure

    return(ktp.pres.data)


# Setters
def set_values(ktp_vals):
    # WORK IN PROGRESSS
    # Sets the KTP values
