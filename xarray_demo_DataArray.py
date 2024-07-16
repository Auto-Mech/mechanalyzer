import xarray
import numpy

temps = [1000, 1500, 2000, 2500]
press = [1, 10, numpy.inf]

# Define the array of rates (note: it makes sense to structure it as a list 
# of rates at a given P based on how the MESS parser is written, though it
# doesn't really matter) 
ktp_vals = [[1e1, 1e2, 1e3, 1e4], [1e5, 1e6, 1e7, 1e8], [1e9, 1e10, 1e11, 1e12]]

# Create the DataArray
ktp = xarray.DataArray(ktp_vals, [("pres", press), ("temp", temps)])
print(ktp)

# Get a slice at a selected pressure or temperature value
print('\nPressure slice by value')
print(ktp.sel(pres=numpy.inf))  # the fact that you can use np.inf is very convenient!
print('\nTemperature slice by value')
print(ktp.sel(temp=1500))

# Get a specific value at a selected temperature value and pressure value
print('\nTemperature slice by value')
print(ktp.sel(temp=1500, pres=1))

# Get a slice at a selected pressure or temperature index
print('\nPressure slice by index')
print(ktp.isel(pres=0))
print('\nTemperature slice by index')
print(ktp.isel(temp=0))

# Can also add some metadata (probably unnecessary for ktp objects, but cool)
ktp.attrs["units"] = "s^-1"  # any number of arbitrary attributes can be added
ktp.name = "N2O=N2+O"

# Can do various arithmetic operations
print('\nVarious operations')
print(ktp + ktp)
print(3 * ktp)

# To get the actual values, use the values OR data attribute
print('\nGet the values, which should be a numpy array')
print(ktp.values)
print('type of ktp.values: ', type(ktp.values))
print('\nTemperature slice by value, as an array')
print(ktp.sel(temp=1500, pres=1).data)
print(type(ktp.sel(temp=1500, pres=1).data))

# Can get all dimensions associated with an array usings the dims attribute
print('\nNames of all dimensions, as a tuple')
print(ktp.dims)

# Can get the coordinates associated with a dimension by using that coordinate name
print('\nPressure values: ', ktp.pres.data)
print('\nTemperature values: ', ktp.temp.data)

# Can do lookup on values with the nearest value
# This is giving an error even though I thought I copied it from the website?
#print('\nTemperature slice by nearest lookup')
#print(ktp.sel(temp=1600), method='nearest')
