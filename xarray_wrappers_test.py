import xarray_wrappers
import xarray
import numpy

temps = [1000, 1500, 2000, 2500]
press = [1, 10, numpy.inf]
ktp_vals = [[1e1, 1e2, 1e3, 1e4], [1e5, 1e6, 1e7, 1e8], [1e9, 1e10, 1e11, 1e12]]

ktp = xarray_wrappers.from_data(temps, press, ktp_vals)
print(ktp)

def test_get_temperatures():
    temp = xarray_wrappers.get_temperatures(ktp)
    print(temp)

def test_get_pressures():
    pres = xarray_wrappers.get_pressures(ktp)
    print(pres)

def test_get_values():
    vals = xarray_wrappers.get_values(ktp)
    print(vals)

def test_get_pslice():
    pslice = xarray_wrappers.get_pslice(ktp)
    print(pslice)

def test_get_tslice():
    tslice = xarray_wrappers.get_tslice(ktp)
    print(tslice)

def test_get_spec_vals():
    vals = xarray_wrappers.get_spec_vals(ktp, temp = 1500, pres = 1)
    print(vals)

def test_get_ipslice():
    ipslice = xarray_wrappers.get_ipslice(ktp, i = 0)
    print(ipslice)

def test_get_itslice():
    itslice = xarray_wrappers.get_itslice(ktp, i = 0)
    print(itslice)

test_get_temperatures()
test_get_pressures()
test_get_values()
test_get_pslice()
test_get_tslice()
test_get_spec_vals()
test_get_ipslice()
test_getitslice()
