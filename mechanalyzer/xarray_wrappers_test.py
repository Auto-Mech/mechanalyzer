import xarray_wrappers
import xarray
import numpy

Temps = [1000, 1500, 2000, 2500]
Press = [1, 10, numpy.inf]
Rates = [[1e1, 1e2, 1e3, 1e4], [1e5, 1e6, 1e7, 1e8], [1e9, 1e10, 1e11, 1e12]]

Ktp = xarray_wrappers.from_data(Temps, Press, Rates)
print(Ktp)

def test_set_rates():
    ktp = xarray_wrappers.set_rates(Ktp, Rates)
    print(ktp)

def test_get_temperatures():
    temp = xarray_wrappers.get_temperatures(Ktp)
    print(temp)


def test_get_pressures():
    pres = xarray_wrappers.get_pressures(Ktp)
    print(pres)


def test_get_values():
    vals = xarray_wrappers.get_values(Rates)
    print(vals)


def test_get_pslice():
    pslice = xarray_wrappers.get_pslice(Ktp, numpy.inf)
    print(pslice)


def test_get_tslice():
    tslice = xarray_wrappers.get_tslice(Ktp, 1500)
    print(tslice)


def test_get_spec_vals():
    vals = xarray_wrappers.get_spec_vals(Ktp, 1500, 1)
    print(vals)


def test_get_ipslice():
    ipslice = xarray_wrappers.get_ipslice(Ktp, 0)
    print(ipslice)


def test_get_itslice():
    itslice = xarray_wrappers.get_itslice(Ktp, 0)
    print(itslice)
    
test_set_rates()
test_get_pressures()
test_get_temperatures()
test_get_values()
test_get_pslice()
test_get_tslice()
test_get_spec_vals()
test_get_ipslice()
test_get_itslice()
