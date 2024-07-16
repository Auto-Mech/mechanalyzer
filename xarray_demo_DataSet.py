import xarray
import numpy

temps = [1000, 1500, 2000, 2500]
enthalpy = [60, 70, 80, 90]
gibbs = [10, 15, 20, 25]
c_p = [1, 1.1, 1.2, 1.3]
entropy = [3, 4, 5, 6]

therm = xarray.Dataset(
    {
        "enthalpy": (["temp"], enthalpy),
        "gibbs": (["temp"], gibbs),
        "c_p": (["temp"], c_p),
        "entropy": (["temp"], entropy),
    },
    coords={
        "temp": temps,
    },
)

print(therm)
print(therm.temp)
print(therm.temp.values)


