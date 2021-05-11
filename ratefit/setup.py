""" Install ratefit for dealing with calculations and fits
    to common rate-constant functional form expressions
"""

from distutils.core import setup


setup(
    name="ratefit",
    version="0.1.0",
    packages=[
        'ratefit',
        'ratefit.calc',
        'ratefit.fit',
        'ratefit.fit.arrhenius',
        'ratefit.fit.chebyshev',
        'ratefit.fit.troe'],
    package_dir={
        'ratefit': 'ratefit',
        'calc': 'calc',
        'fit': 'fit'
    },
    package_data={
        'ratefit': ['tests/data/*',
                    'fit/arrhenius/dsarrfit.mako',
                    'fit/troe/troefit.mako']
    }
)
